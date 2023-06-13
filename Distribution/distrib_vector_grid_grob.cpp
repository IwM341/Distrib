#include "functions.hpp"
#include <cmath>
#include <fstream>

#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/arg_parser.hpp"
#include "func/dens_grid.hpp"
#include "func/matrix_functions.hpp"
#include "grob/csv_io.hpp"
#include "func/move_to_go.hpp"
#include "func/load_histo.hpp"
#include "grob/grid_gen.hpp"
#include "../factors/factors.hpp"
int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}


auto PhiFactorSS = [](double q)->double{
    return 1.0;
};

Function::VectorGrid<double> CreateEGrid_sqrt(size_t N,double Emin,double eps = 1e-4){
    return  DensityGrid(Emin,0,[Emin,eps](double x){
        double t = (x-Emin)/std::abs(Emin);
        return 1/std::abs(Emin)*0.5*( 0.5/sqrt(t+eps)+0.5/sqrt(1-t+eps));
    },N);
}

auto CreateSqrtGrid(size_t N,double a,double b,double eq = 1,
                    double ineq_a = 0,double ineq_b = 0,double eps = 1e-4){
    return grob::density_grid(a,b,[=](double x){
        double t = (x-a)/(b-a);
        return eq + ineq_a*(0.5/sqrt(t+eps))+ineq_b*(0.5/sqrt(1-t+eps));
    },N,20,30);
}

template <typename Functype_E_L,typename Functype_N_E>
auto Create_EL_grid(const Function::VectorGrid<double>&Egrid,
               const Functype_E_L&rhoL,
               const Functype_N_E&N_num){
    return Function::Histogramm(Egrid,
                              Vector(Egrid.size()-1,
                                     [&Egrid,&rhoL,&N_num](size_t i){
                                       return Function::Histogramm<double,Function::VectorGrid<double>>(
                                            DensityGrid(0,1,Function::FunctorM(Egrid[i],rhoL),N_num(Egrid[i]))
                                       );
                                 }
                             ));
}


struct dF_Nuc_naive{
    size_t N,Z;
    double const_factor;
    dF_Nuc_naive(size_t Z = 1,size_t  N = 0,double const_factor = 1):Z(Z),N(N),const_factor(const_factor){}
    MC::MCResult<double> EnergyLoss(double Emax)const{
        return MC::MCResult<double>(0,1);
    }
    double ScatterFactor(double q,double enLoss)const{
        return const_factor;
    }
};



struct dF_Nuc_No_FormFactor{
    dF_Nuc_No_FormFactor(){}
    inline MC::MCResult<double> EnergyLoss(double Emax)const noexcept{
        return MC::MCResult<double>(0,1);
    }
    inline double ScatterFactor(double q,double enLoss)const noexcept{
        return 1;
    }
};

struct random_shuffler{
    size_t N;
    size_t x;
    random_shuffler( size_t N,size_t x):N(N),x(x){
        size_t MaxSize = -1;
        if(x > MaxSize/N){
            std::cout << "warning: x*N > " << MaxSize <<std::endl;
        }
        else{
            if(std::gcd(x,N) != 1){
                std::cout << "gcd(x,N) != 1" << std::endl;
            }
        }
    }
    random_shuffler(size_t N):N(N){
        x = 0.57*N;
        for(;std::gcd(x,N)!=1;++x);

    }
    size_t operator [](size_t j) const{return (j*x) % N;}
};


template <typename HistoType>
struct  grob_histo_adapter{
    HistoType & Histo;
    grob_histo_adapter(HistoType & Histo) noexcept :Histo(Histo){}
    operator HistoType &() noexcept{
        return Histo;
    }
    operator HistoType const&()const noexcept{
        return Histo;
    }

    template <typename...Args>
    inline auto putValue(Args &&...args) noexcept->decltype(Histo.put(std::forward<Args>(args)...)){
        return Histo.put(std::forward<Args>(args)...);
    }
};

template <typename HistoType,typename LE_FuncType>
struct  LE_histo_adapter{
    HistoType & Histo;
    LE_FuncType const & LE_func;
    LE_histo_adapter(HistoType & Histo,LE_FuncType const & LE_func) noexcept :
        Histo(Histo),LE_func(LE_func){}
    operator HistoType &() noexcept{
        return Histo;
    }
    operator HistoType const&()const noexcept{
        return Histo;
    }

    inline bool putValue(double value,double E,double L){
        #ifndef _NDEBUG
        double Lmax = LE_func(E);
        double l = L/Lmax;
        //auto [IsInHisto,Index] =  Histo.Grid.spos(_T(E,l));
        return Histo.put(value,E,l);
        #else
        return Histo.put(value,E,L/LE_func(E));
        #endif
    };
};

/*!
 * \brief FillScatterHisoFromElement
 * \param N_el concentration(r)  function
 * \param Therm Therm(r)  function
 * \param phi potential(r)  function
 * \param Vesc Vesc(r = 1)
 * \param VescR Vesc(r)
 * \param dF nuclei factor
 * \param mk
 * \param dmk
 * \param mp target mass
 * \param S_HL Histo of Histos scatter from L to H
 * \param S_LH Histo of Histos scatter from H to L
 * \param E_H Histo Evaporation from H
 * \param E_L Histo Evaporation from L
 * \param PhiFactor function of scatter
 * \param Nmk_HL MK samples in scatter from L to H
 * \param Nmk_LH MK samples in scatter from H to L
 */
template <typename N_type,typename Therm_type,typename Phi_Type,typename VescFuncType,
          typename dF_Type,typename EvapHisto,typename ScatterHisto,typename LE_FuncType,
          typename PhiFactorType>
inline void FillScatterHistoFromElement(const N_type & N_el,Therm_type const & Therm,Phi_Type const& phi,
                                       double Vesc, const VescFuncType & VescR,
                                       const dF_Type&dF,double mk,double dmk,double mp,
                                ScatterHisto &S_LH,ScatterHisto &S_HL,EvapHisto & E_H,EvapHisto & E_L,
                                       LE_FuncType const & LE_func,
                                const PhiFactorType & PhiFactor,size_t Nmk_HL,size_t Nmk_LH){
    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto G1 = [](){
        constexpr double fc = 1.0/(RAND_MAX+1.0);
        return fc*(rand()+fc*rand());
    };

    auto NH = S_LH.Grid.size();
    auto NL = S_HL.Grid.size();


    std::cout << "calculating scatter hl" <<std::endl;
    show_prog prog_hl(100);
    prog_hl.show(0);
    size_t total_progress_hl = NL;
    size_t curr_progress_hl = 0;

    #pragma omp parallel for
    for(size_t i=0;i<NL;++i){
        auto Index = S_HL.Grid.FromLinear(i);
        const auto & [e0e1,l0l1] = S_HL.Grid[Index];

        if(S_HL.Grid.LinearIndex(Index)!= i){
            print("S_HL.Grid.LinearIndex(Index)!= i");
        }

        double e_nd=  e0e1.center();
        double l_nd=  l0l1.center()*maxLndf(phi,e_nd);

        auto TI = CalculateTrajectory(phi,e_nd,l_nd,100);
        auto phi_min = phi(TI.Trajectory.values[0]);

        /*check*//*
        for(auto j = 0;j<TI.Trajectory.Grid.size()-1;++j){
            double r = 0.5*(TI.Trajectory.values[j+1]+TI.Trajectory.values[j]);
            double v_r = (TI.Trajectory.values[j+1]-TI.Trajectory.values[j])/
                    (TI.Trajectory.Grid[j+1]-TI.Trajectory.Grid[j]);
            double v_phi = sqrt(l_nd*l_nd/(r*r+1e-10)+ v_r*v_r);
            double E = v_phi*v_phi-phi(r);
            if(std::abs((E-e_nd)/(E+e_nd-phi(r))) > 1e-1 && l_nd>1e-10 && j>2){

                PVAR(j);
                COMPARE(E,e_nd);
                PVAR(r);
                PVAR(v_r);
                PVAR(v_phi);
                PVAR(l_nd);
                print("T[j] = ",TI.Trajectory.Grid[j]);
                print("T[j+1] = ",TI.Trajectory.Grid[j+1]);
                print("r[j] = ",TI.Trajectory.values[j]);
                print("r[j+1] = ",TI.Trajectory.values[j+1]);
                auto TI = CalculateTrajectory(phi,e_nd,l_nd,100);
                exit(0);
            }

        }*/
        if(Therm[0] != 0 || (phi_min + e_nd)*(0.5*Vesc*Vesc)*mp*mk/(mp+mk) > dmk){
            LE_histo_adapter S_HL_old(S_HL[Index],LE_func);

                auto HL_int = TrajectoryIntegral1(G1,e_nd,l_nd,TI,
                               Vesc,N_el,Therm,VescR,
                               S_HL_old,
                               E_L[Index],
                               mk,mp,-dmk,dF,PhiFactor,Nmk_HL);
                /*
                auto HL_Th = sum_lambda(NH,[&](size_t ih){return S_HL[Index].Values[ih];});
                if(HL_int == 0 || std::abs(HL_int-HL_Th) > 1e-10*(HL_int+HL_Th)){
                    COMPARE(HL_int,HL_Th);
                }*/
        }
        //auto TI_summ
        curr_progress_hl++;
        #ifndef NDEBUG
        prog_hl.show(curr_progress_hl/((float)total_progress_hl));
        #endif
    }
    prog_hl.end();


    std::cout << "calculating scatter lh" <<std::endl;
    show_prog prog_lh(100);
    prog_lh.show(0);
    size_t total_progress_lh = S_LH.size();
    size_t curr_progress_lh = 0;

    if(dmk == 0){
        prog_lh.end();
        print("skip, dmk = 0");
    } else{
        #pragma omp parallel for
        for(size_t i=0;i<NH;++i){

            auto Index = S_LH.Grid.FromLinear(i);
            auto [e0e1,l0l1] = S_LH.Grid[Index];

            if(i != S_LH.Grid.LinearIndex(Index)){
                print("linear index not matches multiindex");
                PVAR(Index);
                PVAR(i);
            }

            double e_nd=  e0e1.center();
            double l_nd=  l0l1.center()*maxLndf(phi,e_nd);//ELH_H.Lmax(e_nd);

            LE_histo_adapter S_LH_old(S_LH[Index],LE_func);
            auto TI = CalculateTrajectory(phi,e_nd,l_nd,100);
            auto LH_int =  TrajectoryIntegral1(G1,e_nd,l_nd,TI,
                               Vesc,N_el,Therm,VescR,
                               S_LH_old,
                               E_H[Index],
                               mk,mp,dmk,dF,PhiFactor,Nmk_LH);
            /*
            auto LH_Th = sum_lambda(NL,[&](size_t ih){return S_LH[Index].Values[ih];});
            if(LH_int == 0 || std::abs(LH_int-LH_Th) > 1e-10*(LH_int+LH_Th)){
                COMPARE(LH_int,LH_Th);
            }*/

            curr_progress_lh++;
            #ifndef NDEBUG
            prog_lh.show(curr_progress_lh/((float)total_progress_lh));
            #endif

        }
        prog_lh.end();
    }


}




void print_params();
int main(int argc,char **argv){

    boost::property_tree::ptree cmd_params;
    auto cmd_parse_log = parse_command_line_v1(argc,argv,cmd_params);
    if(!cmd_parse_log.empty()){
        std::cout << cmd_parse_log <<std::endl;
        return 0;
    }
    if(ptree_contain(cmd_params,"help")){
        print_params();
    }


    std::map<std::string,std::string> hl_defaut_params =
    {{"lh_out","lh.bmat"},{"hl_out","hl.bmat"}};

    std::map<std::string,std::string> evap_defaut_params =
    {{"eh_out","eh.bmat"},{"el_out","el.bmat"}};

    std::map<std::string,std::string> distrib_defaut_params =
    {{"dh_out","dh.bmat"},{"dl_out","dl.bmat"}};


    for(const auto &[key,val]:hl_defaut_params){
        if(!ptree_contain(cmd_params,key))
            cmd_params.put(key,val);
    }

    bool count_evap = false;
    bool count_distrib= false;
    MatrixFormat MF = MatrixFormat::BINARY;

    if(ptree_gets(cmd_params,"format") == "text")
        MF=MatrixFormat::TEXT;

    if(ptree_contain(cmd_params,"evap") or (ptree_contain(cmd_params,"eh_out") and
                                            ptree_contain(cmd_params,"el_out")) ){
           count_evap = true;
           for(const auto &[key,val]:evap_defaut_params){
               if(!ptree_contain(cmd_params,key))
                   cmd_params.put(key,val);
           }
    }

    if(ptree_contain(cmd_params,"distrib") or (ptree_contain(cmd_params,"dh_out") and
                                            ptree_contain(cmd_params,"dl_out")) ){
           count_distrib = true;
           for(const auto &[key,val]:distrib_defaut_params){
               if(!ptree_contain(cmd_params,key))
                   cmd_params.put(key,val);
           }
    }
    boost::filesystem::path programm_path = ptree_gets(cmd_params,"config_path");



    boost::property_tree::ptree out_params;

    out_params.put("format", (MF == MatrixFormat::BINARY) ? "bin" : "text");
    /*writing out params config file*/
    boost::filesystem::path out_path = programm_path;
    boost::filesystem::path config_out = "";
    if(ptree_contain(cmd_params,"config_out")){
        config_out = config_path_from(cmd_params.pgets("config_out"), programm_path);
        out_path = boost::filesystem::path(config_out).parent_path().lexically_normal();
    }
    boost::filesystem::path lh_out_path =
            config_path_from(cmd_params.pgets("lh_out"),out_path);
    boost::filesystem::path hl_out_path =
            config_path_from(cmd_params.pgets("hl_out"),out_path);
    out_params.put("lh_out",config_path_to(lh_out_path,out_path).string());
    out_params.put("hl_out",config_path_to(hl_out_path,out_path).string());


    boost::filesystem::path eh_out_path = "";
    boost::filesystem::path el_out_path = "";
    boost::filesystem::path dl_out_path = "";
    boost::filesystem::path dh_out_path = "";

    boost::filesystem::path grid_h_out_path = config_path_from(
                ptree_condition<std::string>(cmd_params,"h_grid","h_grid.txt"),
                out_path);
    boost::filesystem::path grid_l_out_path = config_path_from(
                ptree_condition<std::string>(cmd_params,"l_grid","l_grid.txt"),
                out_path);
    out_params.put("h_grid",config_path_to(grid_h_out_path,out_path).string());
    out_params.put("l_grid",config_path_to(grid_l_out_path,out_path).string());

    if(count_evap){
        eh_out_path = config_path_from(cmd_params.pgets("eh_out"),out_path);
        el_out_path = config_path_from(cmd_params.pgets("el_out"),out_path);
        out_params.put("eh_out",config_path_to(eh_out_path,out_path));
        out_params.put("el_out",config_path_to(el_out_path,out_path));

    }
    if(count_distrib){
        dh_out_path = config_path_from(cmd_params.pgets("dh_out"),out_path);
        dl_out_path = config_path_from(cmd_params.pgets("dl_out"),out_path);
        out_params.put("dh_out",config_path_to(dh_out_path,out_path));
        out_params.put("dl_out",config_path_to(dl_out_path,out_path));

    }


    if(config_out.string() != ""){
        std::ofstream out_config_file(config_out.string());
        boost::property_tree::write_json(out_config_file,out_params);
    }
    /*writing out params*/

    size_t Nmk_HL = ptree_condition<int>(cmd_params,"Nmk_HL",10000);
    size_t Nmk_LH = ptree_condition<int>(cmd_params,"Nmk_LH",10000);

    size_t Nmk_H = ptree_condition<int>(cmd_params,"Nmk_H",100000);
    size_t Nmk_L = ptree_condition<int>(cmd_params,"Nmk_L",100000);

    //boost::property_tree::write_json(std::cout,cmd_params);

    double l_fraction = cmd_params.get("l_fraction",-1);
    double h_fraction = cmd_params.get("h_fraction",-1);


    if(ptree_gets(cmd_params,"debug") == "true"){
        std::cout << "cmd params: ";
        boost::property_tree::write_json(std::cout,cmd_params);
        print();
        std::cout << "out params: ";
        boost::property_tree::write_json(std::cout,cmd_params);
        print();
        PVAR(programm_path);
        PVAR(out_path);
        PVAR(lh_out_path);
        PVAR(hl_out_path);
        PVAR(eh_out_path);
        PVAR(el_out_path);
        PVAR(dh_out_path);
        PVAR(dl_out_path);
        PVAR(grid_h_out_path);
        PVAR(grid_l_out_path);
        PVAR(config_out);
        #ifdef _OPENMP
        std::cout << "openmp verstion, thread num = " << omp_thread_count() << std::endl;
        #endif
    }


    std::string bm_filename;
    double Vesc1;
    try {
        bm_filename = cmd_params.get<std::string>("body");
        Vesc1 = cmd_params.get<double>("Vesc");

        if(ptree_gets(cmd_params,"debug") != ""){
            std::cout << "body filename: " << bm_filename<<std::endl;
            std::cout << "Vesc " << Vesc1<<std::endl;
        }
    } catch (std::exception &e) {
        std::cout << "can't get body model" <<std::endl;
        return 0;
    }
    auto BM = BodyModel::fromFile(Vesc1,config_path_from(bm_filename,programm_path).string());
    if(ptree_gets(cmd_params,"debug") == "true"){
        print("loaded body model");
    }


    const auto & phi = BM["phi"];
    //auto phi = Vector(BM["phi"].size(),[N = BM["phi"].size()](size_t i){return 1.5-0.5*pow(i/(N-1.0),2);});
    //PVAR(phi);


    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto G1 = [](){
        const double fc = 1.0/(RAND_MAX+1.0);
        return fc*(rand()+fc*rand());
    };

    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    auto Therm = grob::make_function<grob::linear_extrapolator> (grobGridPointCast(R),BM["Temp"]*1e-13);
    auto Vesc =  grob::make_function<grob::linear_extrapolator> (grobGridPointCast(R),BM["Vesc"]);


    if(ptree_gets(cmd_params,"debug") == "true"){
        print("loaded r functions");
    }

    double Emin = -phi[0];
    double Lmax = maxLndf(PhiC,0);


    size_t NE = ptree_condition<int>(cmd_params,"NE",30+1);
    size_t NL_max = ptree_condition<int>(cmd_params,"NLmax",30+1);

    //auto E_grid = (ptree_condition<std::string>(cmd_params,"Egird","") != "Uniform"? CreateEGrid_sqrt(NE,Emin) :
    //                                     Function::VectorGrid<double>(Emin,0.0,NE));

    auto N_L_func = [&NL_max,Lmax,&BM](double _E){return 2+(NL_max-2)*maxLnd(BM,_E)/Lmax;};
    auto Lgridrho = [](double _E,double _L_nd){return 1.0;};
    //auto HistoGrid = Create_EL_grid(E_grid,Lgridrho,N_L_func);
    //if(ptree_condition<std::string>(cmd_params,"Lgird","") == "Uniform"){
    //    HistoGrid = decltype (HistoGrid)(E_grid,Function::VectorGrid<double>(0.,1.,NL_max));
    //}

    //PVAR(HistoGrid.gridStr());
    //exit(0);

    //EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>> ELH_L
    //        (HistoGrid,Function::FunctorM(BM,maxLnd));
    //auto ELH_H = ELH_L;

    //std::ofstream gh_out(grid_h_out_path.string());
    //PVAR(gh_out.is_open());

    //ELH_H.save_text(gh_out);
    //ELH_L.save_text(std::ofstream(grid_l_out_path.string()));

    auto LE_func = grob::make_function_f(grob::GridUniform<double>(Emin,0,4*NE),
                          [&PhiC](double _E){return maxLndf(PhiC,_E);}
                          );

    WriteObject(LE_func,std::ofstream(
                    config_path_from(cmd_params.pgets("LE_json","le.json.txt"),out_path)
                    ));
    grob::as_csv(LE_func).save(std::ofstream(
                                   config_path_from(
                                       cmd_params.pgets("LE","le.txt"),
                                       out_path)
                                   ));

    // *******************
    // grob implementation
    // *******************
    // Creating 2dim grids
    /*
    auto GridH = grobGetHistoGrid((decltype(ELH_H)::HBase const&)(ELH_H));

    // Creating 2dim grids
    auto GridL = grobGetHistoGrid((decltype(ELH_H)::HBase const&)(ELH_H));

    */
    double ineqE_a = cmd_params.get<double>("Egrid.ineq_a",
                                    cmd_params.get<double>("Egrid.ineq",0.0));

    double ineqE_b = cmd_params.get<double>("Egrid.ineq_b",
                                    cmd_params.get<double>("Egrid.ineq",0.0));

    double ineqL_b = cmd_params.get<double>("Lgrid.ineq",0.0);

    auto GridE = grob::make_histo_grid(CreateSqrtGrid(NE,Emin,0,
                                1,
                                ineqE_a,
                                ineqE_b));
    auto GridH = grob::make_grid_f(GridE,[&](size_t i){
            double _E = GridE[i].left();
        return grob::make_histo_grid(
                CreateSqrtGrid(N_L_func(_E),0.0,1.0,1,0.0,ineqL_b)
                );
    });
    WriteObject(GridH,std::ofstream(
                    config_path_from(cmd_params.pgets("h_grid","h_grid.txt"),out_path)
                    ));
    auto GridL = GridH;
    WriteObject(GridL,std::ofstream(
                    config_path_from(cmd_params.pgets("l_grid","l_grid.txt"),out_path)
                    ));
    auto NH = GridH.size();
    auto NL = GridL.size();
    // Creating 2dim histos for capture
    auto Cap_H = grob::make_histo_ref<double>(GridH);
    auto Cap_L = grob::make_histo_ref<double>(GridL);

    // Creating 2dim histos for evaporation
    auto Evap_H = grob::make_histo_ref<double>(GridH);
    auto Evap_L = grob::make_histo_ref<double>(GridL);

    // Creating matrix layout for MatrixElements
    std::vector<double> HL_Mat(NH*NL,0);
    std::vector<double> LH_Mat(NH*NL,0);

    // make histo of histos, which views

    auto S_LH = grob::make_histo_ref(GridH,grob::make_vector(NH,
                         [&](size_t i){
                            return grob::make_histo_ref(GridL,
                                grob::make_slice(LH_Mat,i*NL,NL)
                             );
                         }
        ));
    auto S_HL =  grob::make_histo_ref(GridL,grob::make_vector(NL,
                         [&](size_t i){
                            return grob::make_histo_ref(GridH,
                                grob::make_slice(HL_Mat,i*NH,NH)
                             );
                         }
        ));

    // *******************
    // grob implementation
    // *******************

    double mk = ptree_condition(cmd_params,"mk",5.0);
    double delta_mk = ptree_condition(cmd_params,"dmk",1e-6);

    std::vector<std::string> ElementList;
    try {
        auto &Elements =
                cmd_params.get_child("elements");
        for(const auto &el : Elements){
            ElementList.push_back(el.second.data());
        }
    }  catch (std::exception &e) {
        print("error in element list");
    }
    PVAR(ElementList);
    auto RhoND = BM["Rho"];
    //Therm.Values*=0;
    bool not_fill_ss = ptree_condition(cmd_params,"not_fill",false);

    double V_body = ptree_condition(cmd_params,"Vbody",0.73e-3);
    double V_disp = ptree_condition(cmd_params,"Vdisp",0.52e-3);

    if(ptree_gets(cmd_params,"debug") != ""){
        PVAR(V_body);
        PVAR(V_disp);
    }

    double FullCap_H = 0;
    double FullCap_L = 0;
    double VescMin = BM.VescMin();
    boost::property_tree::ptree elements_portions;
    for(const auto & element : ElementList){
        auto CalcCapture = [&](auto const & dF_Factor,
                auto const & mPhiFactor,auto const & N_el_conc,
                double m_nuc){
            LE_histo_adapter Cap_H_old(Cap_H,LE_func);
            LE_histo_adapter Cap_L_old(Cap_L,LE_func);
            double numEl_H = 0;
            if(delta_mk != 0){
                numEl_H = SupressFactor_v1(Cap_H_old,m_nuc,mk,-delta_mk,dF_Factor,N_el_conc,BM.VescMin(),Vesc,Therm,
                                 mPhiFactor,G,Nmk_H,1.5,V_disp,V_body);
            }
            std::cout << "H fraction for " << element << " = " << numEl_H << std::endl;
            double numEl_L = SupressFactor_v1(Cap_L_old,m_nuc,mk,delta_mk,dF_Factor,N_el_conc,BM.VescMin(),Vesc,Therm,
                             mPhiFactor,G,Nmk_L,1.5,V_disp,V_body);
            std::cout << "L fraction for " << element << " = " << numEl_L << std::endl;

            FullCap_H += numEl_H;
            FullCap_L += numEl_L;

            elements_portions.put("L."+element,numEl_L);
            elements_portions.put("H."+element,numEl_H);
        };
        if(element == "H_e_p"){
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM["H"]);
            dF_Nuc dF(1,1);
            Phi_Fac_S PhiFac(_mp,mk);
            CalcCapture(dF,PhiFac,Element_N,_mp);
        } else if (element == "H_e_e") {
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM["H"]);
            dF_H_Elastic_Electron dF;
            Phi_Fac_Electron_S PhiFac(_mp,mk);
            CalcCapture(dF,PhiFac,Element_N,_mp);
        } else if (element == "H_m_p") {
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM["H"]);
            dF_H_Migdal_Proton dF(G);
            Phi_Fac_S PhiFac(_mp,mk);
            CalcCapture(dF,PhiFac,Element_N,_mp);
        } else if (element == "H_m_e") {
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM["H"]);
            dF_H_Migdal_Electron dF(G);
            Phi_Fac_Electron_S PhiFac(_mp,mk);
            CalcCapture(dF,PhiFac,Element_N,_mp);
        } else if (element == "H_i_p") {
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM["H"]);
            dF_H_Ion_Proton dF(G);
            Phi_Fac_S PhiFac(_mp,mk);
            CalcCapture(dF,PhiFac,Element_N,_mp);
        } else if (element == "H_i_e") {
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM["H"]);
            dF_H_Ion_Electron dF(G);
            Phi_Fac_Electron_S PhiFac(_mp,mk);
            CalcCapture(dF,PhiFac,Element_N,_mp);
        } else {
            decltype(Therm) Element_N(Therm.Grid,RhoND*BM[element]);

            std::cout << "calculating for element " << element <<std::endl;
            size_t M = ME.at(element);
            size_t Z = QE.at(element);
            double m_nuc = (M)*_mp;
            dF_Nuc dF(M,Z);
            Phi_Fac_S PhiFactor_Nuc(m_nuc,mk);
            //auto PhiFactor_Nuc_55 = [M,mk,VescMin](double q){return fast_pow(M*q/(mk*VescMin),4);};
            if(!not_fill_ss){
                FillScatterHistoFromElement(Element_N,Therm,PhiC,BM.VescMin(),Vesc,dF,
                                       mk,delta_mk,m_nuc,
                                       S_LH,S_HL,Evap_H,Evap_L,LE_func,PhiFactor_Nuc,Nmk_HL,Nmk_LH);
            }

            if(count_distrib){
                CalcCapture(dF,PhiFactor_Nuc,Element_N,m_nuc);
            }
        }
    }
    if(ptree_contain(cmd_params,"fractions")){
        std::ofstream fraction_file(config_path_from(cmd_params.pgets("fractions"),out_path).string());
        boost::property_tree::write_json(fraction_file,elements_portions);
    }

    if(contains_nan(HL_Mat)){
        std::cout << "HL_Mat contains nan" <<std::endl;
    }
    if(contains_nan(LH_Mat)){
        std::cout << "LH_Mat contains nan" <<std::endl;
    }

    if(contains_nan(Evap_H.Values)){
        std::cout << "Ev_H contains nan" <<std::endl;
    }
    if(contains_nan(Evap_L.Values)){
        std::cout << "Ev_L contains nan" <<std::endl;
    }

    if(contains_nan(Cap_H.Values)){
        std::cout << "Cap_H contains nan" <<std::endl;
    }
    if(contains_nan(Cap_L.Values)){
        std::cout << "Cap_L contains nan" <<std::endl;
    }

    double c_plus = std::max(FullCap_H,FullCap_L);
    PVAR(c_plus);

    /*
    Cap_H.Values /= c_plus;
    Cap_L.Values /= c_plus;
    HL_Mat /= c_plus;
    LH_Mat /= c_plus;
    Evap_H.Values /= c_plus;
    Evap_L.Values /= c_plus;
    */
    boost::property_tree::ptree OutParams;

    OutParams.put("cp",c_plus);
    OutParams.put("nuh",FullCap_H/c_plus);
    OutParams.put("nul",FullCap_L/c_plus);

    OutParams.put("eh",max(Evap_H.Values));
    OutParams.put("el",max(Evap_L.Values));

    auto max_val = [](size_t size, auto &&ValLambda){
        auto x0 = ValLambda(0);
        for(size_t i=1;i<size;++i){
            x0 = std::max(x0,ValLambda(i));
        }
        return x0;
    };

    auto S_HL_out = Vector(NL,[&](size_t i){return std::accumulate(&HL_Mat[i*NH],&HL_Mat[(i+1)*NH],0.0);});
    auto S_LH_out = Vector(NH,[&](size_t i){return std::accumulate(&LH_Mat[i*NL],&LH_Mat[(i+1)*NL],0.0);});

    OutParams.put("st_hl",max_val(NL,[&](size_t i){return S_HL_out[i];}));
    OutParams.put("st_lh",max_val(NH,[&](size_t i){return S_LH_out[i];}));

    boost::property_tree::write_json(config_path_from(cmd_params.pgets("pout","out.txt"),out_path).string(),OutParams);

    if(count_distrib){
        saveMatrix(Cap_H.Values.data(),1,NH,dh_out_path.string(),MF);
        saveMatrix(Cap_L.Values.data(),1,NL,dl_out_path.string(),MF);
        OutParams.put("CapH",config_path_to(dh_out_path,out_path));
        OutParams.put("CapL",config_path_to(dl_out_path,out_path));
    }
    if(count_evap){
        saveMatrix(Evap_H.Values.data(),NH,1,eh_out_path.string(),MF);
        saveMatrix(Evap_L.Values.data(),NL,1,el_out_path.string(),MF);
        OutParams.put("EvapH",config_path_to(eh_out_path,out_path));
        OutParams.put("EvapL",config_path_to(el_out_path,out_path));
    }

    if(!not_fill_ss){
        saveMatrix(LH_Mat.data(),NH,NL,lh_out_path.string(),MF);
        saveMatrix(HL_Mat.data(),NL,NH,hl_out_path.string(),MF);
        saveMatrix(S_HL_out.data(),1,NL,config_path_from(cmd_params.pgets("o_hl","o_hl.bmat"),out_path).string(),MF);
        saveMatrix(S_LH_out.data(),1,NH,config_path_from(cmd_params.pgets("o_lh","o_lh.bmat"),out_path).string(),MF);
    }

    return 0;
}


void print_params(){
    printd('\n',
           "mk : [mass of WIMP, GeV]",
           "dmk : [optinal, (default 0) delta mass of WIMP, GeV]",
           "body : [filename of celestial bject]",
           "Vesc : [first consmic velocity of body]",
           "fractions : [optional, filename to put fractions of elements in capture rate]",
           "Vbody : [optional, velocity of object, default is 0.73e-3]",
           "Vdisp : [optional, velocity in halo distribution, default is 0.52e-3]",
           "LE : [output path to save dat file of Lmax(E) function]",
           "LE_json : [output path to save json file of Lmax(E) function]",
           "pout : [output path to save full rates (capture, scatter)]",
           "Egrid.ineq_a : [optional, default 0, makes grid more tight in E ~ Emin]",
           "Egrid.ineq_b : [optional, default 0, makes grid more tight in E ~ Emax]",
           "Lgrid.ineq : [optional, default 0, makes grid more tight in L ~ Lmax]",
           "NE : [optional, default 31, size of E grid",
           "NLmax : [optional, default 31, max size of L grid",
           "elements : [array of elements to consider]",
           "Nmk_H, Nmk_L : [Monte-Carlo steps in capture rate in H and L resp.]",
           "Nmk_HL, Nmk_LL : [Monte-Carlo steps in scatter HL and LH  resp.]");
}
