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



int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}


auto PhiFactorSS = [](double mk,double q)->double{
    return 1.0;
};

Function::VectorGrid<double> CreateEGrid_sqrt(size_t N,double Emin,double eps = 1e-4){
    return  DensityGrid(Emin,0,[Emin,eps](double x){
        double t = (x-Emin)/std::abs(Emin);
        return 1/std::abs(Emin)*0.5*( 0.5/sqrt(t+eps)+0.5/sqrt(1-t+eps));
    },N);
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


struct dF_Nuc{
    size_t N,Z;
    double const_factor;
    dF_Nuc(size_t Z = 1,size_t  N = 0,double const_factor = 1):Z(Z),N(N),const_factor(const_factor){}
    MC::MCResult<double> EnergyLoss(double Emax)const{
        return MC::MCResult<double>(0,1);
    }
    double ScatterFactor(double q,double enLoss)const{
        return const_factor;
    }
};

struct dF_Nuc_M{
    size_t M;
    double const_factor;
    dF_Nuc_M(size_t M=0,double const_factor = 1):M(M),const_factor(const_factor){}
    inline MC::MCResult<double> EnergyLoss(double Emax)const noexcept{
        return MC::MCResult<double>(0,1);
    }
    inline double ScatterFactor(double q,double enLoss)const noexcept{
        return M*M;
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

template <typename N_type,typename Therm_type,typename Phi_Type,typename VescFuncType,
          typename dF_Type,typename EvapHisto,typename ScatterHisto,
          typename PhiFactorType>
inline void FillScatterHisoFromElement(const N_type & N_el,Therm_type const & Therm,Phi_Type const& phi,
                                       double Vesc, const VescFuncType & VescR,
                                       const dF_Type&dF,double mk,double dmk,double mp,
                                ScatterHisto &S_HL,ScatterHisto &S_LH,EvapHisto & E_H,EvapHisto & E_L,
                                const PhiFactorType & PhiFactor,size_t Nmk_HL,size_t Nmk_LH){
    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto G1 = [](){
        double fc = 1.0/(RAND_MAX+1.0);
        return fc*(rand()+fc*rand());
    };
    #ifdef _OPENMP
    random_shuffler rs_hl(S_HL.Grid.size()-1);
    #else
    random_shuffler rs_hl(S_HL.Grid.size()-1,1);
    #endif

    std::cout << "calculating scatter hl" <<std::endl;
    show_prog prog_hl(100);
    prog_hl.show(0);
    size_t total_progress_hl = S_HL.num_of_element();
    size_t curr_progress_hl = 0;

    const size_t NE_HL = S_HL.Grid.size()-1;
    #pragma omp parallel for
    for(size_t _i=0;_i<NE_HL;++_i){
        size_t i = rs_hl[_i];
        double e_nd=  0.5*(S_HL.Grid[i]+S_HL.Grid[i+1]);
        for(size_t j = 0;j<S_HL.values[i].Grid.size()-1;++j){
            double l_nd = 0.5*(S_HL.values[i].Grid[j]+S_HL.values[i].Grid[j])*
                    maxLndf(phi,e_nd);//ELH_H.Lmax(e_nd);

            auto TI = CalculateTrajectory(phi,e_nd,l_nd,100);

            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               Vesc,N_el,Therm,VescR,
                               S_HL.values[i].values[j],
                               E_H.values[i].values[j],
                               mk,mp,dmk,dF,PhiFactor,Nmk_HL);

            curr_progress_hl++;
            #ifndef NDEBUG
            prog_hl.show(curr_progress_hl/((float)total_progress_hl));
            #endif
        }
    }
    prog_hl.end();

    #ifdef _OPENMP
    random_shuffler rs_lh(S_LH.Grid.size()-1);
    #else
    random_shuffler rs_lh(S_LH.Grid.size()-1,1);
    #endif

    std::cout << "calculating scatter lh" <<std::endl;


    show_prog prog_lh(100);
    prog_lh.show(0);
    size_t total_progress_lh = S_LH.num_of_element();
    size_t curr_progress_lh = 0;

    const size_t NE_LH = S_LH.Grid.size()-1;

    #pragma omp parallel for
    for(size_t _i=0;_i<NE_LH;++_i){
        size_t i = rs_lh[_i];
        double e_nd=  0.5*(S_LH.Grid[i]+S_LH.Grid[i+1]);
        for(size_t j = 0;j<S_LH.values[i].Grid.size()-1;++j){
            double l_nd = 0.5*(S_LH.values[i].Grid[j]+S_LH.values[i].Grid[j])*
                    maxLndf(phi,e_nd);//ELH_H.Lmax(e_nd);
            auto TI = CalculateTrajectory(phi,e_nd,l_nd,100);
            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               Vesc,N_el,Therm,VescR,
                               S_LH.values[i].values[j],
                               E_L.values[i].values[j],
                               mk,mp,-dmk,dF,PhiFactor,Nmk_LH);

            curr_progress_lh++;
            #ifndef NDEBUG
            prog_lh.show(curr_progress_lh/((float)total_progress_lh));
            #endif
        }
    }
    prog_lh.end();
}





int main(int argc,char **argv){

    boost::property_tree::ptree cmd_params;
    auto cmd_parse_log = parse_command_line_v1(argc,argv,cmd_params);
    if(!cmd_parse_log.empty()){
        std::cout << cmd_parse_log<<std::endl;
        return 1;
    }


    std::map<std::string,std::string> hl_defaut_params =
    {{"lh_out","lh.bmat"},{"hl_out","hl.mat"}};

    std::map<std::string,std::string> evap_defaut_params =
    {{"eh_out","eh.bmat"},{"el_out","el.mat"}};

    std::map<std::string,std::string> distrib_defaut_params =
    {{"dh_out","dh.bmat"},{"dl_out","dl.mat"}};


    for(const auto &[key,val]:hl_defaut_params){
        if(!ptree_contain(cmd_params,key))
            cmd_params.put(key,val);
    }

    bool count_evap = false;
    bool count_distrib= false;
    MatrixFormat MF = BINARY;

    if(ptree_gets(cmd_params,"format") == "text")
        MF=TEXT;

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

    out_params.put("format", (MF == BINARY) ? "bin" : "text");
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




    if(ptree_gets(cmd_params,"debug") != ""){
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
    auto BM = BodyModel::fromFile(Vesc1,config_path_from(bm_filename,programm_path));



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

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-13);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);


    double Emin = -phi[0];
    double Lmax = maxLndf(PhiC,0);



    size_t NE = ptree_condition<int>(cmd_params,"NE",30+1);
    size_t NL_max = ptree_condition<int>(cmd_params,"NLmax",30+1);
    auto E_grid = (ptree_condition<std::string>(cmd_params,"Egird","") == "Uniform"? CreateEGrid_sqrt(NE,Emin) :
                                          Function::VectorGrid<double>(Emin,0.0,NE));
    auto N_L_func = [&NL_max,Lmax,&BM](double _E){return 2+(NL_max-2)*maxLnd(BM,_E)/Lmax;};
    auto Lgridrho = [](double _E,double _L_nd){return 1.0;};
    auto HistoGrid = Create_EL_grid(E_grid,Lgridrho,N_L_func);
    if(ptree_condition<std::string>(cmd_params,"Lgird","") == "Uniform"){
        HistoGrid = decltype (HistoGrid)(E_grid,Function::VectorGrid<double>(0.,1.,NL_max));
    }


    //PVAR(HistoGrid.gridStr());
    //exit(0);

    EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>> ELH_L
            (HistoGrid,Function::FunctorM(BM,maxLnd));
    auto ELH_H = ELH_L;

    std::ofstream gh_out(grid_h_out_path.string());
    PVAR(gh_out.is_open());
    ELH_H.save_text(gh_out);
    ELH_L.save_text(std::ofstream(grid_l_out_path.string()));

    auto EvapHisto_H = ELH_H;
    auto EvapHisto_L = ELH_L;

    //print("initialising");
    auto S_LH =  Function::GridExtractorType<decltype(ELH_H),decltype(ELH_L)::HBase>::type::sameGrid(ELH_L);
    auto S_HL =  Function::GridExtractorType<decltype(ELH_L),decltype(ELH_H)::HBase>::type::sameGrid(ELH_H);
    //print("initialised");

    for(size_t i=0;i<S_HL.num_of_element();++i){
        S_HL[i] = ELH_L;
    }
    for(size_t i=0;i<S_LH.num_of_element();++i){
        S_LH[i] = ELH_H;
    }

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

    }
    PVAR(ElementList);
    auto RhoND = BM["RhoND"];

    bool not_fill_ss = ptree_condition(cmd_params,"not_fill",false);
    for(const auto & element : ElementList){
        decltype(Therm) Element_N(R,RhoND*BM[element]);
        std::cout << "calculating for element " << element <<std::endl;
        size_t M = ME.at(element);
        double m_nuc = (M)*_mp;
        dF_Nuc_M dF(M);
        if(!not_fill_ss){
            FillScatterHisoFromElement(Element_N,Therm,PhiC,BM.VescMin(),Vesc,dF,
                                   mk,delta_mk,m_nuc,
                                   S_HL,S_LH,EvapHisto_H,EvapHisto_L,PhiFactorSS,Nmk_HL,Nmk_LH);
        }
        if(count_distrib){
            SupressFactor_v1(ELH_H,m_nuc,mk,-delta_mk,dF,Element_N,BM.VescMin(),Vesc,Therm,
                             PhiFactorSS,G,Nmk_H,2,U0,U0);
        }
    }


    size_t N_H = ELH_H.num_of_element();
    size_t N_L = ELH_L.num_of_element();

    std::vector<double> MatT_HL(N_H*N_L);
    std::vector<double> MatT_LH(N_H*N_L);


    if(contains_nan(MatT_HL)){
        std::cout << "MatT_HL contains nan" <<std::endl;
    }
    if(contains_nan(MatT_LH)){
        std::cout << "MatT_LH contains nan" <<std::endl;
    }


    for(size_t i=0;i<N_L;++i){
        S_LH[i].saveIter(MatT_LH.data()+i*N_H,MatT_LH.data()+(i+1)*N_H);
    }
    for(size_t i=0;i<N_H;++i){
        S_HL[i].saveIter(MatT_HL.data()+i*N_L,MatT_HL.data()+(i+1)*N_L);
    }


    if(!not_fill_ss){
        saveMatrix(MatT_HL.data(),N_H,N_L,hl_out_path.string(),MF);
        saveMatrix(MatT_LH.data(),N_L,N_H,lh_out_path.string(),MF);
    }
    std::vector<double> Ev_H = EvapHisto_H.AllValues();
    std::vector<double> Ev_L = EvapHisto_L.AllValues();

    if(contains_nan(Ev_H)){
        std::cout << "Ev_H contains nan" <<std::endl;
    }
    if(contains_nan(Ev_L)){
        std::cout << "Ev_L contains nan" <<std::endl;
    }

    if(count_evap){
        saveMatrix(Ev_H.data(),N_H,1,el_out_path.string(),MF);
        saveMatrix(Ev_L.data(),N_L,1,eh_out_path.string(),MF);
    }

    std::vector<double> D_H;// = ELH_H.AllValues();
    std::vector<double> D_L;// = ELH_L.AllValues();



    if(count_distrib){
        D_H = ELH_H.AllValues();
        D_L = ELH_L.AllValues();
        saveMatrix(D_H.data(),1,N_H,dh_out_path.string(),MF);
        saveMatrix(D_L.data(),1,N_L,dl_out_path.string(),MF);
    }
    if(contains_nan(D_H)){
        std::cout << "D_H contains nan" <<std::endl;
    }
    if(contains_nan(D_L)){
        std::cout << "D_L contains nan" <<std::endl;
    }


    return 0;
}
