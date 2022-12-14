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
    dF_Nuc(size_t Z = 1,size_t  N = 0):Z(Z),N(N){}
    MC::MCResult<double> EnergyLoss(double Emax)const{
        return MC::MCResult<double>(0,1);
    }
    double ScatterFactor(double q,double enLoss)const{
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

template <typename N_type,typename Therm_type,typename Phi_Type,typename VescFuncType,
          typename dF_Type,typename EvapHisto,typename ScatterHisto,
          typename PhiFactorType>
inline void FillScatterHisoFromElement(const N_type & N_el,Therm_type const & Therm,Phi_Type const& phi,
                                       double Vesc, const VescFuncType & VescR,
                                       const dF_Type&dF,double mk,double dmk,double mp,
                                ScatterHisto &S_HL,ScatterHisto &S_LH,EvapHisto & E_H,EvapHisto & E_L,
                                const PhiFactorType & PhiFactor,size_t Nmk_HL=10000,size_t Nmk_LH = 10000){
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

    #pragma omp parallel for
    for(size_t _i=0;_i<S_HL.Grid.size()-1;++_i){
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

        }
    }

    #ifdef _OPENMP
    random_shuffler rs_lh(S_LH.Grid.size()-1);
    #else
    random_shuffler rs_lh(S_LH.Grid.size()-1,1);
    #endif

    #pragma omp parallel for
    for(size_t _i=0;_i<S_LH.Grid.size()-1;++_i){
        size_t i = rs_lh[i];
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

        }
    }

}





int main(int argc,char **argv){

    boost::property_tree::ptree cmd_params;
    auto cmd_parse_log = parse_command_line_v1(argc,argv,cmd_params);


    boost::property_tree::write_json(std::cout,cmd_params);
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

    if(ptree_gets(cmd_params,"mode") == "text")
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


    /*writing out params config file*/
    boost::filesystem::path out_path = programm_path;
    boost::filesystem::path config_out = "";
    if(ptree_contain(cmd_params,"config_out")){
        config_out = config_path_from(cmd_params.pgets("config_out"), programm_path);
        out_path = boost::filesystem::path(config_out).parent_path();
    }
    boost::filesystem::path lh_out_path =
            config_path_from(cmd_params.pgets("lh_out"),out_path);
    boost::filesystem::path hl_out_path =
            config_path_from(cmd_params.pgets("lh_out"),out_path);
    out_params.put("lh_out",config_path_to(lh_out_path,out_path));
    out_params.put("hl_out",config_path_to(hl_out_path,out_path));

    boost::filesystem::path eh_out_path = "";
    boost::filesystem::path el_out_path = "";
    boost::filesystem::path dl_out_path = "";
    boost::filesystem::path dh_out_path = "";

    boost::filesystem::path grid_h_out_path = config_path_from(
                ptree_contain(cmd_params,"h_grid") ? "h_grid.txt" : cmd_params.pgets("h_grid"),
                out_path);
    boost::filesystem::path grid_l_out_path = config_path_from(
                ptree_contain(cmd_params,"l_grid") ? "l_grid.txt" : cmd_params.pgets("l_grid"),
                out_path);
    out_params.put("h_grid",config_path_to(grid_h_out_path,out_path));
    out_params.put("l_grid",config_path_to(grid_l_out_path,out_path));

    if(count_evap){
        eh_out_path = config_path_from(cmd_params.pgets("eh_out"),out_path);
        el_out_path = config_path_from(cmd_params.pgets("el_out"),out_path);
        out_params.put("eh_out",config_path_to(eh_out_path,out_path));
        out_params.put("el_out",config_path_to(el_out_path,out_path));

    }
    if(count_distrib){
        dh_out_path = config_path_from(cmd_params.pgets("dh_out"),out_path);
        dh_out_path = config_path_from(cmd_params.pgets("dl_out"),out_path);
        out_params.put("dh_out",config_path_to(dh_out_path,out_path));
        out_params.put("dl_out",config_path_to(dl_out_path,out_path));

    }


    if(config_out.string() != ""){
        std::ofstream out_config_file(config_out.string());
        boost::property_tree::write_json(out_config_file,out_params);
    }
    /*writing out params*/

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
        PVAR(config_out);
    }


    std::string bm_filename;
    double Vesc1= 2.056e-3;
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



    size_t NE = 30+1;
    size_t NL_max = 30+1;
    auto E_grid = CreateEGrid_sqrt(31,Emin);
    auto N_L_func = [&NL_max,Lmax,&BM](double _E){return 2+(NL_max-2)*maxLnd(BM,_E)/Lmax;};
    auto Lgridrho = [](double _E,double _L_nd){return 1.0;};
    auto HistoGrid = Create_EL_grid(E_grid,Lgridrho,N_L_func);


    PVAR(HistoGrid.gridStr());
    //exit(0);

    EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>> ELH_L
            (HistoGrid,Function::FunctorM(BM,maxLnd));
    auto ELH_H = ELH_L;

    ELH_H.save_text(std::ifstream(grid_h_out_path));
    ELH_L.save_text(std::ifstream(grid_l_out_path));

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

    double mk = 5;
    double delta_mk = 1e-6;


    double f_p_av = 4e-7;
    double min_atom_r = 0.7;
    /*
    decltype(Therm) F_e_r(R,Vector(R.size(),[&BM,&Therm,f_p_av](size_t i){
        double rho = BM["RhoND"][i];
        double s = pow(2*M_PI*Therm[i]/me,1.5)*exp(-Rd/Therm[i]);
        return f_p_av*rho*Ion_Deg(f_p_av*rho,s);
    }));

    decltype(Therm) atom_degree(R,Vector(R.size(),[&BM,&Therm,f_p_av](size_t i){
        double rho = BM["RhoND"][i];
        double s = pow(2*M_PI*Therm[i]/me,1.5)*exp(-Rd/Therm[i]);
        return Atom_ion_Deg(f_p_av*rho,s);
    }));

    decltype(Therm) ion_degree(R,Vector(R.size(),[&BM,&Therm,f_p_av](size_t i){
        double rho = BM["RhoND"][i];
        double s = pow(2*M_PI*Therm[i]/me,1.5)*exp(-Rd/Therm[i]);
        return Ion_Deg(f_p_av*rho,s);
    }));

    decltype(Therm) H_N(R,BM["RhoND"]*BM["H1"]);
    auto Atom_N =  H_N*atom_degree;
    */
    size_t Nmk_HL = 10000;
    size_t Nmk_LH = 10000;

    size_t Nmk_H = 1000000;
    size_t Nmk_L = 1000000;

    std::vector<std::string> ElementList{"H","He3","He"};
    for(const auto & element : ElementList){
        decltype(Therm) Element_N(R,BM["RhoND"]*BM[element]);
        size_t Z = 1;
        size_t N = 0;
        double m_nuc = (Z+N)*_mp;
        dF_Nuc dF(Z,N);
        FillScatterHisoFromElement(Element_N,Therm,PhiC,BM.VescMin(),Vesc,dF,
                                   mk,delta_mk,m_nuc,
                                   S_HL,S_LH,EvapHisto_H,EvapHisto_L,PhiFactorSS,Nmk_HL,Nmk_LH);
        if(count_distrib){
            SupressFactor_v1(ELH_H,m_nuc,mk,-delta_mk,dF,Element_N,BM.VescMin(),Vesc,Therm,
                             PhiFactorSS,G1,Nmk_H,U0,U0);
        }
    }


    size_t N_H = ELH_H.num_of_element();
    size_t N_L = ELH_L.num_of_element();

    std::vector<double> MatT_HL(N_H*N_L);
    std::vector<double> MatT_LH(N_H*N_L);

    for(size_t i=0;i<N_L;++i){
        S_LH[i].saveIter(MatT_LH.data()+i*N_H,MatT_LH.data()+(i+1)*N_H);
    }
    for(size_t i=0;i<N_H;++i){
        S_HL[i].saveIter(MatT_HL.data()+i*N_L,MatT_HL.data()+(i+1)*N_L);
    }

    saveMatrix(MatT_HL.data(),N_H,N_L,hl_out_path.string(),MF);
    saveMatrix(MatT_LH.data(),N_L,N_H,lh_out_path.string(),MF);

    std::vector<double> Ev_H = EvapHisto_H.AllValues();
    std::vector<double> Ev_L = EvapHisto_L.AllValues();

    if(count_evap){
        saveMatrix(Ev_H.data(),N_H,1,el_out_path.string(),MF);
        saveMatrix(Ev_L.data(),N_L,1,eh_out_path.string(),MF);
    }

    std::vector<double> D_H;// = ELH_H.AllValues();
    std::vector<double> D_L;// = ELH_L.AllValues();

    if(count_distrib){
        D_H = ELH_H.AllValues();
        D_L = ELH_L.AllValues();
        saveMatrix(D_H.data(),N_H,1,dl_out_path.string(),MF);
        saveMatrix(D_L.data(),N_L,1,dh_out_path.string(),MF);
    }


    return 0;
}
