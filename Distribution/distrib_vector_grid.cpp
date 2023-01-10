#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/arg_parser.hpp"

const double mp = 0.938;
const double me = 0.51e-3;
auto PhiFactor55p = [](double mk,double q)->double{
    return fast_pow(q*q/(mp*mk),2);
};
auto PhiFactor55e = [](double mk,double q)->double{
    return fast_pow(q*q/(me*mk),2);
};
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
    MC::MCResult<double> EnergyLoss(double Emax){
        return MC::MCResult<double>(0,1);
    }
    double ScatterFactor(double q){
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
          typename dF_Type,typename EvapHisto,
          typename DistribHisto,typename ScatterHisto,
          typename PhiFactorType>
inline void FillCaptureHisoFromElement(const N_type & N_el,Therm_type const & Therm,Phi_Type const& phi,
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
    //out_params.put("out_path",cmd_params.get<std::string>("config_path"));
    /*
    for(const auto & item:out_defaut_params){
        std::regex file_param("(.+)_out");
        std::smatch par;
        if(std::regex_match(item.first,par,file_param)){
            out_params.put(std::string(par[1])+"_in",item.second);
        }
    }
    */

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

    bool sem = true;
    double workSize=S_HL.num_of_element();
    show_prog PrgBar(100,"histo fill hl progress");
    size_t work_progress_tmp = 0;
    PrgBar.show(0.0);
    using namespace std::string_literals;
    //print("end\t","l_nd\t","T_in\t","T_out\t","T_all\t","T_all_th");
    for(size_t i=0;i<S_HL.Grid.size()-1;++i){
        double e_nd=  0.5*(S_HL.Grid[i]+S_HL.Grid[i+1]);
        for(size_t j = 0;j<S_HL.values[i].Grid.size()-1;++j){
            double l_nd = 0.5*(S_HL.values[i].Grid[j]+S_HL.values[i].Grid[j])*
                    maxLnd(BM,e_nd);//ELH_H.Lmax(e_nd);
            auto TI = CalculateTrajectory(PhiC,e_nd,l_nd,100);

            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               S_HL.values[i].values[j],
                               EvapHisto_H.values[i].values[j],
                               mk,mp,delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);

            ++work_progress_tmp;
            PrgBar.show(work_progress_tmp/((float)workSize));
        }
    }
    PrgBar.end();
    show_prog PrgBar1(100,"histo fill lh progress");
    PrgBar1.show(0);

    double workSize1=S_LH.num_of_element();
    size_t work_progress_tmp1 = 0;

    for(size_t i=0;i<S_LH.Grid.size()-1;++i){
        double e_nd=  0.5*(S_LH.Grid[i]+S_LH.Grid[i+1]);
        for(size_t j = 0;j<S_LH.values[i].Grid.size()-1;++j){
            double l_nd = 0.5*(S_LH.values[i].Grid[j]+S_LH.values[i].Grid[j])*
                    maxLnd(BM,e_nd);//ELH_H.Lmax(e_nd);
            auto TI = CalculateTrajectory(PhiC,e_nd,l_nd,100);
            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               S_LH.values[i].values[j],
                               EvapHisto_L.values[i].values[j],
                               mk,mp,-delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);
            ++work_progress_tmp1;
            PrgBar1.show(work_progress_tmp1/((float)workSize1));
        }
    }
    PrgBar1.end();

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

    std::vector<double> Ev_H = EvapHisto_H.AllValues();
    std::vector<double> Ev_L = EvapHisto_L.AllValues();






    //PVAR(ScatterPreHisto_LH.Composition([](auto x){return x.summ();}).toString());
    //PVAR(ScatterPreHisto1_LH.Composition([](auto x){return x.summ();}).toString());
    /*
    auto HTV_LH = ScatterPreHisto_LH.Composition([](auto x){return x.AllValues();});
    auto HTV_HL = ScatterPreHisto_HL.Composition([](auto x){return x.AllValues();});

    auto HTV1_LH = ScatterPreHisto1_LH.Composition([](auto x){return x.AllValues();});
    auto HTV1_HL = ScatterPreHisto1_HL.Composition([](auto x){return x.AllValues();});

    auto VecHisto_LH = Function::GridExtractorType<std::vector<double>,decltype (ELH_L)::HBase>::type::sameGrid(ELH_L);
    auto VecHisto_HL = Function::GridExtractorType<std::vector<double>,decltype (ELH_H)::HBase>::type::sameGrid(ELH_H);

    auto VecHisto1_LH = Function::GridExtractorType<std::vector<double>,decltype (ELH_L)::HBase>::type::sameGrid(ELH_L);
    auto VecHisto1_HL = Function::GridExtractorType<std::vector<double>,decltype (ELH_H)::HBase>::type::sameGrid(ELH_H);

    for(size_t i=0;i<VecHisto_LH.Grid.size()-1;++i){
        double e_nd = 0.5*(VecHisto_LH.Grid[i]+VecHisto_LH.Grid[i+1]);
        for(size_t j=0;j<VecHisto_LH.values[i].Grid.size()-1;++j){
            double l_norm = 0.5*(VecHisto_LH.values[i].Grid[j]+VecHisto_LH.values[i].Grid[j+1]);
            VecHisto_LH.values[i].values[j] = HTV_LH(e_nd,l_norm);
        }
    }
    for(size_t i=0;i<VecHisto_HL.Grid.size()-1;++i){
        double e_nd = 0.5*(VecHisto_HL.Grid[i]+VecHisto_HL.Grid[i+1]);
        for(size_t j=0;j<VecHisto_HL.values[i].Grid.size()-1;++j){
            double l_norm = 0.5*(VecHisto_HL.values[i].Grid[j]+VecHisto_HL.values[i].Grid[j+1]);
            VecHisto_HL.values[i].values[j] = HTV_HL(e_nd,l_norm);
        }
    }

    auto VecHisto_LH_mapped = VecHisto_LH;
    VecHisto_LH_mapped.map(HTV_LH);






    auto VecHisto_HL_mapped = VecHisto_HL;
    VecHisto_HL_mapped.map(HTV_HL);


    for(size_t i=0;i<VecHisto1_LH.Grid.size()-1;++i){
        double e_nd = 0.5*(VecHisto1_LH.Grid[i]+VecHisto1_LH.Grid[i+1]);
        for(size_t j=0;j<VecHisto1_LH.values[i].Grid.size()-1;++j){
            double l_norm = 0.5*(VecHisto1_LH.values[i].Grid[j]+VecHisto1_LH.values[i].Grid[j+1]);
            double l_nd = maxLnd(BM,e_nd)*l_norm;
            VecHisto1_LH.values[i].values[j] = HTV1_LH(e_nd,l_nd);
        }
    }
    for(size_t i=0;i<VecHisto1_HL.Grid.size()-1;++i){
        double e_nd = 0.5*(VecHisto1_HL.Grid[i]+VecHisto1_HL.Grid[i+1]);
        for(size_t j=0;j<VecHisto1_HL.values[i].Grid.size()-1;++j){
            double l_norm = 0.5*(VecHisto1_HL.values[i].Grid[j]+VecHisto1_HL.values[i].Grid[j+1]);
            double l_nd = maxLnd(BM,e_nd)*l_norm;
            VecHisto1_HL.values[i].values[j] = HTV1_HL(e_nd,l_nd);
        }
    }


    auto VecHisto1_LH_mapped = VecHisto1_LH;
    VecHisto1_LH_mapped.map([&HTV1_LH,&BM](double e_nd,double l_norm){
        return HTV1_LH(e_nd,l_norm*maxLnd(BM,e_nd));
    });

    auto VecHisto1_HL_mapped = VecHisto1_HL;
    VecHisto1_HL_mapped.map([&HTV1_HL,&BM](double e_nd,double l_norm){
        return HTV1_HL(e_nd,l_norm*maxLnd(BM,e_nd));
    });

    auto Mat0_LH = VecHisto_LH.AllValues();
    auto Mat0_HL = VecHisto_HL.AllValues();





    auto MT_LH = MatrixTranspose(Mat0_LH);
    auto MT_HL = MatrixTranspose(Mat0_HL);

    auto MT_LH_out = vmap([](const auto &V){return vector_sum(V);},Mat0_LH);
    auto MT_HL_out = vmap([](const auto &V){return vector_sum(V);},Mat0_HL);

    auto Mat0_LH_mapped = VecHisto_LH_mapped.AllValues();

    auto Mat1_LH = VecHisto1_LH.AllValues();
    auto Mat1_LH_mapped = VecHisto1_LH_mapped.AllValues();


    SupressFactor(ELH_H,mk,-delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",100000);
    //SupressFactor(ELH_L,mk,-delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",100000);


    double c_capt = ELH_H.summ();
    PVAR(c_capt);
    auto EH_Distrib = ELH_H.AllValues()/c_capt;
    auto EL_Distrib = EH_Distrib*0;

    auto EvapHisto_H = ELH_H;
    EvapHisto_H.map(EvaporationPreHisto_H);
    auto Ev_H = EvapHisto_H.AllValues();

    auto EvapHisto_L = ELH_L;
    EvapHisto_L.map(EvaporationPreHisto_L);
    auto Ev_L = EvapHisto_L.AllValues();





    auto vector_min = [](const auto &Vec){
        auto x = *Vec.begin();
        for(auto it = Vec.begin();it != Vec.end();++it){
            x = std::max(x,*it);
        }
        return x;
    };
    //PVAR(MT_HL_out);
    //PVAR(MT_LH_out);
    //PVAR(Ev_H);
    //PVAR(Ev_L);



    auto E_Histo = MergeHisto(ELH_H);


    Gnuplot hs;
    hs.plotd(E_Histo.toFunction().toString(),"with steps");
    hs.show();

    std::ofstream Mat0_LH_file("D:/Desktop/del_files/matrix/Mat0_LH.matrix",std::ios::binary);
    SaveMatrixBinary(Mat0_LH,Mat0_LH_file);

    std::ofstream Mat0_HL_file("D:/Desktop/del_files/matrix/Mat0_HL.matrix",std::ios::binary);
    SaveMatrixBinary(Mat0_HL,Mat0_HL_file);

    std::ofstream Ev_H_file("D:/Desktop/del_files/matrix/Ev_H.matrix",std::ios::binary);
    SaveVectorBinary(Ev_H,Ev_H_file);

    std::ofstream Ev_L_file("D:/Desktop/del_files/matrix/Ev_L.matrix",std::ios::binary);
    SaveVectorBinary(Ev_L,Ev_L_file);

    std::ofstream H_capt_file("D:/Desktop/del_files/matrix/H_capt.matrix",std::ios::binary);
    SaveVectorBinary(EH_Distrib,H_capt_file);


    dmatrix LLmat,HLmat;
    dmatrix LHmat,HHmat;

    auto E1 = E_Matrix<double>(Ev_L.size());
    double Nu_Max = std::max(mnorm_1(MT_LH),mnorm_1(MT_HL));
    double tau = std::min(0.01,0.1/Nu_Max);

    _R(LLmat,HLmat,LHmat,HHmat) = smart_pow(_T(E1-tau*MatrixDiagonal(MT_LH_out),tau*MT_HL,
                                               tau*MT_LH,E1-tau*MatrixDiagonal(MT_HL_out)),100000,
                                            [](const auto &T1,const auto &T2){
        const auto & A11 = std::get<0>(T1);
        const auto & A12 = std::get<1>(T1);
        const auto & A21 = std::get<2>(T1);
        const auto & A22 = std::get<3>(T1);

        const auto & B11 = std::get<0>(T2);
        const auto & B12 = std::get<1>(T2);
        const auto & B21 = std::get<2>(T2);
        const auto & B22 = std::get<3>(T2);

        return  _T(MatrixMult(A11,B11)+MatrixMult(A12,B21),MatrixMult(A11,B12)+MatrixMult(A12,B22),
                   MatrixMult(A21,B11)+MatrixMult(A22,B21),MatrixMult(A21,B12)+MatrixMult(A22,B22));
    });

    auto LfinV = vmap([](auto x){return x[0];},LLmat);
    auto HfinV = vmap([](auto x){return x[0];},HHmat);

    auto LfinH = ELH_L;
    LfinH.loadIter(LfinV.begin(),LfinV.end());
    auto HfinH = ELH_H;
    HfinH.loadIter(HfinV.begin(),HfinV.end());

    Gnuplot fdistrib;
    fdistrib.plotd(MergeHisto(LfinH).toFunction().toString(),"with steps title \"light\"");
    fdistrib.plotd(MergeHisto(HfinH).toFunction().toString(),"with steps title \"heavy\"");
    fdistrib.show();


    wait();
    */
    return 0;
}
