#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
#include "func/el_histo.hpp"

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


template <typename HistoType>
struct HBaseType;

template <typename ...Args>
struct HBaseType<Function::Histogramm<Args...>>{
    typedef Function::Histogramm<Args...> type;
};


template <typename ...Args>
struct HBaseType<EL_Histo<Args...>>{
    typedef typename EL_Histo<Args...>::HBase type;
};
template<typename Histo,typename...Args>
std::vector<std::vector<double>> Func_Of_Histo_To_Matrix(const Function::GridFunction<Args...> &PreHisto,const Histo& GridH){
    auto VectorFunc = PreHisto.Composition([](const auto &x){return x.AllValues();});
    auto HistoMat = Function::GridExtractorType<std::vector<double>,typename HBaseType<Histo>::type>::type::sameGrid(GridH);
    HistoMat.map(VectorFunc);
    return HistoMat.AllValues();
}

template<typename Histo,typename...Args>
std::vector<double> Func_Histo_To_Vector(const Function::GridFunction<Args...> &PreHisto,const Histo& GridH){
    auto HistoMat = Function::GridExtractorType<double,typename HBaseType<Histo>::type>::type::sameGrid(GridH);
    HistoMat.map(PreHisto);
    return HistoMat.AllValues();
}

int main0(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    double Emin = -BM["phi"][0];




    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-13);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);

    size_t NE = 30+1;
    size_t NLmax = 30 + 1;
    Function::UniformGrid<double> E_grid(Emin,0,NE);
    EL_Histo<double,Function::UniformGrid<double>,Function::UniformGrid<double>> ELH_L(E_grid,Function::FunctorM(BM,maxLnd));

    double Lmax = max(ELH_L.Lmax.values);
    for(size_t i=0;i<ELH_L.values.size();++i){
        ELH_L.values[i] = Function::Histogramm<double,Function::UniformGrid<double>>(
                    Function::UniformGrid<double>(0,1,2 + ELH_L.Lmax.values[i]/Lmax*NLmax));
    }
    auto ELH_H = ELH_L;
    //PVAR(ELH_L.gridStr());

    Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto_LH(E_grid,Vector(E_grid.size(),[&ELH_L](size_t i){
                                size_t N = ELH_L.values[std::min(i,ELH_L.values.size()-1)].Grid.size();
                                return Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,1.0/*std::max(ELH_L.Lmax.values[i],1e-10)*/,N),
                                std::vector<decltype(ELH_L)>(N,ELH_L));
                            }));
    auto ScatterPreHisto_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_Ion_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_AntiIon_LH = ScatterPreHisto_LH;

    auto EvaporationPreHisto_L = ScatterPreHisto_LH.Composition([](auto x)->double{return 0.0;});
            /*Function::GridFunctionCreator2<Function::LinearInterpolator>::Create(E_grid,
        Vector(E_grid.size(),[&ELH_L](size_t i){
            size_t N = ELH_L.values[std::min(i,ELH_L.values.size()-1)].Grid.size();
            return Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(
                Function::UniformGrid<double>(0,ELH_L.Lmax.values[i],N),std::vector<double>(N,0));
        }));*/
    auto EvaporationPreHisto_H = EvaporationPreHisto_L;



    double mk = 5;
    double delta_mk = 3e-5;


    double f_p_av = 4e-7;
    decltype(Therm) F_e_r(R,Vector(R.size(),[&BM,&Therm,f_p_av](size_t i){
        double rho = BM["RhoND"][i];
        double s = pow(2*M_PI*Therm[i]/me,1.5)*exp(-Rd/Therm[i]);
        return f_p_av*rho*Ion_Deg(f_p_av*rho,s);
    }));


    for(size_t i=0;i<ScatterPreHisto_LH.Grid.size();++i){
        double e_nd=  ScatterPreHisto_LH.Grid[i];
        for(size_t j = 0;j<ScatterPreHisto_LH.values[i].Grid.size();++j){           
            double l_nd = ScatterPreHisto_LH.values[i].Grid[j]*ELH_H.Lmax[i];

            auto TI = CalculateTrajectory(PhiC,e_nd,l_nd,10);

            TrajectoryIntegral(G,e_nd,l_nd,TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_LH.values[i].values[j],EvaporationPreHisto_L.values[i].values[j],
                               mk,-delta_mk,ELASTIC,PROTON,PhiFactorSS,100000);
            TrajectoryIntegral(G,e_nd,l_nd,TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_HL.values[i].values[j],EvaporationPreHisto_H.values[i].values[j],
                               mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,100000);
            TrajectoryIntegral(G,e_nd,l_nd,TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_Ion_HL.values[i].values[j],EvaporationPreHisto_H.values[i].values[j],
                               mk,delta_mk,IONIZATION,PROTON,PhiFactorSS,100000);

            TrajectoryIntegral_AntiIonization(G,e_nd,l_nd,TI,[](double){return 1.0;},F_e_r,Therm,Vesc,BM.VescMin(),
                                ScatterPreHisto_AntiIon_LH.values[i].values[j],EvaporationPreHisto_L.values[i].values[j],PROTON,
                                mk,mp,-delta_mk,PhiFactorSS,100000);
        }
    }

    SupressFactor(ELH_H,mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",100000);
    SupressFactor(ELH_L,mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",100000);

    //PVAR(vmap([](const auto &x){return x.summ();},ScatterPreHisto_LH.AllValues()));
    //PVAR(ScatterPreHisto_LH.gridStr());

    auto MatHisto_LH = Function::GridExtractorType<std::vector<double>,decltype(ScatterPreHisto_LH)>::type::sameGrid(ScatterPreHisto_LH);
    auto MatHisto_HL = MatHisto_LH;
    auto MatHisto_Ion_HL = MatHisto_LH;

    for(size_t i=0;i<ScatterPreHisto_LH.Grid.size();++i){
        for(size_t j=0;j<ScatterPreHisto_LH.values[i].Grid.size();++j){
            MatHisto_LH.values[i].values[j] = ScatterPreHisto_LH.values[i].values[j].AllValues();
            MatHisto_HL.values[i].values[j] = ScatterPreHisto_HL.values[i].values[j].AllValues();
            MatHisto_Ion_HL.values[i].values[j] = ScatterPreHisto_Ion_HL.values[i].values[j].AllValues();
        }
    }
    auto HMat_LH = Function::GridExtractorType<std::vector<double>,decltype(ELH_L)::HBase>::type::sameGrid(ELH_L);
    auto HMat_HL = HMat_LH;
    auto HMat_Ion_HL = HMat_LH;

    auto HEvMat_H = decltype(ELH_L)::HBase::sameGrid(ELH_L);
    auto HEvMat_L = HEvMat_H;

    for(size_t i=0;i<HMat_LH.Grid.size()-1;++i){
        double E_av = 0.5*(HMat_LH.Grid[i]+HMat_LH.Grid[i+1]);
        for(size_t j=0;j<HMat_LH.values[i].Grid.size()-1;++j){
            double L_av = 0.5*(HMat_LH.values[i].Grid[j]+HMat_LH.values[i].Grid[j+1])*ELH_L.Lmax.values[i];
            HMat_LH.values[i].values[j] = MatHisto_LH(E_av,L_av);
            HMat_HL.values[i].values[j] = MatHisto_HL(E_av,L_av);
            HMat_Ion_HL.values[i].values[j] = MatHisto_Ion_HL(E_av,L_av);
            HEvMat_H.values[i].values[j] = EvaporationPreHisto_H(E_av,L_av);
            HEvMat_L.values[i].values[j] = EvaporationPreHisto_L(E_av,L_av);
        }
    }
    auto MT_LH_T = HMat_LH.AllValues();
    auto MT_HL_T = HMat_HL.AllValues();
    auto MT_Ion_HL_T = HMat_Ion_HL.AllValues();




    auto MT_LH_out = vmap([](const auto &V){return vector_sum(V);},MT_LH_T);
    auto MT_HL_out = vmap([](const auto &V){return vector_sum(V);},MT_HL_T);
    auto MT_Ion_HL_out = vmap([](const auto &V){return vector_sum(V);},MT_Ion_HL_T);

    PVAR(MT_LH_out);
    PVAR(MT_HL_out);
    PVAR(MT_Ion_HL_out);

    auto Ev_H = HEvMat_H.AllValues();
    auto Ev_L = HEvMat_L.AllValues();

    auto MT_LH = MatrixTranspose(MT_LH_T);
    auto MT_HL = MatrixTranspose(MT_HL_T);
    auto MT_Ion_HL = MatrixTranspose(MT_Ion_HL_T);

    auto R_Func = [&MT_LH,&MT_HL,&Ev_H,&Ev_L,&MT_LH_out,&MT_HL_out](auto L,auto H,double h){
        return _T(L + h*(dot(MT_HL,H)-Ev_L*L-MT_LH_out*L),H + h*(dot(MT_LH,L) - Ev_H*H-MT_HL_out*H));
    };

    auto R_Func_Ion = [&MT_LH,&MT_HL,&MT_Ion_HL,&Ev_H,&Ev_L,&MT_LH_out,&MT_HL_out,&MT_Ion_HL_out](auto L,auto H,double h){
        return _T(L + h*(dot(MT_HL,H)+dot(MT_Ion_HL,H)-Ev_L*L-MT_LH_out*L),
                  H + h*(dot(MT_LH,L) - Ev_H*H-MT_HL_out*H-MT_Ion_HL_out*H));
    };

    double n_fact = exp(-2*delta_mk/(mk*U0*U0));
    PVAR(n_fact);



    Function::UniformGrid<double> T_grid(0,100,101);
    Function::GridFunction<double,decltype(T_grid),Function::LinearInterpolator> N_L(T_grid);
    decltype(N_L) N_H(T_grid);

    decltype(N_L) N_L_Ion(T_grid);
    decltype(N_L) N_H_Ion(T_grid);

    auto H_D = ELH_H.AllValues()/ELH_H.summ();
    auto H_D_Ion = H_D;
    auto L_D = ELH_L.AllValues()/ELH_L.summ()*n_fact;
    auto L_D_Ion = L_D;

    N_L[0] = vector_sum(L_D);
    N_H[0] = vector_sum(H_D);
    N_L_Ion[0] = vector_sum(H_D_Ion);
    N_H_Ion[0] = vector_sum(L_D_Ion);

    for(size_t i=1;i<T_grid.size();++i){
        _R(L_D,H_D) = R_Func(L_D,H_D,0.1);
        _R(L_D_Ion,H_D_Ion) = R_Func_Ion(L_D_Ion,H_D_Ion,0.1);
        N_L[i] = vector_sum(L_D);
        N_H[i] = vector_sum(H_D);

        N_L_Ion[i] = vector_sum(H_D_Ion);
        N_H_Ion[i] = vector_sum(L_D_Ion);
    }

    double eh = vector_sum(Ev_H*H_D)/vector_sum(H_D);
    double el = vector_sum(Ev_L*L_D)/vector_sum(L_D);

    double h_m = vector_sum(MT_HL_out*H_D)/vector_sum(H_D);
    double h_m_ion = vector_sum(MT_Ion_HL_out*H_D)/vector_sum(H_D);
    double l_m = vector_sum(MT_LH_out*L_D)/vector_sum(L_D);

    double hl = vector_sum(dot(MT_HL,H_D))/vector_sum(H_D);
    double hl_ion = vector_sum(dot(MT_Ion_HL,H_D))/vector_sum(H_D);
    double lh = vector_sum(dot(MT_LH,L_D))/vector_sum(L_D);

    std::cout << std::scientific;
    std::cout << "h\' = " << "h + " << lh <<"*l" <<" - " << h_m << "*h"<<" - " << h_m_ion << "*h(ion)" << " - " << eh << "*h" <<std::endl;
    std::cout << "l\' = " << "l + " << hl <<"*h" <<" + " << hl_ion <<"*h(ion)"<<" - " << l_m << "*l" << " - " << el << "*l" <<std::endl;

    Gnuplot plt;
    plt.plotd(N_L.toString(),"with lines title \"light\"");
    plt.plotd(N_H.toString(),"with lines title \"heavy\"");

    plt.plotd(N_L_Ion.toString(),"with lines title \"light with ion\"");
    plt.plotd(N_H_Ion.toString(),"with lines title \"heavy with ion\"");

    plt.show();
    std::cin.get();
    return 0;
}

int main1(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    double Emin = -BM["phi"][0];




    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-13);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);

    size_t NE = 30+1;
    size_t NLmax = 30 + 1;
    Function::UniformGrid<double> E_grid(Emin,0,NE);
    EL_Histo<double,Function::UniformGrid<double>,Function::UniformGrid<double>> ELH_L(E_grid,Function::FunctorM(BM,maxLnd));

    double Lmax = max(ELH_L.Lmax.values);
    for(size_t i=0;i<ELH_L.values.size();++i){
        ELH_L.values[i] = Function::Histogramm<double,Function::UniformGrid<double>>(
                    Function::UniformGrid<double>(0,1,2 + ELH_L.Lmax.values[i]/Lmax*NLmax));
    }
    auto ELH_H = ELH_L;
    //PVAR(ELH_L.gridStr());

    Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto_LH(E_grid,Vector(E_grid.size(),[&ELH_L](size_t i){
                                size_t N = ELH_L.values[std::min(i,ELH_L.values.size()-1)].Grid.size();
                                return Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,1.0,N),
                                std::vector<decltype(ELH_L)>(N,ELH_L));
                            }));
    auto ScatterPreHisto_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_Ion_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_AntiIon_LH = ScatterPreHisto_LH;

    auto EvaporationPreHisto_L = ScatterPreHisto_LH.Composition([](auto x)->double{return 0.0;});
    auto EvaporationPreHisto_H = EvaporationPreHisto_L;



    double mk = 0.5;
    double delta_mk = 3e-6;


    decltype(Therm) H_N(R,BM["RhoND"]*BM["H1"]);
    double f_p_av = 4e-7;
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

    auto Atom_N =  H_N*atom_degree;
    PVAR(Atom_N.toString());

    Gnuplot ion_plot;
    ion_plot.plotd(atom_degree.toString(),"with lines title \"atom percent\"");
    ion_plot.plotd(ion_degree.toString(),"with lines title\"eletron percent\"");
    ion_plot.command("set logscale y");
    ion_plot.show();

    for(size_t i=0;i<ScatterPreHisto_LH.Grid.size();++i){
        double e_nd=  ScatterPreHisto_LH.Grid[i];
        for(size_t j = 0;j<ScatterPreHisto_LH.values[i].Grid.size();++j){
            double l_nd = ScatterPreHisto_LH.values[i].Grid[j]*ELH_H.Lmax[i];

            auto TI = CalculateTrajectory(PhiC,e_nd,l_nd,10);

            TrajectoryIntegral1_Opt(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               ScatterPreHisto_LH.values[i].values[j],
                               EvaporationPreHisto_L.values[i].values[j],
                               mk,mp,-delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);
            TrajectoryIntegral1_Opt(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               ScatterPreHisto_HL.values[i].values[j],
                               EvaporationPreHisto_H.values[i].values[j],
                               mk,mp,delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);
/*
            TrajectoryIntegral1_Opt(G,e_nd,l_nd,TI,
                               BM.VescMin(),Atom_N,Therm,Vesc,ScatterPreHisto_Ion_HL.values[i].values[j],EvaporationPreHisto_H.values[i].values[j],
                               mk,mp,delta_mk,dF_H<IONIZATION,PROTON,decltype(G)>(G),PhiFactorSS,100000);

            TrajectoryIntegral_AntiIonization(G,e_nd,l_nd,TI,H_N,F_e_r,Therm,Vesc,BM.VescMin(),
                                ScatterPreHisto_AntiIon_LH.values[i].values[j],EvaporationPreHisto_L.values[i].values[j],PROTON,
                                mk,mp,-delta_mk,PhiFactorSS,100000);
*/
        }
    }

    SupressFactor(ELH_H,mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",100000);
    SupressFactor(ELH_L,mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",100000);

    //PVAR(vmap([](const auto &x){return x.summ();},ScatterPreHisto_LH.AllValues()));
    //PVAR(ScatterPreHisto_LH.gridStr());


    auto MT_LH_T = Func_Of_Histo_To_Matrix(ScatterPreHisto_LH,ELH_L);
    auto MT_HL_T = Func_Of_Histo_To_Matrix(ScatterPreHisto_HL,ELH_L);
    auto MT_Ion_HL_T =  Func_Of_Histo_To_Matrix(ScatterPreHisto_Ion_HL,ELH_L);
    auto MT_Anti_Ion_LH_T =  Func_Of_Histo_To_Matrix(ScatterPreHisto_AntiIon_LH,ELH_L);



    auto MT_LH_out = vmap([](const auto &V){return vector_sum(V);},MT_LH_T);
    auto MT_HL_out = vmap([](const auto &V){return vector_sum(V);},MT_HL_T);
    auto MT_Ion_HL_out = vmap([](const auto &V){return vector_sum(V);},MT_Ion_HL_T);
    auto MT_AI_LH_Out = vmap([](const auto &V){return vector_sum(V);},MT_Anti_Ion_LH_T);

    PVAR(MT_LH_out);
    PVAR(MT_HL_out);
    PVAR(MT_Ion_HL_out);
    PVAR(MT_AI_LH_Out);

    auto Ev_H = Func_Histo_To_Vector(EvaporationPreHisto_H,ELH_H);
    auto Ev_L = Func_Histo_To_Vector(EvaporationPreHisto_L,ELH_L);

    auto MT_LH = MatrixTranspose(MT_LH_T);
    auto MT_AI_LH = MatrixTranspose(MT_Anti_Ion_LH_T);

    auto MT_HL = MatrixTranspose(MT_HL_T);
    auto MT_Ion_HL = MatrixTranspose(MT_Ion_HL_T);

    auto R_Func = [&MT_LH,&MT_HL,&Ev_H,&Ev_L,&MT_LH_out,&MT_HL_out](auto L,auto H,double h){
        return _T(L + h*(dot(MT_HL,H)-Ev_L*L-MT_LH_out*L),H + h*(dot(MT_LH,L) - Ev_H*H-MT_HL_out*H));
    };

    auto R_Func_Ion = [&MT_LH,&MT_HL,&MT_Ion_HL,&Ev_H,&Ev_L,&MT_LH_out,&MT_HL_out,&MT_Ion_HL_out](auto L,auto H,double h){
        return _T(L + h*(dot(MT_HL,H)+dot(MT_Ion_HL,H)-Ev_L*L-MT_LH_out*L),
                  H + h*(dot(MT_LH,L) - Ev_H*H-MT_HL_out*H-MT_Ion_HL_out*H));
    };

    auto R_Func_Anti_Ion = [&MT_LH,&MT_HL,&MT_Ion_HL,&Ev_H,&Ev_L,
            &MT_LH_out,&MT_HL_out,&MT_Ion_HL_out,&MT_AI_LH,&MT_AI_LH_Out]
            (auto L,auto H,double h){
        return _T(L + h*(dot(MT_HL,H)+dot(MT_Ion_HL,H)-MT_AI_LH_Out*L-Ev_L*L-MT_LH_out*L),
                  H + h*(dot(MT_LH,L) + dot(MT_AI_LH,L)- Ev_H*H-MT_HL_out*H-MT_Ion_HL_out*H));
    };

    double n_fact = exp(-2*delta_mk/(mk*U0*U0));
    PVAR(n_fact);



    Function::UniformGrid<double> T_grid(0,1,1001);
    Function::GridFunction<double,decltype(T_grid),Function::LinearInterpolator> N_L(T_grid);
    decltype(N_L) N_H(T_grid);

    decltype(N_L) N_L_Ion(T_grid);
    decltype(N_L) N_H_Ion(T_grid);

    decltype(N_L) N_L_AI(T_grid);
    decltype(N_L) N_H_AI(T_grid);

    auto H_D = ELH_H.AllValues()/ELH_H.summ();
    auto H_D_Ion = H_D;
    auto H_D_AI = H_D;

    auto L_D = ELH_L.AllValues()/ELH_L.summ()*n_fact;
    auto L_D_Ion = L_D;
    auto L_D_AI = L_D;

    N_L[0] = vector_sum(L_D);
    N_H[0] = vector_sum(H_D);

    N_L_Ion[0] = vector_sum(L_D_Ion);
    N_H_Ion[0] = vector_sum(H_D_Ion);

    N_L_AI[0] = vector_sum(L_D_AI);
    N_H_AI[0] = vector_sum(H_D_AI);

    for(size_t i=1;i<T_grid.size();++i){
        _R(L_D,H_D) = R_Func(L_D,H_D,0.01);
        _R(L_D_Ion,H_D_Ion) = R_Func_Ion(L_D_Ion,H_D_Ion,0.01);

        _R(L_D_AI,H_D_AI) = R_Func_Anti_Ion(L_D_AI,H_D_AI,0.01);

        N_L[i] = vector_sum(L_D);
        N_H[i] = vector_sum(H_D);

        N_L_Ion[i] = vector_sum(L_D_Ion);
        N_H_Ion[i] = vector_sum(H_D_Ion);

        N_L_AI[i] = vector_sum(L_D_AI);
        N_H_AI[i] = vector_sum(H_D_AI);
    }

    double eh = vector_sum(Ev_H*H_D)/vector_sum(H_D);
    double el = vector_sum(Ev_L*L_D)/vector_sum(L_D);

    double h_m = vector_sum(MT_HL_out*H_D)/vector_sum(H_D);
    double h_m_ion = vector_sum(MT_Ion_HL_out*H_D)/vector_sum(H_D);
    double l_m = vector_sum(MT_LH_out*L_D)/vector_sum(L_D);

    double hl = vector_sum(dot(MT_HL,H_D))/vector_sum(H_D);



    double hl_ion = vector_sum(dot(MT_Ion_HL,H_D))/vector_sum(H_D);
    double lh = vector_sum(dot(MT_LH,L_D))/vector_sum(L_D);

    double lh_ai = vector_sum(dot(MT_AI_LH,L_D))/vector_sum(L_D);

    std::cout << std::scientific;
    std::cout << "h\' = " << "h + " << lh <<"*l" <<" - " <<
                 h_m << "*h"<<" - " << h_m_ion << "*h(ion)" << " + " << lh_ai << "*l(ai)" <<
                 " - " << eh << "*h" <<std::endl;
    std::cout << "l\' = " << "l + " << hl <<"*h" <<" + " << hl_ion <<"*h(ion)"<<
                 " - " << l_m << "*l" <<  " - " << lh_ai << "*l(ai)" << " - " << el << "*l" <<std::endl;

    Gnuplot plt;
    plt.plotd(N_L.toString(),"with lines title \"light\"");
    plt.plotd(N_H.toString(),"with lines title \"heavy\"");

    plt.plotd(N_L_Ion.toString(),"with lines title \"light with ion\"");
    plt.plotd(N_H_Ion.toString(),"with lines title \"heavy with ion\"");

    plt.plotd(N_L_AI.toString(),"with lines title \"light with antiion\"");
    plt.plotd(N_H_AI.toString(),"with lines title \"heavy with antiion\"");

    plt.show();
    std::cin.get();
    return 0;
}

int main3(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");


    double Emin = -BM["phi"][0];
    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-13);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);

    size_t NE = 30+1;
    size_t NLmax = 30 + 1;
    Function::UniformGrid<double> E_grid(Emin,0,NE);
    EL_Histo<double,Function::UniformGrid<double>,Function::UniformGrid<double>> ELH_L(E_grid,Function::FunctorM(BM,maxLnd));

    double Lmax = max(ELH_L.Lmax.values);
    for(size_t i=0;i<ELH_L.values.size();++i){
        ELH_L.values[i] = Function::Histogramm<double,Function::UniformGrid<double>>(
                    Function::UniformGrid<double>(0,1,2 + ELH_L.Lmax.values[i]/Lmax*NLmax));
    }
    auto ELH_H = ELH_L;
    //PVAR(ELH_L.gridStr());

    Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto_LH(E_grid,Vector(E_grid.size(),[&ELH_L](size_t i){
                                size_t N = ELH_L.values[std::min(i,ELH_L.values.size()-1)].Grid.size();
                                return Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,std::max(ELH_L.Lmax.values[i],1e-10),N),
                                std::vector<decltype(ELH_L)>(N,ELH_L));
                            }));
    auto ScatterPreHisto_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_Ion_HL = ScatterPreHisto_LH;

    auto ScatterPreHisto_Ion_HL_1 = ScatterPreHisto_LH;
    auto ScatterPreHisto_LH_1 = ScatterPreHisto_LH;

    auto EvaporationPreHisto_L = Function::GridFunctionCreator2<Function::LinearInterpolator>::Create(E_grid,
        Vector(E_grid.size(),[&ELH_L](size_t i){
            size_t N = ELH_L.values[std::min(i,ELH_L.values.size()-1)].Grid.size();
            return Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(
                Function::UniformGrid<double>(0,ELH_L.Lmax.values[i],N),std::vector<double>(N,0));
        }));
    auto EvaporationPreHisto_H = EvaporationPreHisto_L;

    auto EvaporationPreHisto_H_1 = EvaporationPreHisto_L;
    auto EvaporationPreHisto_L_1 = EvaporationPreHisto_L;


    double mk = 2;
    double delta_mk = 1e-5;




    for(size_t i=0;i<ScatterPreHisto_LH.Grid.size();++i){
        for(size_t j = 0;j<ScatterPreHisto_LH.values[i].Grid.size();++j){
            auto TI = CalculateTrajectory(PhiC,ScatterPreHisto_LH.Grid[i],ScatterPreHisto_LH.values[i].Grid[j],10);

            print("E,L=",ScatterPreHisto_LH.Grid[i],", ",ScatterPreHisto_LH.values[i].Grid[j]);
            TrajectoryIntegral(G,ScatterPreHisto_LH.Grid[i],ScatterPreHisto_LH.values[i].Grid[j],TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_LH.values[i].values[j],EvaporationPreHisto_L.values[i].values[j],
                               mk,-delta_mk,ELASTIC,PROTON,PhiFactorSS,100000);

            TrajectoryIntegral1_Opt(G,ScatterPreHisto_LH_1.Grid[i],ScatterPreHisto_LH.values[i].Grid[j],TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_LH_1.values[i].values[j],EvaporationPreHisto_L_1.values[i].values[j],
                               mk,mp,-delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);

            TrajectoryIntegral(G,ScatterPreHisto_HL.Grid[i],ScatterPreHisto_HL.values[i].Grid[j],TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_HL.values[i].values[j],EvaporationPreHisto_H.values[i].values[j],
                               mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,10000);
            TrajectoryIntegral(G,ScatterPreHisto_HL.Grid[i],ScatterPreHisto_HL.values[i].Grid[j],TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto_Ion_HL.values[i].values[j],EvaporationPreHisto_H.values[i].values[j],
                               mk,delta_mk,IONIZATION,PROTON,PhiFactorSS,10000);

        }
    }

    PVAR(EvaporationPreHisto_L.toString());
    PVAR(EvaporationPreHisto_L_1.toString());




    auto MatHisto_LH = ScatterPreHisto_LH.Composition([](auto x){return x.AllValues();});
    auto MatHisto_LH_Opt = ScatterPreHisto_LH_1.Composition([](auto x){return x.AllValues();});


    for(size_t i=0;i<10;++i){
        double E = Emin*G();
        double L = ELH_H.Lmax(Emin)*G();

        print("(E,L) = (", E,", ",L,")");
        PVAR(MatHisto_LH(E,L));
        PVAR(MatHisto_LH_Opt(E,L));
    }
    return 0;
}

int main(){
    return main1();
}

