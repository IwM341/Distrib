#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
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

int main(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    auto & phi = BM["phi"];
    //auto phi = Vector(BM["phi"].size(),[N = BM["phi"].size()](size_t i){return 1.5-0.5*pow(i/(N-1.0),2);});



    auto G = [](){return rand()/(RAND_MAX+1.0);};


    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-13);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);



    double Emin = -phi[0];
    double Lmax = maxLnd(BM,0);

    size_t NE = 30+1;
    size_t NL_max = 30+1;
    auto E_grid = CreateEGrid_sqrt(31,Emin);
    auto N_L_func = [&NL_max,Lmax,&BM](double _E){return 2+(NL_max-2)*maxLnd(BM,_E)/Lmax;};
    auto Lgridrho = [](double _E,double _L_nd){return 1.0;};
    auto HistoGrid = Create_EL_grid(E_grid,Lgridrho,N_L_func);


    PVAR(HistoGrid.gridStr());
    exit(0);

    size_t NLmax = 20 + 1;




    PVAR(E_grid);

    EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>> ELH_L
            (HistoGrid,Function::FunctorM(BM,maxLnd));

    auto ELH_H = ELH_L;
    //PVAR(ELH_L.gridStr());


    Function::UniformGrid<double> FEgrid(Emin,0,21);
    size_t NLmaxF = 21;
    Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto_LH(FEgrid,Vector(FEgrid.size(),[&ELH_L,&FEgrid,&BM,Lmax,NLmaxF](size_t i){
                                size_t N = 2 + NLmaxF*maxLnd(BM,FEgrid[i])/Lmax;
                                return Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,1.0,N),
                                std::vector<decltype(ELH_L)>(N,ELH_L));
                            }));

    Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto1_LH(FEgrid,Vector(FEgrid.size(),[&ELH_L,&FEgrid,&BM,Lmax,NLmaxF](size_t i){
                                double L = maxLnd(BM,FEgrid[i]);
                                size_t N = 2 + NLmaxF*L/Lmax;
                                return Function::GridFunction<decltype(ELH_L),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,L,N),
                                std::vector<decltype(ELH_L)>(N,ELH_L));
                            }));

    auto ScatterPreHisto_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_Ion_HL = ScatterPreHisto_LH;
    auto ScatterPreHisto_AntiIon_LH = ScatterPreHisto_LH;

    auto ScatterPreHisto1_HL = ScatterPreHisto1_LH;


    auto EvaporationPreHisto_L = ScatterPreHisto_LH.Composition([](auto x)->double{return 0.0;});
            /*Function::GridFunctionCreator2<Function::LinearInterpolator>::Create(E_grid,
        Vector(E_grid.size(),[&ELH_L](size_t i){
            size_t N = ELH_L.values[std::min(i,ELH_L.values.size()-1)].Grid.size();
            return Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(
                Function::UniformGrid<double>(0,ELH_L.Lmax.values[i],N),std::vector<double>(N,0));
        }));*/
    auto EvaporationPreHisto_H = EvaporationPreHisto_L;



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
    //print("end\t","l_nd\t","T_in\t","T_out\t","T_all\t","T_all_th");
    for(size_t i=0;i<ScatterPreHisto_LH.Grid.size();++i){

        double e_nd=  ScatterPreHisto_LH.Grid[i];
        //PVAR(maxLnd(BM,e_nd));
        //PVAR(e_nd);
        for(size_t j = 0;j<ScatterPreHisto_LH.values[i].Grid.size();++j){
            double l_nd = ScatterPreHisto_LH.values[i].Grid[j]*maxLnd(BM,e_nd);//ELH_H.Lmax(e_nd);
            //PVAR(l_nd);
            auto TI = CalculateTrajectory(PhiC,e_nd,l_nd,10);
            //print(e_nd,"\t",l_nd,"\t",TI.T_in,"\t",TI.T_out,"\t",TI.T_in+TI.T_out,"\t",M_PI/2*pow(-e_nd,-1.5));
            //PVAR(TI.Trajectory.Grid);
            //PVAR(TI.Trajectory.values);
            /*
            if(e_nd>-0.5 and l_nd/ELH_H.Lmax[i] > 0.5){
                Gnuplot gp;
                gp.plotd(TI.Trajectory.toString());
                gp.show();
                std::cin.get();
            }
            */
            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               ScatterPreHisto_LH.values[i].values[j],
                               EvaporationPreHisto_L.values[i].values[j],
                               mk,mp,-delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);

            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               ScatterPreHisto_HL.values[i].values[j],
                               EvaporationPreHisto_H.values[i].values[j],
                               mk,mp,delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);




        }
    }


    for(size_t i=0;i<ScatterPreHisto_LH.Grid.size();++i){

        double e_nd=  ScatterPreHisto1_LH.Grid[i];
        //PVAR(maxLnd(BM,e_nd));
        //PVAR(e_nd);
        for(size_t j = 0;j<ScatterPreHisto_LH.values[i].Grid.size();++j){
            double l_nd = ScatterPreHisto1_LH.values[i].Grid[j];
            //PVAR(l_nd);
            auto TI = CalculateTrajectory(PhiC,e_nd,l_nd,10);
            double x;
            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               ScatterPreHisto1_LH.values[i].values[j],
                               x,
                               mk,mp,-delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);

            TrajectoryIntegral1(G,e_nd,l_nd,TI,
                               BM.VescMin(),H_N,Therm,Vesc,
                               ScatterPreHisto1_HL.values[i].values[j],
                               x,
                               mk,mp,delta_mk,dF_H<ELASTIC,PROTON>(),PhiFactorSS,100000);
        }
    }

    //PVAR(ScatterPreHisto_LH.Composition([](auto x){return x.summ();}).toString());
    //PVAR(ScatterPreHisto1_LH.Composition([](auto x){return x.summ();}).toString());

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
    return 0;
}
