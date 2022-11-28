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






#define VAR_1
//#define VAR_1
int main(void)
{

    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    double Emin = -BM["phi"][0];
    auto G = [](){return rand()/(RAND_MAX+1.0);};
    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-13);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);

    size_t NE = 10+1;
    size_t NLmax = 10 + 1;
    Function::UniformGrid<double> E_grid(Emin,0,NE);
    EL_Histo<double,Function::UniformGrid<double>,Function::UniformGrid<double>> ELH(E_grid,Function::FunctorM(BM,maxLnd));
    double Lmax = max(ELH.Lmax.values);
    for(size_t i=0;i<ELH.values.size();++i){
        ELH.values[i] = Function::Histogramm<double,Function::UniformGrid<double>>(
                    Function::UniformGrid<double>(0,1,2 + ELH.Lmax.values[i]/Lmax*NLmax));
    }
    //PVAR(ELH.gridStr());

    Function::GridFunction<decltype(ELH),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto(E_grid,Vector(E_grid.size(),[&ELH](size_t i){
                                size_t N = ELH.values[std::min(i,ELH.values.size()-1)].Grid.size();
                                return Function::GridFunction<decltype(ELH),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,std::max(ELH.Lmax.values[i],1e-10),N),
                                std::vector<decltype(ELH)>(N,ELH));
                            }));

    auto EvaporationPreHisto = Function::GridFunctionCreator2<Function::LinearInterpolator>::Create(E_grid,
        Vector(E_grid.size(),[&ELH](size_t i){
            size_t N = ELH.values[std::min(i,ELH.values.size()-1)].Grid.size();
            return Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(
                Function::UniformGrid<double>(0,ELH.Lmax.values[i],N),std::vector<double>(N,0));
        }));


    double mk = 1;
    double delta_mk = 0;//1e-4;

    auto InelasticFactor = dF_H<ELASTIC,PROTON>();
    /*
    struct dF_New{

    };*/
    //InelasticFactor.
    for(size_t i=0;i<ScatterPreHisto.Grid.size();++i){
        for(size_t j = 0;j<ScatterPreHisto.values[i].Grid.size();++j){
            //PVAR(ScatterPreHisto.Grid[i]);
            //PVAR(ScatterPreHisto.values[i].Grid[j]);
            auto TI = CalculateTrajectory(PhiC,ScatterPreHisto.Grid[i],ScatterPreHisto.values[i].Grid[j],10);
            TrajectoryIntegral(G,ScatterPreHisto.Grid[i],ScatterPreHisto.values[i].Grid[j],TI,
                               BM.VescMin(),[](double){return 1.0;},Therm,Vesc,ScatterPreHisto.values[i].values[j],EvaporationPreHisto.values[i].values[j],
                               mk,/*mp,*/delta_mk,ELASTIC,PROTON,PhiFactorSS,400000);
        }
    }

    SupressFactor(ELH,/*mp,*/mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",10000);

    auto ScatterFunction = ScatterPreHisto.Composition([](auto x){return x.AllValues();});
    auto/*histo*/ HMat = Function::GridExtractorType<std::vector<double>,decltype(ELH)::HBase>::type::sameGrid(ELH);
    HMat.map([&ScatterFunction,&ELH](double e,double l){return ScatterFunction(e,l*ELH.Lmax(e));});
    auto Mat_In_Out = HMat.AllValues();
    auto S_matrix = MatrixTranspose(Mat_In_Out)-MatrixDiagonal(Vector(Mat_In_Out.size(),[&Mat_In_Out](size_t i){
                                                                   return vector_sum(Mat_In_Out[i]);
                                                             }));

    auto/*histo*/ EvFunction= Function::GridExtractorType<double,decltype(ELH)::HBase>::type::sameGrid(ELH);
    EvFunction.map([&ELH,&EvaporationPreHisto](double e,double l){return EvaporationPreHisto(e,ELH.Lmax(e)*l);});

    PVAR(EvaporationPreHisto.toString());
    //PVAR(EvFunction.toString());

    auto Evap_matrix = MatrixDiagonal(EvFunction.AllValues());

    auto I_matrix = E_Matrix<double>(S_matrix.size());


    double tau = 1;
    auto R_matrix_scatter = I_matrix + tau*S_matrix;
    auto R_matrix = I_matrix + tau*(S_matrix-Evap_matrix);

    //PVAR(R_matrix);
    //PVAR(R_matrix_scatter);

    PVAR(S_matrix);

    auto Markov_result_matrix = MatrixPow(R_matrix_scatter,10001);

    //PVAR(Markov_result_matrix);

    auto EquilibriumVector = vmap([](auto x){return x[0];},Markov_result_matrix);

    //PVAR(EquilibriumVector);

    EquilibriumVector /= vector_sum(EquilibriumVector);
    double EvapSpeed = vector_sum(dot(Evap_matrix,EquilibriumVector));

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator>
            RealTrajectoryFull(Function::UniformGrid<double>(0,1,101));



    double T_v = tau*(RealTrajectoryFull.Grid.size()-1);
    auto tmpVector = ELH.AllValues()/ELH.summ();
    PVAR(tmpVector);
    RealTrajectoryFull[0] = 1;
    for(size_t i=1;i<RealTrajectoryFull.Grid.size();++i){
        tmpVector = dot(R_matrix,tmpVector);
        RealTrajectoryFull[i] = vector_sum(tmpVector);
        if(i!=0)
            print("del: ",(1 - RealTrajectoryFull[i]/RealTrajectoryFull[i-1]));
    }


    PVAR(tmpVector/vector_sum(tmpVector));
    PVAR(EquilibriumVector);

    Function::UniformGrid<double> Tgrid(0,1,101);

    decltype (RealTrajectoryFull) Real_traj(Tgrid);
    decltype (RealTrajectoryFull) Simple_traj(Tgrid);

    Simple_traj[0] = 1;
    double m_fact = pow(1-EvapSpeed*tau,(RealTrajectoryFull.Grid.size()-1)/(Tgrid.size()-1));
    for(size_t i=0;i<Real_traj.Grid.size();++i){
        Real_traj[i] = RealTrajectoryFull(Real_traj.Grid[i]);
        if(i!=0){
            Simple_traj[i] =  Simple_traj[i-1]*m_fact;

        }
    }



    PVAR(EvapSpeed);

    PVAR(RealTrajectoryFull.toString());
    PVAR(Simple_traj.toString());

    Gnuplot plt;
    plt.plotd(Real_traj.toString(),"with lines title \"Real-traj\"");
    plt.plotd(Simple_traj.toString(),"with lines title \"Simle-traj\"");
    plt.show();

    std::cin.get();



    return 0;
}


