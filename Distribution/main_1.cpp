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






int main(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    double Emin = -BM["phi"][0];
    auto G = [](){return rand()/(RAND_MAX+0.);};
    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-14);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);

    size_t NE = 30+1;
    size_t NLmax = 30 + 1;
    Function::UniformGrid<double> E_grid(Emin,0,NE);
    EL_Histo<double,Function::UniformGrid<double>,Function::UniformGrid<double>> ELH(E_grid,Function::FunctorM(BM,maxLnd));
    double Lmax = max(ELH.Lmax.values);
    for(size_t i=0;i<ELH.values.size();++i){
        ELH.values[i] = Function::Histogramm<double,Function::UniformGrid<double>>(
                    Function::UniformGrid<double>(0,1,2 + ELH.Lmax.values[i]/Lmax*NLmax));
    }
    PVAR(ELH.gridStr());

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



    double mk = 3;
    double delta_mk = 1e-4;



    for(size_t i=0;i<ScatterPreHisto.Grid.size();++i){
        for(size_t j = 0;j<ScatterPreHisto.values[i].Grid.size();++j){
            auto TI = CalculateTrajectory(PhiC,ScatterPreHisto.Grid[i],ScatterPreHisto.values[i].Grid[j],10);
            TrajectoryIntegral(G,ScatterPreHisto.Grid[i],ScatterPreHisto.values[i].Grid[j],TI,
                               BM.VescMin(),Therm,[](double){return 1.0;},Vesc,ScatterPreHisto.values[i].values[j],EvaporationPreHisto.values[i].values[j],
                               mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,10000);
        }
    }

    SupressFactor(ELH,mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",10000);

    PVAR(vmap([](const auto &x){return x.summ();},ScatterPreHisto.AllValues()));
    PVAR(ScatterPreHisto.gridStr());

    auto MatHisto = Function::GridExtractorType<std::vector<double>,decltype(ScatterPreHisto)>::type::sameGrid(ScatterPreHisto);
    for(size_t i=0;i<ScatterPreHisto.Grid.size();++i){
        for(size_t j=0;j<ScatterPreHisto.values[i].Grid.size();++j){
            MatHisto.values[i].values[j] = ScatterPreHisto.values[i].values[j].AllValues();
        }
    }
    auto HMat = Function::GridExtractorType<std::vector<double>,decltype(ELH)::HBase>::type::sameGrid(ELH);
    for(size_t i=0;i<HMat.Grid.size()-1;++i){
        double E_av = 0.5*(HMat.Grid[i]+HMat.Grid[i+1]);
        for(size_t j=0;j<HMat.values[i].Grid.size()-1;++j){
            double L_av = 0.5*(HMat.values[i].Grid[j]+HMat.values[i].Grid[j+1])*ELH.Lmax.values[i];
            HMat.values[i].values[j] = MatHisto(E_av,L_av);
        }
    }
    auto SmatT = HMat.AllValues();
    auto Rmat = Rmatrix(SmatT,10.0);


    Gnuplot distrib;
    distrib.show_cmd = "splot";
    distrib.plotd(ELH.toFunction().toString(),"with lines title \"distrib\"");


    auto ResDistrib = ELH;

    auto ResVector = dot(MatrixPow(Rmat,1000),ELH.AllValues());
    ResDistrib.loadIter(ResVector.begin(),ResVector.end());

    PVAR(ELH.summ());
    PVAR(vector_sum(ResVector));

    distrib.plotd(ResDistrib.toFunction().toString(),"with lines title \"final distrib\"");



    distrib.show();
    std::cin.get();



    //PVAR(SmatT);
    //PVAR(EvaporationPreHisto.toString());
    //PVAR(ScatterPreHisto.gridStr());
    return 0;
}
