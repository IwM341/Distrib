#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/dens_grid.hpp"

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



int main(int argc,char **argv){

    double Emin = -3;
    double Lmax = 1.13;
    auto maxLnd_l = [Emin,Lmax](double e){return (e-Emin)*Lmax/(-Emin);};
    size_t NE = 30+1;
    size_t NL_max = 30+1;
    auto E_grid = CreateEGrid_sqrt(31,Emin);
    auto N_L_func = [&NL_max,Lmax,&maxLnd_l](double _E){return 2+(NL_max-2)*maxLnd_l(_E)/Lmax;};
    auto Lgridrho = [](double _E,double _L_nd){return 1.0;};
    auto HistoGrid = Create_EL_grid(E_grid,Lgridrho,N_L_func);

    EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>> my_ELH
            (HistoGrid,maxLnd_l);

    PVAR(my_ELH.Lmax.toString());
    PVAR(my_ELH.gridStr());

    my_ELH.save_text(std::ofstream("../tests/my_grid.txt"));
    auto loaded_ELH = decltype(my_ELH)::load_text(std::ifstream("../tests/my_grid.txt"));

    std::ofstream f1("../tests/before.txt");
    std::ofstream f2("../tests/after.txt");
    f1<<my_ELH.Lmax.toString()<<std::endl << "**********" <<my_ELH.gridStr()<<std::endl;
    f2<<loaded_ELH.Lmax.toString()<<std::endl << "**********" <<loaded_ELH.gridStr()<<std::endl;

    return 0;
}
