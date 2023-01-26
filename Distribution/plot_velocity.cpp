#include <iostream>


#include "functions.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/arg_parser.hpp"
#include "func/dens_grid.hpp"
#include "func/matrix_functions.hpp"


int main(int argc,char**argv){
    double mU0 = 0.73e-3;
    double mUeff = 1.43e-3;
    double mD0 = 0.53e-3;
    double Vesc = 2e-3;

    std::string S;


    while(S != "q"){
        try {
            mUeff = std::stod(S);
        }  catch (std::exception & e) {
            std::cout << "can't understand" << std::endl;
        }
        Function::Histogramm<double,Function::UniformGrid<double>> HV(
                        Function::UniformGrid<double>(0,2*(mD0+mU0),101));
        auto G = [](){return rand()/(RAND_MAX+1.0);};
        auto G1 = [&G](){return (rand()+G())/(RAND_MAX+1.0);};
        const size_t N = 10000000;
        for(size_t i=0;i<N;++i){
            auto[v,dens] = Velocity(G1,Vesc,mD0,mU0);
            double u = sqrt(v.quad()-Vesc*Vesc);
            HV.putValue(dens/N,u);
        }
        std::cout << HV.summ() <<std::endl;

        const auto F = HV.toFunction();
        auto F1 = F;
        F1.map([=](double u){
                if(u==0){
                    return 0.0;
                }
                return sqrt(u*u+Vesc*Vesc)/(sqrt(2*M_PI)*mU0*mD0)*(exp(-0.5*(u-mU0)*(u-mU0)/(mD0*mD0))-
                                                 exp(-0.5*(u+mU0)*(u+mU0)/(mD0*mD0)));
            });
        Gnuplot gp("D:\\Soft\\gnuplot\\bin\\gnuplot.exe");
        gp.plotd(F.toString(),"with lines title \"histo\"");
        gp.plotd(F1.toString(),"with lines title \"true\"");
        gp.show();
        std::cin>>S;
    }
	return 0;
}
