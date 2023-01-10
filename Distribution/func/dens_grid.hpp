#ifndef DENSE_GRID_H
#define DENSE_GRID_H
#include <utils>

template <typename BodyFuncType>
struct normalised_function{
    BodyFuncType func;
    double Norm;
    normalised_function(BodyFuncType && F,double a,double b,size_t N_int = 1000000):func(std::forward(F)){
        Norm = integrateAB(F,a,b,N_int);
    }
    auto operator()(double x)const{
        return func(x)/Norm;
    }
};


template <typename DensFunc>
Function::VectorGrid<double> DensityGrid(double a,double b,const DensFunc & rho,size_t N,size_t iter_num = 10){
    std::vector<double> grid(N);
    grid[0] = 0;

    std::vector<double> X0tiks = Vector(N,[N](size_t i){return i/(N-1.0);});
    std::vector<double> X1tiks(N);
    X1tiks[0] = X0tiks[0] = 0;

    for(size_t iter = 0;iter<iter_num;++iter){
//        PVAR(X0tiks);
        for(size_t i=1;i<N;++i){
            double a_x = X0tiks[i-1];
            double b_x = X0tiks[i];
            X1tiks[i] = X1tiks[i-1] + integrateAB([&rho,a,b,a_x,b_x](double t){
                return 1.0/rho(a + (b-a)*(a_x + t*(b_x-a_x)));
            },0,1,20);
        }
        X1tiks *=(1.0/X1tiks.back());
        std::swap(X0tiks,X1tiks);
    }
    return Function::VectorGrid<double>(a + (b-a)*X0tiks);
}
#endif//DENSE_GRID_H
