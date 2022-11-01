#ifndef DIFF_SOLVE_HPP
#define DIFF_SOLVE_HPP

#include <utils>

template <typename T,typename FuncType>
T step(T x,double h,const FuncType & F){
    T x_mid = x + 0.5*h*F(x);
    return x + h*F(x_mid);
}

template <typename T,typename FuncType>
std::vector<T> solve_diff(T x0,const FuncType & F,size_t N_steps,double h){
    std::vector<T> ret(N_steps);
    ret[0] = x0;
    for(size_t i=1;i< N_steps;++i){
        ret[i] = step(ret[i-1],h,F);
    }
    return ret;
}

#endif // DIFF_SOLVE_HPP
