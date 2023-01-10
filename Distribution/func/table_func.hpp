#ifndef TABLE_FUNC_H
#define TABLE_FUNC_H
#include <utils>
#include <random>
#include <cmath>
#include <tuple>
#include <type_traits>



template <typename T,typename...Args>
auto __table(const T &first,const Args&...Other){
    std::vector<T> tab(1+ arg_count<Args...>::N);
    __fill__iterator(tab.begin(),first,Other...);
    return tab;
}

template <typename T>
void printTab(const std::vector<std::vector<T>> &Tab){
    std::cout << std::scientific;
    for(size_t i=0;i<Tab[0].size();++i){
        for(size_t j=0;j<Tab.size();++j){
            std::cout << Tab[j][i] << "\t";
        }
        std::cout << std::endl;
    }
}
template <typename...Args>
void printTabVectors(const Args&...args){
    printTab(__table(args...));
}


#define JSON_VAR(x) (std::string("\"") + #x + "\"" + " : " + debugdefs::to_debug_string(x))
#define jend ""//",\n"
#define jnext ", "
template <typename FuncType2>
std::string SaveTable(double xmin,double xmax,size_t Nx,double ymin,double ymax,size_t Ny,const FuncType2 & F){
    std::ostringstream res;
    res << "{" << jend;
    {
        res << JSON_VAR(xmin) << jnext;
        res << JSON_VAR(xmax) << jnext;
        res << JSON_VAR(ymin) << jnext;
        res << JSON_VAR(ymax) << jnext;

        res << "\"F\": " << "[ ";
        {
            for(size_t i=0;i<Nx;++i){
                double bx = i/(double(Nx-1));
                double xtmp = xmin*(1-bx) + xmax*bx;
                res << "[ ";
                for(size_t j=0;j<Ny;++j){
                    double by = j/(double(Ny-1));
                    res << F(xtmp,ymin*(1-by) + ymax*by);
                    if(j < Ny-1)
                        res << ", ";
                }
                if(i < Nx-1)
                    res << "],";
                else
                    res << "]";
            }
        }
        res << "]";
    }
    res << "}" << jend;
    return res.str();
}



#endif
