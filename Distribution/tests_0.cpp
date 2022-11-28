#include <utils>
#include "functions.hpp"

bool eq(double x,double y){
    return std::abs(x-y) < 1e-15;
}
bool is_basis(const vec3 &v0,const vec3 &v1,const vec3 &v2){
    return (eq(v0*v0,1) && eq(v1*v1,1) && eq(v2*v2,1) && eq(v1*v0,0) && eq(v2*v0,0) && eq(v1*v2,0));
}

int main(void){
    auto G = [](){return rand()/(RAND_MAX+1.0);};
    vec3 v1,v2;
    for(size_t i=0;i<100;++i){
        vec3 V = vec3::PolarCos(G(),RandomCos(G),RandomPhi(G));

        _R(v1,v2) = Perps(V);

        if(!is_basis(V/V.norm(),v1,v2)){
            PVAR(V);
            PVAR(v1);
            PVAR(v1.norm());
            PVAR(v2);
            PVAR(v2.norm());
            PVAR(v1*V);
            PVAR(v2*V);
            PVAR(v1*v2);
            print();
        }


    }
};
