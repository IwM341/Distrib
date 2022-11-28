#ifndef RANDOM_DEFS_H
#define RANDOM_DEFS_H
#include "../functions.hpp"
template <class Generator>
inline double RandomCos(Generator && G){
    return G()*2.0-1;
}

template <class Generator>
inline MC::MCResult<double> RandomCos_c(Generator && G,double cosTheta_min,double cosTheta_max){
    return MC::MCResult<double>(cosTheta_min + G()*(cosTheta_max-cosTheta_min),0.5*(cosTheta_max-cosTheta_min));
}

template <class Generator>
inline double RandomPhi(Generator && G){
    return G()*2*M_PI;
}

template <class Generator>
inline double RandomPhi(Generator && G,double phi0,double phi1){
    return phi0 + G()*(phi1-phi0);
}

template <class Generator>
inline vec3 RandomN(Generator const & G){
    return vec3::PolarCos(1.0,RandomCos(G),RandomPhi(G));
}

template <class Generator>
inline MC::MCResult<vec3> Boltsman(Generator&& G, double Vdisp,double Vmin = 0){
    double V = sqrt(Vmin*Vmin-2*Vdisp*Vdisp*log(1-G()));
    return MC::MCResult<vec3>(vec3::PolarCos(V,RandomCos(G),RandomPhi(G)),
                exp(-Vmin*Vmin/(2*Vdisp*Vdisp))*sqrt(2/M_PI)*V/Vdisp);
}

template <class Generator>
inline double Gauss2Norm(Generator && G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(1-G()));
    return V;
}

template <class Generator>
inline vec3 Gauss2(Generator && G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(1-G()));
    return vec3::PolarCos(V,0,RandomPhi(G));
}

template <class Generator>
inline double Gauss(Generator && G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(1-G()));
    return V*cos(RandomPhi(G));
}
template <class Generator>
inline vec3 Gauss3(Generator && G, double Vdisp){
    double V1 = Vdisp*sqrt(-2*log(1-G()));
    double V2 = Vdisp*sqrt(-2*log(1-G()));
    double phi1 = RandomPhi(G);
    double phi2 = RandomPhi(G);
    return vec3(V1*cos(phi1),V1*sin(phi1),V2*cos(phi2));
}

template <class Generator>
inline MC::MCResult<vec3> Gauss3_min(Generator && G, double Vdisp, double Vmin){
    Vmin = (Vmin >0 ? Vmin : 0.0 );
    double a0 = exp(-Vmin*Vmin/(2*Vdisp*Vdisp));
    double v_nd = sqrt(-2*log(a0*(1.0-G())));
    return MC::MCResult<vec3>(vec3::PolarCos(v_nd*Vdisp,RandomCos(G),RandomPhi(G)),sqrt(2.0/M_PI)*v_nd*a0);
}

template <class Generator>
inline MC::MCResult<double> Gauss3_min_abs(Generator && G, double Vdisp, double Vmin){
    Vmin = (Vmin >0 ? Vmin : 0.0 );
    double a0 = exp(-Vmin*Vmin/(2*Vdisp*Vdisp));
    double v_nd = sqrt(-2*log(a0*(1.0-G())));
    return MC::MCResult<double>(v_nd*Vdisp,sqrt(2.0/M_PI)*v_nd*a0);
}

#endif
