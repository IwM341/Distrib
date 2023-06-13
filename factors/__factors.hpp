#ifndef __FACTORS_HPP
#define __FACTORS_HPP
#include <cmath>
constexpr double fermi_GeV = 5;
constexpr double __mp__ = 0.938;
constexpr double __me__ = 0.51e-3;

inline double BesselX(double x) noexcept{
    if(x<0.01){
        return 1.0/3-x*x*(1-x*x/28)/10;
    }else{
        return (sin(x)-x*cos(x))/(x*x*x);
    }
}
struct dF_Nuc_M_SI{
    size_t M;
    double const_factor;
    double s,R;
    dF_Nuc_M_SI(size_t A,size_t Z,size_t const_factor = 1):M(A),const_factor(const_factor){
        s =  fermi_GeV*0.9;
        double b = (1.23*pow(M,1.0/3)-0.6)*fermi_GeV;
        double a = 0.52*fermi_GeV;
        R = sqrt(b*b+7*M_PI*M_PI*a*a/3-5*s*s);
    }
    inline MC::MCResult<double> EnergyLoss(double Emax)const noexcept{
        return MC::MCResult<double>(0,1);
    }
    inline double ScatterFactor(double q,double enLoss)const noexcept{
        double bf = 3*BesselX(q*R)*exp(-q*q*s*s/2);
        return const_factor*bf*bf;
    }
};


struct Phi_Factor_Scatter_SI{
   double fac;
   Phi_Factor_Scatter_SI(double mi,double mk){
       fac = (mi/__mp__);
       fac *= fac;
       fac *= fac;
   }
   inline double operator ()(double q) const noexcept{
       return fac;
   }
};
struct Phi_Factor_Scatter_SD{
    double q_fac;
    constexpr static double v_0 = 2e-3;
    Phi_Factor_Scatter_SD(double mi,double mk){
        q_fac = (mi/__mp__);
        q_fac *= q_fac;
        q_fac *= q_fac;
        q_fac /= (v_0*v_0*v_0*v_0*mk*mk*mk*mk);
        if(mi > 1.5){
            q_fac = 0;
        }
    }
    inline double operator ()(double q) const noexcept{
        double f = (q*q);
        return q_fac*f*f;
    }
};

struct Phi_Factor_Scatter_Electron_SD {

    constexpr static double v_0 = 2e-3;
    double e_fac;
    Phi_Factor_Scatter_Electron_SD(double mi,double mk){
        e_fac = mi/__me__;
        e_fac *= e_fac;
        double q_fac = (mi/__mp__);
        q_fac *= q_fac;
        q_fac *= q_fac;
        q_fac /= (v_0*v_0*v_0*v_0*mk*mk*mk*mk);
        e_fac *= q_fac;
        if(mi > 1.5){
            e_fac = 0;
        }
    }
    inline double operator ()(double q) const noexcept{
        double f = (q*q);
        return e_fac*f*f;
    }
};
struct Phi_Factor_Annihilation_1{
    inline double operator ()(double v_diff) const noexcept{
        return 1;
    }
};

struct Phi_Factor_Annihilation_2{
    constexpr static double v_0 = 2e-3;
    inline double operator ()(double v_diff) const noexcept{
        double s = (v_diff/v_0);
        s *= s;
        s *= s;
        return s;
    }
};

#endif//FACTORS_HPP
