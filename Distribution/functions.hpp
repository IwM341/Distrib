#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <utils>
#include <random>
#include <cmath>
#include <tuple>
#include <type_traits>

#include "func/random_defs.hpp"
#include "func/traj.hpp"



#define U0 0.7667e-3
//#define alpha 0.0073



template <class Generator>
/*MK generator of input velocity*/
inline MC::MCResult<vec3> Velocity(Generator && G,double VescTmp,
                        double Vdisp = 1.1*U0,double mU0 = U0){
    auto ksi = sqrt(-2*log(1-G()));
    auto phi = RandomPhi(G,0,M_PI);


    auto sinPhi = sin(phi);
    auto cosPhi = cos(phi);

    auto u0 = mU0/Vdisp;
    auto ve = VescTmp/Vdisp;

    auto u = sqrt(u0*u0+ksi*ksi+2*u0*ksi*cosPhi);


    auto sinTheta = (u != 0.0 ? ksi*sinPhi/u : 0);
    auto v = sqrt(u*u+ve*ve);
    auto n = RandomN(G);
    /*
    if(std::isinf((n*(v*Vdisp)).norm())){
        PVAR(ksi);
        PVAR(u);
        PVAR(v);
        PVAR(n);
    }*/
    return MC::MCResult<vec3>(n*(v*Vdisp),sinTheta*v*sqrt(M_PI/2));
}


template <class Generator>
/*MK generator of input velocity*/
inline MC::MCResult<vec3> VelocityConstrained(Generator && G,double VescTmp,double Vmin,double Vmax,
                        double Vdisp = 1.1*U0,double mU0 = U0){

    double umin = (Vmin > VescTmp ? sqrt(Vmin*Vmin - VescTmp*VescTmp)/Vdisp :  0.0);
    double umax = ( (Vmax > Vmin && Vmax > VescTmp) ? sqrt(Vmax*Vmax - VescTmp*VescTmp)/Vdisp : umin);

    double u0 = mU0/Vdisp;

    double ksi_max = umax + u0;
    double ksi_min = (mU0 - umax > 0 ? u0 - umax : ( umin-u0 > 0 ?  umin-u0 : 0 ) );

    double ksi_rd = exp(-ksi_min*ksi_min/2)-exp(-ksi_max*ksi_max/2);

    double ksi = sqrt(-2*log(exp(-ksi_max*ksi_max/2) + G()*ksi_rd));

    double cosThMin = std::min(1.0,std::max(-1.0,(umin*umin-ksi*ksi-u0*u0)/(2*ksi*u0) ));
    double cosThMax = std::max(-1.0,std::min(1.0,(umax*umax-ksi*ksi-u0*u0)/(2*ksi*u0)));

    double max_theta = acos(cosThMin);
    double min_theta = acos(cosThMax);

    double rd_th = (max_theta-min_theta)/M_PI;

    double theta = min_theta + (max_theta-min_theta)*G();

    double ve = VescTmp/Vdisp;

    double u = sqrt(u0*u0+ksi*ksi+2*u0*ksi*cos(theta));

    double sinTheta = (u != 0.0 ? ksi*sin(theta)/u : 0);

    auto v = sqrt(u*u+ve*ve);
    auto n = RandomN(G);

    return MC::MCResult<vec3>(n*(v*Vdisp),ksi_rd*rd_th*sinTheta*v*sqrt(M_PI/2));
}

template <class Generator> 
/*MK generator of output nu'*/
inline MC::MCResult<vec3> NuOut(Generator && G,const vec3& Vcm,const vec3&Nu,
                        double Vesc,double VescMin,double mp,double mk,double deltaE = 0){
	double VcmN = Vcm.norm();
	
    vec3  n_v = Vcm/VcmN;

    double cosThetaVcm = n_v.z;
    double sinThetaVcm = sqrt(n_v.x*n_v.x+n_v.y*n_v.y);

    double cosPhiVcm = 1.0;
    double sinPhiVcm = 0.0;

    if(sinThetaVcm > 1e-10){
        cosPhiVcm = n_v.x/sinThetaVcm;
        sinPhiVcm =  n_v.y/sinThetaVcm;
    }

    vec3 n_1(cosThetaVcm*cosPhiVcm,cosThetaVcm*sinPhiVcm,-sinThetaVcm);
    vec3 n_2(-sinPhiVcm,cosPhiVcm,0);
	
    /*
    std::cout << n_1*n_1 << "\t" << n_1*n_2 << "\t" << n_1*n_v << std::endl;
    std::cout << n_2*n_1 << "\t" << n_2*n_2 << "\t" << n_2*n_v << std::endl;
    std::cout << n_v*n_1 << "\t" << n_v*n_2 << "\t" << n_v*n_v << std::endl<< std::endl;
    */

    double Nu1_squared =Nu.quad()-deltaE*2*mp/(mk*(mp+mk));
	if(Nu1_squared<=0.0)
		return MC::MCResult<vec3>(vec3(0,0,0),0);
	
	double Nu1 = sqrt(Nu1_squared);
	
    double cosTh1max = (Vesc*Vesc-Nu1_squared-VcmN*VcmN)/(2*VcmN*Nu1);

    if(!(cosTh1max > -1))
		return MC::MCResult<vec3>(vec3(0,0,0),0);
	else if(cosTh1max >= 1){
		cosTh1max = 1;
	}
	
	double cosTh1 = (1+cosTh1max)*G()-1;
    double sinTh1 = sqrt(1.0-cosTh1*cosTh1);
    double phi1 = RandomPhi(G);

    const vec3 vNu1 = Nu1*(n_v*cosTh1+n_1*sinTh1*cos(phi1)+n_2*sinTh1*sin(phi1));

    /*
    if( (vNu1 + Vcm).norm() >Vesc*(1+1e-10)){
        PVAR(n_v);
        PVAR(n_1);
        PVAR(n_2);

        PVAR(cosTh1);
        PVAR(cosTh1max);
        PVAR(vNu1);
        PVAR(Vesc);

    }*/

    return MC::MCResult<vec3>(vNu1,0.5*(1.0+cosTh1max)*Nu1/VescMin);
}


/*
template <typename ScatterFuncType,typename VescRFuncType,typename N_FuncType,typename TempRFuncType,typename Generator>
inline MC::MCResult<std::tuple<vec3,double,double>> Vout(double mk,double delta_mk,ScatteringType ST,Target T,
                              ScatterFuncType const & PhiFactor, VescRFuncType const & VescR,
                             N_FuncType const & nR ,TempRFuncType const & TempR, Generator && G,
                               double Vdisp, double mU0){

    double mp = 0.938;
    double me = 0.52e-3;
    double factor = 1.0;

    //generate radius
    double r_nd = pow(G(),1.0/3.0);

    //gain escape velocity from redius
    double Vesc = VescR(r_nd);

    //random input velocity
    auto VelocityMk = Velocity(G,Vesc,Vdisp,mU0);
    auto V_wimp = VelocityMk.Result;
    factor *= VelocityMk.RemainDensity;


    double n_nd = nR(r_nd);//TODO n_nd as a function of radius
    vec3 V1 = Gauss3(G,sqrt(TempR(r_nd)*mp/2));;//TODO: add thermal distribution of nuclei velocity

    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;

    // this factor considers inelastic scaterring
    double Inelastic_rd = 1.0;

    //deltaE - enegry to ionization
    double deltaE = 0.0;
    int nMigdal;

    if(ST == IONIZATION){
        //generating ionization deltaE
        if(E_cm+delta_mk > Rd){
            //dE = max energy loss by ionization - min energyloss od ionization
            double dE = E_cm+delta_mk-Rd;
            deltaE = Rd + G()*dE;
            Inelastic_rd *= dE/Rd;
        }
        else{
            Inelastic_rd = 0;
        }
    }
    else if(ST == MIGDAL){
        double dnMx = nMax(E_cm+delta_mk)-2;
        if(dnMx > 0){
            nMigdal = 2 + G()*(int)(dnMx+1);
            deltaE = deltaEMgd(nMigdal);
        }
        else{
            deltaE = 0;
            nMigdal = -1;
            Inelastic_rd = 0;
        }
        Inelastic_rd *= (dnMx+1);
    }

    factor *= Inelastic_rd;

    // Generating out velocity
    auto Numk = NuOut(G,Vcm,Nu,Vesc,mp,mk,deltaE-delta_mk);
    vec3 Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();

    // s - appeared in exponent in scalar product of states
    double s = q/(_alpha*me);
    if(T == PROTON){
        s *= me/mp;
    }

    //Integrating factor from scalar product of states
    if(ST == IONIZATION && deltaE > Rd){
        factor *= IonFactor(s,phiMax(deltaE),dE_Rd);
    }
    else if(ST == MIGDAL && nMigdal >= 2){
        factor *= MigdalFactor(s,nMigdal);
    }
    else {
        factor *= ElasticFactor(s);
    }

    //factor from matrix element
    factor *= PhiFactor(mk,q);


    if( (Nu1+Vcm).norm() > Vesc*(1.+1e-10) && factor != 0){
        std::cout <<
            SVAR((Nu1+Vcm).norm()) + "\n" +
            SVAR(factor)+ "\n" +
            SVAR(Vesc)+ "\n" +
            SVAR(Nu1+Vcm)+ "\n\n";
    }


    return MC::MCResult<std::tuple<vec3,double,double>>(
                                                           std::tuple<vec3,double,double>(Nu1+Vcm,r_nd,Vesc),
                                                           factor);
}

*/
/*
template <typename FuncType,typename HType>
auto SupressFactor(HType & H, double mk,double delta_mk,ScatteringType ST,Target T,
                    FuncType PhiFactor,
                    const BodyModel& BM,const std::string &element,
                    size_t Nmk,
                    double Vdisp = U0,double mU0 = U0){

    //std::default_random_engine stdGen;
    //std::uniform_real_distribution<double> stdUniform(1.0,0.0);
    //auto G = [&stdGen,&stdUniform]()->double{return stdUniform(stdGen);};

    srand (time(NULL));
    auto G = [](){return  (1.0 +rand())/(RAND_MAX+1);};

    double sum = 0;
    double sum2 = 0;


    auto VescR = BM.UniformRadFunc("Vesc");
    auto TempR = BM.UniformRadFunc("Temp");
    auto NR = [](double){return 1.0;};//BM.UniformRadFunc(element);
    double VescMin = BM.VescMin();
    for(size_t i=0;i<Nmk;++i){
        auto mk_res = Vout(mk,delta_mk,ST,T,PhiFactor,VescR,NR,TempR,G,Vdisp,mU0);
        auto v_nd = std::get<0>(mk_res.Result)/VescMin;
        auto r_nd = std::get<1>(mk_res.Result);
        auto v_esc_nd = std::get<2>(mk_res.Result)/VescMin;

        double E_nd = (v_nd*v_nd - v_esc_nd*v_esc_nd);
        double L_nd = r_nd*sqrt(v_nd.x*v_nd.x + v_nd.y*v_nd.y);
        H.putValue(mk_res.RemainDensity/Nmk,E_nd,L_nd);
        /*
        if(!H.putValue(mk_res.RemainDensity/Nmk,E_nd,L_nd)){
            PVAR(v_nd.norm());
            PVAR(r_nd);
            PVAR(v_esc_nd);
            PVAR(E_nd);
            PVAR(L_nd);
            PVAR(maxLnd(BM,E_nd));
            PVAR(H.Grid[H.Grid.pos(E_nd)]);
            PVAR(H.values[H.Grid.pos(H.Grid[H.Grid.pos(E_nd)])].Grid[H.values[H.Grid.pos(E_nd)].Grid.size()-1]);
            print();
        }
        //
        sum += mk_res.RemainDensity/Nmk;
    }
    return sum;

}
*/

template <typename dF_Type,typename ScatterFuncType,typename VescRFuncType,typename N_FuncType,typename TempRFuncType,typename Generator>
inline MC::MCResult<std::tuple<vec3,double,double>> Vout1(double mp,double mk,double delta_mk,dF_Type dF,
                              ScatterFuncType const & PhiFactor, VescRFuncType const & VescR,double VescMin,
                             N_FuncType const & nR ,TempRFuncType const & TempR, Generator  && G,
                               double Vdisp, double mU0,double pow_r=1){

    double factor = 1.0;

    //generate radius
    double r_nd = pow(G(),pow_r);//pow(G(),1.0/3.0);
    factor *= (3*pow_r* pow(r_nd,(3*pow_r-1.0)/pow_r));
    //gain escape velocity from redius
    //if(r_nd < 0.3){
    //    r_nd+=0.0;
    //}
    double Vesc = VescR(r_nd);

    //random input velocity
    auto VelocityMk = Velocity(G,Vesc,Vdisp,mU0);
    auto V_wimp = VelocityMk.Result;//vec3::PolarCos(sqrt(VelocityMk.Result*VelocityMk.Result+Vesc*Vesc),
                                 //RandomCos(G),RandomPhi(G));
    factor *= VelocityMk.RemainDensity;


    double n_nd = nR(r_nd);//TODO n_nd as a function of radius
    factor *=  n_nd;
    vec3 V1 = Gauss3(G,sqrt(TempR(r_nd)/mp));//TODO: add thermal distribution of nuclei velocity

    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;



    auto EnLoss = dF.EnergyLoss(E_cm + delta_mk);
    factor *= EnLoss.RemainDensity;

    // Generating out velocity
    auto Numk = NuOut(G,Vcm,Nu,Vesc,VescMin,mp,mk,EnLoss.Result-delta_mk);
    vec3 Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();

    factor *= dF.ScatterFactor(q,EnLoss.Result);

    //factor from matrix element
    factor *= PhiFactor(q);


    #ifdef VOUT_DEBUG
    if( (Nu1+Vcm).norm() > Vesc*(1.+1e-8) && factor != 0){
        std::cout << "V out of bund:\n"+
            SVAR((Nu1+Vcm).norm()) + "\n" +
            SVAR(factor)+ "\n" +
            SVAR(Vesc)+ "\n" +
            SVAR(Nu1+Vcm)+ "\n" +
            SVAR(V1.norm()) + "\n"+
            SVAR(V_wimp.norm())+ "\n\n";
    }
    /**/
    if(factor !=0){
        PVAR(mp);
        PVAR(mk);
        PVAR(Vesc);
        PVAR(V_wimp.norm());
        PVAR(Vcm.norm());
        PVAR(Nu.norm());
        PVAR((Nu1+Vcm).norm());
        PVAR((((Nu1+Vcm).quad())-Vesc*Vesc)/4.227e-6);
        std::cin.get();
    }
    #endif
    /**/

    return MC::MCResult<std::tuple<vec3,double,double>>(
                                                           std::tuple<vec3,double,double>(Nu1+Vcm,r_nd,Vesc),
                                                           factor);
}
/*
template <typename FuncType,typename HType,typename dF_Type>
auto SupressFactor1(HType & H, double mp,double mk,double delta_mk,dF_Type dF,
                    FuncType PhiFactor,
                    const BodyModel& BM,const std::string &element,
                    size_t Nmk,
                    double Vdisp = U0,double mU0 = U0){

    //std::default_random_engine stdGen;
    //std::uniform_real_distribution<double> stdUniform(1.0,0.0);
    //auto G = [&stdGen,&stdUniform]()->double{return stdUniform(stdGen);};

    srand (time(NULL));
    auto G = [](){return  (1.0 +rand())/(RAND_MAX+1);};

    double sum = 0;
    double sum2 = 0;


    auto VescR = BM.UniformRadFunc("Vesc");
    auto TempR = BM.UniformRadFunc("Temp");
    auto NR = [](double){return 1.0;};//BM.UniformRadFunc(element);
    double VescMin = BM.VescMin();
    for(size_t i=0;i<Nmk;++i){
        auto mk_res = Vout1(mp,mk,delta_mk,dF,PhiFactor,VescR,VescMin,NR,TempR,G,Vdisp,mU0);
        auto v_nd = std::get<0>(mk_res.Result)/VescMin;
        auto r_nd = std::get<1>(mk_res.Result);
        auto v_esc_nd = std::get<2>(mk_res.Result)/VescMin;

        double E_nd = (v_nd*v_nd - v_esc_nd*v_esc_nd);
        double L_nd = r_nd*sqrt(v_nd.x*v_nd.x + v_nd.y*v_nd.y);
        if(H.putValue(mk_res.RemainDensity/Nmk,E_nd,L_nd))
            sum += mk_res.RemainDensity/Nmk;
        else{
           if(mk_res.RemainDensity!=0){
               PVAR(v_nd.norm());
               PVAR(r_nd);
               PVAR(v_esc_nd);
               PVAR(E_nd);
               PVAR(L_nd);
               PVAR(maxLnd(BM,E_nd));
           }
        }
        /*
        if(!H.putValue(mk_res.RemainDensity/Nmk,E_nd,L_nd)){
            PVAR(v_nd.norm());
            PVAR(r_nd);
            PVAR(v_esc_nd);
            PVAR(E_nd);
            PVAR(L_nd);
            PVAR(maxLnd(BM,E_nd));
            PVAR(H.Grid[H.Grid.pos(E_nd)]);
            PVAR(H.values[H.Grid.pos(H.Grid[H.Grid.pos(E_nd)])].Grid[H.values[H.Grid.pos(E_nd)].Grid.size()-1]);
            print();
        }
        //

    }
    return sum;

}
*/
/*!
 * \brief SupressFactor_v1
 * \param H histo of output values
 * \param mi nuclei mass
 * \param mk DM mass
 * \param delta_mk in_mass - out_mass
 * \param dF form factor
 * \param NR concentration(r)
 * \param VescMin
 * \param VescR
 * \param TempR
 * \param PhiFactor form factor
 * \param G generator from 0 to 1 (not including)
 * \param Nmk
 * \param pow_r (double pow, increase leads to generation r in perepherial of sphere)
 * \param Vdisp
 * \param mU0
 */
template <typename VescFuncType,typename ThermFuncType,typename NFuncType,
          typename PhiFuncType,typename HType,typename dF_Type,typename Generator>
auto SupressFactor_v1(HType & H, double mi,double mk,double delta_mk,dF_Type dF,
                      NFuncType const &NR,double VescMin,VescFuncType const & VescR,ThermFuncType const & TempR,
                    PhiFuncType PhiFactor,Generator const&G,
                    size_t Nmk,double pow_r,
                    double Vdisp = U0,double mU0 = U0){

    double sum = 0;
    double sum2 = 0;
    double const_fact_rd = (_mp*fast_pow((_mp+mk),2))/(mi*mi*(mi+mk))/Nmk;//mk/(mk+mi)/Nmk;
    for(size_t i=0;i<Nmk;++i){
        auto mk_res = Vout1(mi,mk,delta_mk,dF,PhiFactor,VescR,VescMin,NR,TempR,G,Vdisp,mU0,pow_r);
        auto v_nd = std::get<0>(mk_res.Result)/VescMin;
        auto r_nd = std::get<1>(mk_res.Result);
        auto v_esc_nd = std::get<2>(mk_res.Result)/VescMin;

        double E_nd = (v_nd*v_nd - v_esc_nd*v_esc_nd);
        double L_nd = r_nd*sqrt(v_nd.x*v_nd.x + v_nd.y*v_nd.y);
        double dens = mk_res.RemainDensity*const_fact_rd;


        //double Ewas = H.values[1].values[0];
        /*
        auto l = L_nd/H.LE_func(E_nd);
        auto [b,MI] = H.Histo.Grid.spos(E_nd,l);
        auto i_h = H.Histo.Grid.LinearIndex(MI);
        if(b){
            auto Elem = H.Histo.Grid[MI];
            auto H_i = H.Histo.Values[i_h];
        }*/
        /*if(E_nd < -4.5){
            dens = dens +0.0;
            print("E<-4.5");
        }*/
        if(H.putValue(dens,E_nd,L_nd)){
            /*
            auto H_i_1 = H.Histo.Values[i_h];
            */
            sum += dens;
        }/*
        else{
            auto l = L_nd/H.LE_func(E_nd);
        }*/
         /*
        if(Ewas != H.values[1].values[0]){
            PVAR(dens);
            PVAR(v_nd);
            PVAR(r_nd);
            PVAR(v_esc_nd);
            PVAR(E_nd);
            PVAR(L_nd);
            PVAR(H.values[1].values[1]);
            print();
        }*/
            /*
        if(std::isnan(E_nd) or std::isnan(L_nd)){
            PVAR(r_nd);
            PVAR(v_nd);
            PVAR(v_esc_nd);
            PVAR(v_esc_nd);
        }*/
    }
    return sum;

}


template <class Generator>
/*MK generator of output nu'*/
inline auto NuOutTherm(Generator const & G,const vec3& Vcm,const vec3&Nu,double VescMin,
                                     double mp,double mk,double deltaE = 0){
    //double VcmN = Vcm.norm();

    vec3  n_v = Vcm.normalized();

    double cosThetaVcm = n_v.z;
    double sinThetaVcm = sqrt(n_v.x*n_v.x+n_v.y*n_v.y);

    double cosPhiVcm = 1.0;
    double sinPhiVcm = 0.0;

    if(sinThetaVcm > 1e-10){
        cosPhiVcm = n_v.x/sinThetaVcm;
        sinPhiVcm =  n_v.y/sinThetaVcm;
    }

    vec3 n_1(cosThetaVcm*cosPhiVcm,cosThetaVcm*sinPhiVcm,-sinThetaVcm);
    vec3 n_2(-sinPhiVcm,cosPhiVcm,0);

    /*
    std::cout << n_1*n_1 << "\t" << n_1*n_2 << "\t" << n_1*n_v << std::endl;
    std::cout << n_2*n_1 << "\t" << n_2*n_2 << "\t" << n_2*n_v << std::endl;
    std::cout << n_v*n_1 << "\t" << n_v*n_2 << "\t" << n_v*n_v << std::endl<< std::endl;
    */

    double Nu1_squared =Nu.quad()-deltaE*2*mp/(mk*(mp+mk));
    if(Nu1_squared<=0.0)
        return MC::MCResult<vec3>(vec3(0,0,0),0);

    double Nu1 = sqrt(Nu1_squared);



    double cosTh1 = RandomCos(G);
    double sinTh1 = sqrt(1.0-cosTh1*cosTh1);
    double phi1 = RandomPhi(G);

    const vec3 vNu1 = Nu1*(n_v*cosTh1+n_1*sinTh1*cos(phi1)+n_2*sinTh1*sin(phi1));

    /*
    if( (vNu1 + Vcm).norm() >Vesc*(1+1e-10)){
        PVAR(n_v);
        PVAR(n_1);
        PVAR(n_2);

        PVAR(cosTh1);
        PVAR(cosTh1max);
        PVAR(vNu1);
        PVAR(Vesc);

    }*/

    return MC::MCResult<vec3>(vNu1,Nu1/VescMin);
}





/*
template <typename ScatterFuncType,typename Generator>
inline MC::MCResult<vec3> VoutTherm(double mk,double delta_mk,ScatteringType ST,Target T,
                              ScatterFuncType const &PhiFactor, const vec3& V_wimp, double n_r,double Therm, Generator && G){

    double mp = 0.938;
    double me = 0.52e-3;
    double factor = 1.0;

    factor *= n_r;
    vec3 V1 = Gauss3(G,sqrt(2*Therm/mp));

    /*
    PVAR(Therm);
    if(V1.norm()>2e-3){
        PVAR(V1);
        PVAR(V_wimp);
    }
    //
    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;

    // this factor considers inelastic scaterring
    double Inelastic_rd = 1.0;

    //deltaE - enegry to ionization
    double deltaE = 0.0;
    int nMigdal;

    if(ST == IONIZATION){
        //generating ionization deltaE
        if(E_cm+delta_mk > Rd){
            //dE = max energy loss by ionization - min energyloss od ionization
            double dE = E_cm+delta_mk-Rd;
            deltaE = Rd + G()*dE;
            Inelastic_rd *= dE/Rd;
        }
        else{
            Inelastic_rd = 0;
        }
    }
    else if(ST == MIGDAL){
        double dnMx = nMax(E_cm+delta_mk)-2;
        if(dnMx > 0){
            nMigdal = 2 + G()*(int)(dnMx+1);
            deltaE = deltaEMgd(nMigdal);
        }
        else{
            deltaE = 0;
            nMigdal = -1;
            Inelastic_rd = 0;
        }
        Inelastic_rd *= (dnMx+1);
    }

    factor *= Inelastic_rd;

    // Generating out velocity
    auto Numk = NuOutTherm(G,Vcm,Nu,mp,mk,deltaE-delta_mk);
    vec3 Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();

    // s - appeared in exponent in scalar product of states
    double s = q/(_alpha*me);
    if(T == PROTON){
        s *= me/mp;
    }

    //Integrating factor from scalar product of states
    if(ST == IONIZATION && deltaE > Rd){
        factor *= IonFactor(s,phiMax(deltaE),dE_Rd);
    }
    else if(ST == MIGDAL && nMigdal >= 2){
        factor *= MigdalFactor(s,nMigdal);
    }
    else {
        factor *= ElasticFactor(s);
    }

    //factor from matrix element
    factor *= PhiFactor(mk,q);


    /*
    if( (Nu1+Vcm).norm() > Vesc*(1.+1e-10) && factor != 0){
        std::cout <<
            SVAR((Nu1+Vcm).norm()) + "\n" +
            SVAR(factor)+ "\n" +
            SVAR(Vesc)+ "\n" +
            SVAR(Nu1+Vcm)+ "\n\n";
    }
    //
    if(std::isnan(factor)){
        PVAR(V_wimp);
        PVAR(V1);
        PVAR(Therm);
        //PVAR(PhiFactor(mk,q));
        //PVAR(Numk.RemainDensity);
        //PVAR(Inelastic_rd);
        //PVAR(ElasticFactor(s));
        //PVAR(Vcm);
        //PVAR(Nu);
        //PVAR(E_cm);
        //PVAR(Numk.Result);
        //PVAR(q);
        print();
    }
    return MC::MCResult<vec3>(Nu1+Vcm,factor);
}
*/
template <typename ScatterFuncType,typename Generator,typename dF_Factor_Type>
inline MC::MCResult<vec3> VoutTherm1(double mk,double mp,double delta_mk,dF_Factor_Type dF,
                              ScatterFuncType const &PhiFactor, const vec3& V_wimp,
                                     double VescMin,double n_r,double Therm, Generator  && G){

    double factor = 1.0;

    factor *= n_r;
    vec3 V1 = Gauss3(G,sqrt(Therm/mp));//TODO: add thermal distribution of nuclei velocity

    //PVAR(Therm);
//    if(V1.norm()>2e-3){
//        PVAR(V_wimp);
//        PVAR(V1);
//        PVAR(Therm);
//    }

    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;

    auto EnLoss = dF.EnergyLoss(E_cm + delta_mk);
    factor *= EnLoss.RemainDensity;

    // Generating out velocity
    auto Numk = NuOutTherm(G,Vcm,Nu,VescMin,mp,mk,EnLoss.Result-delta_mk);
    vec3 const&Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();
    factor *= dF.ScatterFactor(q,EnLoss.Result);
    factor *=  PhiFactor(q);
    /*
    if(std::isnan(factor)){
        int __var__ = 10;
    }*/
    return MC::MCResult<vec3>(Nu1+Vcm,factor);
}


auto Perps(const vec3 & V){
    double ax = std::abs(V.x);
    double ay = std::abs(V.y);
    double az = std::abs(V.z);

    vec3 v1,v2;

    if(ax >= ay && ax >= az){
        double n1 = sqrt(ay*ay+az*az);
        if(n1 == 0){
            return _T(vec3(0,1,0),vec3(0,0,1));
        }
        v1.y = V.z/n1;
        v1.z = -V.y/n1;

        v2 = vec3(-(ay*ay+az*az)/V.x,V.y,V.z);
        v2 /= v2.norm();

        return _T(v1,v2);
    }
    else if(ay >= ax && ay >= az){
        double n1 = sqrt(ax*ax+az*az);
        if(n1 == 0){
            return _T(vec3(1,0,0),vec3(0,0,1));
        }
        v1.x = V.z/n1;
        v1.z = -V.x/n1;

        v2 = vec3(V.x,-(ax*ax+az*az)/V.y,V.z);
        v2 /= v2.norm();
        return _T(v1,v2);
    }
    else{
        double n1 = sqrt(ax*ax+ay*ay);
        if(n1 == 0){
            return _T(vec3(1,0,0),vec3(0,1,0));
        }
        v1.x = V.y/n1;
        v1.y = -V.x/n1;

        v2 = vec3(V.x,V.y,-(ax*ax+ay*ay)/V.z);
        v2 /= v2.norm();
        return _T(v1,v2);
    }
}

//generation of therm velocity to cross the treshold
template <typename ScatterFuncType,typename Generator,typename dF_Factor_Type>
inline MC::MCResult<vec3> VoutTherm_Opt(double mk,double mp,double delta_mk,dF_Factor_Type dF,
                              ScatterFuncType const &PhiFactor, const vec3& V_wimp,
                                        double VescMin,double n_r,double Therm, Generator  && G){

    double factor = 1.0;

    factor *= n_r;
    vec3 V1;
    if(delta_mk >= 0)
        V1 = Gauss3(G,sqrt(2*Therm/mp));
    else{
        double dVmin = sqrt(-2*delta_mk*(mp+mk)/(mk*mp));
        double V_w_n = V_wimp.norm();

        double V1min = (dVmin > V_w_n ? dVmin-V_w_n : 0);

        auto V1_mk = Gauss3_min_abs(G,sqrt(2*Therm/mp),V1min);
        double V1_abs = V1_mk.Result;
        factor *= V1_mk.RemainDensity;

        double cosTheta_max = (V1_abs*V1_abs+V_w_n*V_w_n-dVmin*dVmin)/(2*V_w_n*V1_abs);
        if(cosTheta_max >= 1){
            cosTheta_max = 1;
        }
        else if(cosTheta_max <= -1){
            cosTheta_max = -1.0;
        }

        auto cosTheta_mk = RandomCos_c(G,-1.0,cosTheta_max);
        factor *= cosTheta_mk.RemainDensity;

        vec3 V1_proto = vec3::PolarCos(1,cosTheta_mk.Result,RandomPhi(G));
        auto v1v2 = Perps(V_wimp);
        V1 = V1_abs*(V_wimp*(V1_proto.z/V_w_n) + std::get<0>(v1v2)*V1_proto.x+std::get<1>(v1v2)*V1_proto.y);

        if(std::isnan(V1.norm()))
        {
            PVAR(V_wimp);
            PVAR(Therm);

            //PVAR(Vesc);

            PVAR(dVmin);
            PVAR(V_w_n);
            PVAR(V1_abs);
            PVAR(cosTheta_max);
            PVAR(V1_proto);
            PVAR(std::get<0>(v1v2));
            PVAR(std::get<1>(v1v2));
        }

    }



    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;


    auto EnLoss = dF.EnergyLoss(E_cm + delta_mk);
    factor *= EnLoss.RemainDensity;

    // Generating out velocity
    auto Numk = NuOutTherm(G,Vcm,Nu,VescMin,mp,mk,EnLoss.Result-delta_mk);
    vec3 Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();
    factor *= dF.ScatterFactor(q,EnLoss.Result);
    factor *= PhiFactor(q);
    return MC::MCResult<vec3>(Nu1+Vcm,factor);
}


/*
template <typename HistoType,typename Generator,typename ThermFuncType,typename NFuncType,typename VescFuncType,typename ScatterFuncType>
inline void TrajectoryIntegral(Generator const & G,
                        double e_nd,double l_nd,const TrajectoryInfo & TI,double VescMin,NFuncType const& nR,const ThermFuncType & ThermR,
                               const VescFuncType & Vesc,
                        HistoType &Out,double & EvaporationOut,double mk,double delta_mk,
                               ScatteringType ST,Target T,ScatterFuncType const & PhiFactor,size_t Nmk)
{
    if(TI.T_in == 0.0){
        return;
    }

    for(size_t sch = 0;sch<Nmk;++sch){
        double factor = TI.T_in/(TI.T_in+TI.T_out);
        double t = G()*TI.T_in;
        double r = TI.Trajectory(t);
        double vesc = Vesc(r);
        double Therm = ThermR(r);
        double n_r = nR(r);

//        if(std::isnan(r+vesc+Therm+n_r)){
//            PVAR(r);
//            PVAR(vesc);
//            PVAR(Therm);
//            PVAR(n_r);
//            PVAR(TI.Trajectory.toString());
//            exit(0);
//        }

        double v2 = e_nd*VescMin*VescMin+vesc*vesc;
        if(v2 < 0)
            v2 = 0;
        double v = sqrt(v2);
        double v_tau  = l_nd/r;
        double v_r = (v_tau < v) ? sqrt(v*v-v_tau*v_tau) : 0;
        vec3 Vin = vec3(v_tau,0,v_r);//RandomN(G)*v;
        auto Vout = VoutTherm(mk,delta_mk,ST,T,PhiFactor,Vin,n_r,Therm,G);

        double e_nd_1 = (Vout.Result*Vout.Result - vesc*vesc)/(VescMin*VescMin);
        double l_nd_1  = r*sqrt(Vout.Result.x*Vout.Result.x+Vout.Result.y*Vout.Result.y)/VescMin;


        double dens = Vout.RemainDensity/Nmk;
        if(!Out.putValue(dens,e_nd_1,l_nd_1)){
            EvaporationOut += dens;
        }
        else{
            if(Vout.Result.norm()>vesc){
                PVAR(Therm);
                PVAR(vesc);
                PVAR(Vout.Result);
            }
        }
    }
}
*/

/*!
 * \brief TrajectoryIntegral1
 * \param G generator of values from 0 to 1 (1 is not included)
 * \param e_nd
 * \param l_nd
 * \param TI
 * \param VescMin
 * \param nR - concentration of nuclie(r)
 * \param ThermR T(t)
 * \param Vesc
 * \param Out histo to put values of scaterred particles
 * \param EvaporationOut double & value to add evaporaton count
 * \param mk mass of WIMP
 * \param mp mass of nuclei
 * \param delta_mk = in_mass - out_mass
 * \param dF  form factor
 * \param PhiFactor form factor
 * \param Nmk MK samples
 */
template <typename HistoType,typename Generator,typename ThermFuncType,typename NFuncType,typename VescFuncType,typename dF_Type,typename ScatterFuncType>
inline double TrajectoryIntegral1(Generator const & G,
                        double e_nd,double l_nd,const TrajectoryInfo & TI,double VescMin,NFuncType const& nR,const ThermFuncType & ThermR,
                               const VescFuncType & Vesc,
                        HistoType &Out,double & EvaporationOut,double mk,double mi,double delta_mk,
                               dF_Type dF,ScatterFuncType const & PhiFactor,size_t Nmk)
{
    if(TI.T_in == 0.0){
        return 0.0;
    }
    double const_fact_rd = (_mp*fast_pow((_mp+mk),2))/(mi*mi*(mi+mk))/Nmk*TI.T_in/(TI.T_in+TI.T_out);
    double summ =0;
    for(size_t sch = 0;sch<Nmk;++sch){
        double factor = const_fact_rd;
        double t = G()*TI.T_in;
        double r = TI.Trajectory(t);
        double vesc = Vesc(r);
        double Therm = ThermR(r);
        double n_r = nR(r);

        double v2 = e_nd*VescMin*VescMin+vesc*vesc;
        double v_xy = (l_nd/r)*VescMin;
        if(!(v_xy >= 0))
            v_xy = 0;
        double v2_z = v2 - v_xy*v_xy;
        if(!(v2_z >= 0))
            v2_z = 0;
        double Phi = RandomPhi(G);
        vec3 Vin(v_xy*cos(Phi),v_xy*sin(Phi),sqrt(v2_z));
//        if(std::isnan(Vin.norm())){
//            throw std::runtime_error("nan at Vin");
//        }
        auto Vout = VoutTherm1(mk,mi,delta_mk,dF,PhiFactor,Vin,VescMin,n_r,Therm,G);

        double e_nd_1 = (Vout.Result*Vout.Result - vesc*vesc)/(VescMin*VescMin);
        double l_nd_1  = r*sqrt(Vout.Result.x*Vout.Result.x+Vout.Result.y*Vout.Result.y)/VescMin;


        double dens = Vout.RemainDensity*factor;//*const_fact_rd;
        if(!Out.putValue(dens,e_nd_1,l_nd_1)){
            //PVAR(dens);
            EvaporationOut += dens;
        }else{
            summ += dens;
        }
    }
    return summ;
}

template <typename HistoType,typename Generator,typename ThermFuncType,typename NFuncType,typename VescFuncType,typename dF_Type,typename ScatterFuncType>
inline void TrajectoryIntegral1_Opt(Generator const & G,
                        double e_nd,double l_nd,const TrajectoryInfo & TI,double VescMin,NFuncType const& nR,const ThermFuncType & ThermR,
                               const VescFuncType & Vesc,
                        HistoType &Out,double & EvaporationOut,double mk,double mp,double delta_mk,
                               dF_Type dF,ScatterFuncType const & PhiFactor,size_t Nmk)
{
    if(TI.T_in == 0.0){
        return;
    }

    for(size_t sch = 0;sch<Nmk;++sch){
        double factor = TI.T_in/(TI.T_in+TI.T_out);
        double t = G()*TI.T_in;
        double r = TI.Trajectory(t);
        double vesc = Vesc(r);
        double Therm = ThermR(r);
        double n_r = nR(r);

        double v2 = e_nd*VescMin*VescMin+vesc*vesc;
        if(v2 < 0)
            v2 = 0;
        double v = sqrt(v2);
        vec3 Vin = RandomN(G)*v;
        auto Vout = VoutTherm_Opt(mk,mp,delta_mk,dF,PhiFactor,Vin,n_r,Therm,G);

        double e_nd_1 = (Vout.Result*Vout.Result - vesc*vesc)/(VescMin*VescMin);
        double l_nd_1  = r*sqrt(Vout.Result.x*Vout.Result.x+Vout.Result.y*Vout.Result.y)/VescMin;


        double dens = Vout.RemainDensity/Nmk;
        if(!Out.putValue(dens,e_nd_1,l_nd_1)){
            //PVAR(dens);
            EvaporationOut += dens;
        }
    }
}


double Ion_Deg(double f_p,double T_fact){
    return 2/(1+sqrt(1+4*f_p/T_fact));
}

double Atom_ion_Deg(double f_p,double T_fact){
    double sqrt_denom = Ion_Deg(f_p,T_fact);
    return sqrt_denom*sqrt_denom*(f_p/T_fact);
}

template <typename ScatterFuncType,typename Generator>
inline MC::MCResult<vec3> VoutTherm_AntiIonization(double mk,double mp,double delta_mk,Target T,
                              ScatterFuncType const &PhiFactor, const vec3& V_wimp, double n_r,double f_e,double Therm, Generator  && G){

    double factor = 1.0;

    factor *= n_r;
    vec3 V1 = Gauss3(G,sqrt(2*Therm/mp));

    //generation of external electron
    factor *= 2/(_alpha*_alpha);
    factor *= f_e;
    double Ve = Gauss2Norm(G,sqrt(2*Therm/_me));

    double dE_ion = Rd + 0.5*_me*Ve*Ve;
    double phi = Ve/_alpha;
    //

    //Vcm - is a vrlocity of momentum center
    vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
    //Nu is input velocity of WIMP in cm coordinatesd
    vec3 Nu = mp/(mp+mk)*(V_wimp-V1);

    //Ecm - is kinetic energy in cm
    double E_cm = mk*(mp+mk)/mp*Nu.quad()/2;

    // Generating out velocity
    auto Numk = NuOutTherm(G,Vcm,Nu,mp,mk,-dE_ion-delta_mk);
    vec3 Nu1 = Numk.Result;
    factor*=Numk.RemainDensity;

    // q - exchange momentum
    double q = mk*(Nu-Nu1).norm();
    double s = q/(_alpha*(T == PROTON ? _mp : _me));
    factor *= IonFactor(s,phi,dE_Rd);

    return MC::MCResult<vec3>(Nu1+Vcm,factor);
}
//f_e_R = (concentration)/(m_e^3) - dimentionnless parametr
template <typename HistoType,typename Generator,typename ThermFuncType,typename NFuncType,typename f_e_FuncType,typename VescFuncType,typename ScatterFuncType>
inline void TrajectoryIntegral_AntiIonization(Generator const & G,
                        double e_nd,double l_nd,const TrajectoryInfo & TI,NFuncType const& n_H_R,f_e_FuncType const& f_e_R,
                                    const ThermFuncType & ThermR,
                                    const VescFuncType & Vesc,double VescMin,
                                    HistoType &Out,double & EvaporationOut,Target T,double mk,double mp,double delta_mk,
                                    ScatterFuncType const & PhiFactor,size_t Nmk)
{
    if(TI.T_in == 0.0){
        return;
    }

    for(size_t sch = 0;sch<Nmk;++sch){
        double factor = TI.T_in/(TI.T_in+TI.T_out);
        double t = G()*TI.T_in;
        double r = TI.Trajectory(t);
        double vesc = Vesc(r);
        double Therm = ThermR(r);
        double n_r = n_H_R(r);
        double f_e_r = f_e_R(r);

        double v2 = e_nd*VescMin*VescMin+vesc*vesc;
        if(v2 < 0)
            v2 = 0;
        double v = sqrt(v2);
        vec3 Vin = RandomN(G)*v;
        auto Vout = VoutTherm_AntiIonization(mk,mp,delta_mk,T,PhiFactor,Vin,n_r,f_e_r,Therm,G);

        double e_nd_1 = (Vout.Result*Vout.Result - vesc*vesc)/(VescMin*VescMin);
        double l_nd_1  = r*sqrt(Vout.Result.x*Vout.Result.x+Vout.Result.y*Vout.Result.y)/VescMin;


        double dens = Vout.RemainDensity/Nmk;
        if(!Out.putValue(dens,e_nd_1,l_nd_1)){
            //PVAR(dens);
            EvaporationOut += dens;
        }
    }
}
#endif
