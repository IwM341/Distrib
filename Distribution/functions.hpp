#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <utils>
#include <random>
#include <cmath>
#include <tuple>
#include <type_traits>

#define U0 0.7667e-3
#define alpha 0.0073

template <class Generator> 
inline double RandomCos(Generator const & G){
	return G()*2.0-1;
}
template <class Generator> 
inline double RandomPhi(Generator const & G){
	return G()*2*M_PI;
}

template <class Generator>
inline double RandomPhi(Generator const & G,double phi0,double phi1){
    return phi0 + G()*(phi1-phi0);
}

template <class Generator>
inline vec3 RandomN(Generator const & G){
    return vec3::PolarCos(1.0,RandomCos(G),RandomPhi(G));
}

template <class Generator> 
inline MC::MCResult<vec3> Boltsman(Generator G, double Vdisp,double Vmin = 0){
	double V = sqrt(Vmin*Vmin-2*Vdisp*Vdisp*log(G()));
    return MC::MCResult<vec3>(vec3::PolarCos(V,RandomCos(G),RandomPhi(G)),
				exp(-Vmin*Vmin/(2*Vdisp*Vdisp))*sqrt(2/M_PI)*V/Vdisp);
}

template <class Generator>
inline double Gauss2Norm(Generator const & G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(G()));
    return V;
}

template <class Generator>
inline vec3 Gauss2(Generator const & G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(G()));
    return vec3::PolarCos(V,0,RandomPhi(G));
}

template <class Generator>
inline double Gauss(Generator const & G, double Vdisp){
    double V = Vdisp*sqrt(-2*log(G()));
    return V*cos(RandomPhi(G));
}
template <class Generator>
inline vec3 Gauss3(Generator G, double Vdisp){
    double V1 = Vdisp*sqrt(-2*log(G()));
    double V2 = Vdisp*sqrt(-2*log(G()));
    double phi1 = RandomPhi(G);
    double phi2 = RandomPhi(G);
    return vec3(V1*cos(phi1),V1*sin(phi1),V2*cos(phi2));
}

template <class Generator>
/*MK generator of input velocity*/
inline MC::MCResult<vec3> Velocity(Generator const & G,double VescTmp,
                        double Vdisp = 1.1*U0,double mU0 = U0){
    auto ksi = sqrt(-2*log(G()));
    auto phi = RandomPhi(G,0,M_PI);


    auto sinPhi = sin(phi);
    auto cosPhi = cos(phi);

    auto u0 = mU0/Vdisp;
    auto ve = VescTmp/Vdisp;

    auto u = sqrt(u0*u0+ksi*ksi+2*u0*ksi*cosPhi);


    auto sinTheta = (u != 0.0 ? ksi*sinPhi/u : 0);
    auto v = sqrt(u*u+ve*ve);
    auto n = RandomN(G);
    return MC::MCResult<vec3>(n*(v*Vdisp),sinTheta*v*sqrt(M_PI/2));
}


template <class Generator>
/*MK generator of input velocity*/
inline MC::MCResult<vec3> VelocityConstrained(Generator const & G,double VescTmp,double Vmin,double Vmax,
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
inline MC::MCResult<vec3> NuOut(Generator const & G,const vec3& Vcm,const vec3&Nu,
						double Vesc,double mp,double mk,double deltaE = 0){
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
	
	if(cosTh1max <= -1)
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

    return MC::MCResult<vec3>(vNu1,0.5*(1.0+cosTh1max)*Nu1);
}

enum ScatteringType{
	ELASTIC,IONIZATION,MIGDAL
};
enum Target{
    ELECTRON,PROTON
};

template <typename ScatterFuncType,typename VescRFuncType,typename TempRFuncType,typename Generator>
inline MC::MCResult<std::tuple<vec3,double,double>> Vout(double mk,double delta_mk,ScatteringType ST,Target T,
                              ScatterFuncType const & PhiFactor, VescRFuncType const & VescR,
                               TempRFuncType const & TempR, Generator const & G,
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


    double n_nd = 1.0;//TODO n_nd as a function of radius
    vec3 V1(0,0,0);//TODO: add thermal distribution of nuclei velocity

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
    double s = q/(alpha*me);
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

double maxLnd(const BodyModel& ,double );
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
    auto TempR = 0;//BM.UniformRadFunc("Temp");
    double VescMin = BM.VescMin();
	for(size_t i=0;i<Nmk;++i){
        auto mk_res = Vout(mk,delta_mk,ST,T,PhiFactor,VescR,TempR,G,Vdisp,mU0);
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
        */
        sum += mk_res.RemainDensity/Nmk;
	}
    return sum;
	
}

double maxLnd(const BodyModel & BM, double End){
    if(End > -0.5){
        return sqrt(2+2*End);
    }
    std::vector<double> phi = -BM["phi"];

    auto & R = BM["Radius"];

    //auto V = R*vmap(sqrt,2*(End-phi));

    size_t i1 = Function::find_less(phi,End);
    size_t i0 = 0;

    double L0 = 0;
    double L1 = 0;
    while(i0 + 2 < i1){
        size_t i_l = (i0+i1)/2;
        size_t i_r = i_l + 1;
        double L_l = R[i_l]*sqrt(End-phi[i_l]);
        double L_r = R[i_r]*sqrt(End-phi[i_r]);

        if(L_l > L_r){
            i1 = i_r;
            L1 = L_r;
        }
        else{
            i0 = i_l;
            L0 = L_l;
        }

    }
    return std::max(std::max(L0,L1),R[i0+1]*sqrt(End-phi[i0+1]));
}

template <typename T,typename U>
auto dot(const std::vector<T> & X,const std::vector<U> & Y){
    decltype(std::declval<T>()*std::declval<U>()) summ = 0;
    for(size_t i=0;i<X.size();++i){
        summ += X[i]*Y[i];
    }
    return summ;
}

template <typename T,typename U>
auto dot(const std::vector<std::vector<T>> & X,const std::vector<U> & Y){
    std::vector<decltype(std::declval<T>()*std::declval<U>())> summ(X.size(),0);
    for(size_t i=0;i<X.size();++i){
        summ[i] += dot(X[i],Y);
    }
    return summ;
}

template <typename T,typename U>
auto dot(const std::vector<U> & Y,const std::vector<std::vector<T>> & X){
    std::vector<decltype(std::declval<T>()*std::declval<U>())> sum(Y.size(),0);
    for(size_t i=0;i<X.size();++i){
        for(size_t j=0;j<Y.size();++j){
            sum[i] += X[j][i]*Y[j];
        }
    }
    return sum;
}

template <typename  T>
void matrix_dot(std::vector<T> &result,const std::vector<std::vector<T>> &Mat,const std::vector<T> &vec){
    for(size_t i=0;i<result.size();++i){
        result[i] = dot(Mat[i],vec);
    }
}

template <typename T>
std::vector<std::vector<T>> E_Matrix(size_t N){
    return Vector(N,[N](size_t i){
        return Vector(N,[i](size_t j){
            return (T)(i==j);
        });
    });
}

template <typename T>
std::vector<std::vector<T>> MatrixDiagonal(const std::vector<T> &Diag){
    return Vector(Diag.size(),[N = Diag.size(),&Diag](size_t i){
        return Vector(N,[i,Diag](size_t j){
            return (i==j ? Diag[j] : 0);
        });
    });
}

template <typename T>
auto MatrixTranspose(const std::vector<std::vector<T>> &A){
    std::vector<std::vector<T>> Result(A[0].size());
    for(size_t i=0;i<Result.size();++i){
        Result[i].resize(A.size());
        for(size_t j=0;j<A.size();++j){
            Result[i][j] = A[j][i];
        }
    }
    return Result;
}


template <typename T>
void MatrixMult(const std::vector<std::vector<T>> &A,const std::vector<std::vector<T>> &B,std::vector<std::vector<T>> &Result){
    for(size_t i=0;i<A.size();++i){
        for(size_t j=0;j<B[0].size();++j){
            Result[i][j] = 0;
            for(size_t k=0;k<B.size();++k){
                Result[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

template <typename T>
std::vector<std::vector<T>> MatrixMult(const std::vector<std::vector<T>> &A,const std::vector<std::vector<T>> &B){
    auto BT = MatrixTranspose(B);
    std::vector<std::vector<T>> Ret(A.size());
    for(size_t i=0;i<A.size();++i){
        Ret[i].resize(B[0].size());
        for(size_t j=0;j<B[0].size();++j){
            Ret[i][j] = dot(A[i],BT[j]);
        }
    }
    return Ret;
}


template <typename T>
std::vector<std::vector<T>> MatrixPow(const std::vector<std::vector<T>> &A,size_t N){

    if(N == 0){
        return E_Matrix<T>(A.size());
    }
    std::vector<std::vector<T>> Result = E_Matrix<T>(A.size());
    size_t max2pow = 0;
    size_t pow_of2 = 1;
    while(pow_of2*2<=N){
        max2pow++;
        pow_of2*=2;
    }

    Result = A;
    N %= pow_of2;
    pow_of2 /= 2;
    //max2pow --;

    while(max2pow){
        Result = MatrixMult(Result,Result);
        if(N & pow_of2){
            Result = MatrixMult(Result,A);
        }

        N %=pow_of2;
        max2pow --;
        pow_of2 /= 2;
    }
    return Result;
}

template <typename T>
std::vector<std::vector<T>> Rmatrix(std::vector<std::vector<T>> ScatterMatrix,T step){
    size_t N = ScatterMatrix.size();
    return E_Matrix<T>(N) +
            step*(MatrixTranspose(ScatterMatrix) - MatrixDiagonal(Vector(N,[&ScatterMatrix](size_t i){
                                                                      return vector_sum(ScatterMatrix[i]);
                                                                  })));
}







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


struct TrajectoryInfo{
    //x = r, y = v_r, z = v_perp
    Function::GridFunction<double,Function::VectorGrid<double>,Function::LinearInterpolator> Trajectory;
    double T_in,T_out;
};


#define SAME_SIGN(a,b) ( (a > 0 && b > 0) || (a < 0 && b < 0))
template <typename index_func>
inline size_t find_zero_index(const index_func& f,size_t i0,size_t i1){
    auto f0 = f(i0);
    auto f1 = f(i1);
    if(SAME_SIGN(f0,f1)){
        return i0;
    }
    while(i0 + 1 < i1){
        if(f0 == 0){
            return i0;
        }
        if(f1 == 0){
            return i1;
        }
        size_t im = (i0+i1)/2;
        auto fm = f(im);
        if(SAME_SIGN(fm,f0)){
            i0 = im;
            f0 = fm;
        }
        else{
            i1 = im;
            f1 =fm;
        }
    }
    return i0;
}

template <typename FuncType>
inline auto find_two_ranges_zero(const FuncType & F, size_t N){


    size_t i0 = 1;
    double f0 = F(i0);

    size_t i1 = N-1;
    double f1 = F(i1);

    size_t i_mid;
    double f_mid,f_mid_1;
    while(i0 + 2 < i1){
        i_mid = (i0+i1)/2;
        f_mid = F(i_mid);
        if(f_mid > 0){
            return _T(i0,i_mid,i1);
        }
        f_mid_1 = F(i_mid+1);
        if(f_mid_1 > 0){
            return _T(i_mid,i_mid+1,i1);
        }
        if(f_mid_1 > f_mid){
            i0 = i_mid;
        }
        else{
            i1 = i_mid +1;
        }
    }
    return _T(i0,(i0+i1)/2,i1);
}

struct quadric_interpolator{
    double a,b;
    quadric_interpolator(){}
    quadric_interpolator(double x0,double x1,double f0,double f1){
        setAB(x0,x1,f0,f1);
    }
    inline void setAB(double x0,double x1,double f0,double f1){
        double dx2 = x1*x1-x0*x0;
        b = -(f1-f0)/dx2;
        a = (f0*x1*x1-f1*x0*x0)/dx2;
    }
    template <typename...GF_ARGS>
    inline static quadric_interpolator fromX(const Function::GridFunction<GF_ARGS...> & F, double x){
        return quadric_interpolator(F,F.Grid.pos(x));
    }

    template <typename...GF_ARGS>
    inline quadric_interpolator(const Function::GridFunction<GF_ARGS...> & F, size_t ix){
        setAB(F.Grid[ix],F.Grid[ix+1],F.values[ix],F.values[ix+1]);
    }

    inline double operator ()(double x)const{
        return a - b*x*x;
    }
    operator std::tuple<double,double>() const {
        return _T(a,b);
    }
};

//integral of function 1/(a-e-b*x^2-lp^2/x^2)
template <typename T>
inline auto bound(T x,T a,T b){
    //print("bound ",x);
    return ( x > a ? (x < b ? x : b) : a);
}
inline double d_asin(double x,double dx){
    if(x >=1.0)
        return 0;
    if(x <= -1.0){
        if(x + dx<= -1.0)
            return 0;
        else{
            dx += 1.0+x;
            x = -1;
        }
    }
    if(x+dx >= 1){
        dx = 1-x;
    }
    double c_x = sqrt(1.0-x*x);
    double c_2 = sqrt(1.0-(x+dx)*(x+dx));
    return asin( dx*(c_x + x*(2*x+dx)/(c_x+c_2)) );

}
inline double sqrt_integral(double x0,double dx, quadric_interpolator qi,double ep,double lp){
    double sqrt_a_e_bl = sqrt((qi.a-ep)*(qi.a-ep)-4*qi.b*lp*lp);
    double x1 = x0+dx;
    double y = -(qi.a-ep-2*qi.b*x0*x0)/sqrt_a_e_bl;
    double dy = 2*qi.b*dx*(2*x0+dx)/sqrt_a_e_bl;
    //print("qi = ", qi.a,", ",qi.b);
    return (d_asin(y,dy)/(2*sqrt(qi.b)));
    /*
    PVAR((asin(y+dy)-asin(y))/(2*sqrt(qi.b)));
    debug_return( (
           asin( bound((qi.a-ep-2*qi.b*x0*x0)/sqrt_a_e_bl,-1.,1.)) - asin( bound((qi.a-ep-2*qi.b*x1*x1)/sqrt_a_e_bl,-1.,1.))
                )/(2*sqrt(qi.b))
                  );
    */
}

inline auto quad_solve(double a,double b, double c){
    double S1 = -b/(2*a);
    double D = sqrt(b*b-4*a*c)/(2*a);
    if(D>=0)
        return _T(S1-D,S1+D);
    else
        return _T(S1+D,S1-D);
}
inline auto quad_extr(double a,double b){
    return -b/(2*a);
}

inline double interval_function(double x,double a,double b){
    if(x<a)
        return a + b -2*x;
    else if (x<=b)
        return b-a;
    else
        return 2*x-a-b;
}

template <typename...GF_ARGS>
inline double r_zero_precise(const Function::GridFunction<GF_ARGS...> & phi,size_t ix,double ep,double lp){
    quadric_interpolator qi(phi.Grid[ix],phi.Grid[ix+1],phi.values[ix],phi.values[ix+1]);
    double r0,r1;
    _R(r0,r1) = quad_solve(qi.b,ep-qi.a,lp*lp);

    r0 = sqrt(r0);
    r1 = sqrt(r1);
    //print(r0,"\t",r1);
    //PVAR(r0<=phi.Grid[ix+1]);
    if(interval_function(r0,phi.Grid[ix],phi.Grid[ix+1]) <= interval_function(r1,phi.Grid[ix],phi.Grid[ix+1])){
        return r0;
    }
    else{
        return r1;
    }
}

template <typename FuncType>
size_t argmax(const FuncType & F,size_t i0,size_t i1){
    //cout << i0 <<", " <<i1 <<endl;
    auto f0 = F(i0);
    auto f1 = F(i1);
    const double phib = (sqrt(5.0)-1.0)/2.0;
    const double phi1 = (1-phib);
    size_t il = (size_t)(phib*i0+phi1*i1);
    size_t ir = (size_t)(phib*i1+phi1*i0+0.5);
    auto fl =  F(il);
    auto fr =  F(ir);
    while(i0+4<i1){
        //cout << "[" <<i0 <<", "<< il <<", "<< ir<< ", "<< i1 << "]" <<endl;
        //cout << "[" <<f0 <<", "<< fl <<", "<< fr<< ", "<< f1 << "]" <<endl;
        if(fl >= fr){
            i1 = ir;
            f1 = fr;
            ir = il;
            fr = fl;
            il = (size_t)(phib*i0+phi1*i1);
            fl = F(il);
        }else{
            i0 = il;
            f0 = fl;
            il = ir;
            fl = fr;
            ir = (size_t)(phib*i1+phi1*i0+0.5);
            fr =  F(ir);
        }
    }
    auto mx = f0;
    size_t im = i0;
    while(i0<i1){
        ++i0;
        f0 = F(i0);
        if(mx < f0){
            mx = f0;
            im = i0;
        }
    }
    return im;
}

template <typename...GF_ARGS>
inline auto r_zeros_precise(const Function::GridFunction<GF_ARGS...> & phi,size_t ix,double ep,double lp){
    quadric_interpolator qi(phi.Grid[ix],phi.Grid[ix+1],phi.values[ix],phi.values[ix+1]);
    double r0,r1;
    _R(r0,r1) = quad_solve(qi.b,ep-qi.a,lp*lp);

    r0 = sqrt(r0);
    r1 = sqrt(r1);
    if(r0<r1)
        return _T(r0,r1);
    else
        return _T(r1,r0);
}

template <typename...GF_ARGS>
TrajectoryInfo CalculateTrajectory(const Function::GridFunction<GF_ARGS...> & phi,double E_nd,double L_nd, size_t N_bins){
    const double ep = -E_nd;
    const double lp = L_nd;
    TrajectoryInfo TI;

    double sqrt_e= sqrt(ep);
    double sqrt_1_e_l2 =sqrt(1 - ep - lp*lp);




    auto sqrt_function = [&phi,ep,lp](size_t i){
        return phi.values[i] - ep - lp*lp/(phi.Grid[i]*phi.Grid[i]);
    };



    //PVAR(sqrt_function.toString());
    double rmin,rmax;
    _R(rmin,rmax) = quad_solve(ep,-1,L_nd*L_nd);

    //PVAR(rmin);
    //PVAR(rmax);
    if(rmin >= 1){
        TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,0}),std::vector<double>({1.,1.,}));
        TI.T_in = 0;
        TI.T_out = M_PI/(2*sqrt_e);
        //print("all out");
        return TI;
    }
    else if(rmax >= 1 && rmin < 1){
        //PVAR(rmax);
        //PVAR(rmin);
        rmax = 1;
        //print("partly out");
        if(lp == 0){
            rmin = 0;
        }
        else{
            rmin = r_zero_precise(phi,find_zero_index(sqrt_function,0,phi.Grid.size()-1),ep,lp);
        }
        TI.T_out = (M_PI-2*asin((2*ep-1)/sqrt(1-4*ep*lp*lp)))/(4*ep*sqrt_e) + sqrt_1_e_l2/ep;
//        if(std::isnan(rmin) or std::isnan(rmax)){
//            print("rmin: ",rmin,", rmax: ",rmax);
//            PVAR(ep);
//            PVAR(sqrt_function(phi.Grid.size()-1));
//            PVAR(phi.values[phi.Grid.size()-1]);
//        }
    }
    else{
        //print("all in");
        TI.T_out = 0;
        if(lp == 0){
            if(phi[0] - ep >0){
                rmin = 0;
                rmax = r_zero_precise(phi,find_zero_index(sqrt_function,0,phi.num_of_element()-1),ep,lp);
            }
            else{
                TI.T_in = M_PI/(2*sqrt(quadric_interpolator(phi,0).b));
                TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,TI.T_in}),std::vector<double>({0.,0.,}));
                return TI;
            }
        }
        else{
            size_t i0 = 0;
            size_t i1 = phi.num_of_element()-1;
            size_t im = argmax(sqrt_function,i0,i1);
            if(sqrt_function(im) <=0){
                double rt = 1;
                quadric_interpolator qi;
                if(im >= phi.num_of_element()-1){
                    qi= quadric_interpolator(phi,phi.num_of_element()-2);
                }
                else{
                    qi = quadric_interpolator(phi,im);
                    rt =  sqrt(sqrt(lp*lp/qi.b));
                }

                //COMPARE(rt,phi.Grid[im]);
                //PVAR(qi.b);
                TI.T_in = M_PI/(2*sqrt(qi.b));
                TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,TI.T_in}),std::vector<double>({0.,0.,}));
                return TI;
            }
            else{
                rmin = r_zero_precise(phi,find_zero_index(sqrt_function,i0,im),ep,lp);
                rmax = r_zero_precise(phi,find_zero_index(sqrt_function,im,i1),ep,lp);
                if(std::isnan(rmin) or std::isnan(rmax))
                    printd(", ",i0,im,i1,sqrt_function(i0),sqrt_function(im),sqrt_function(i1));
            }
        }

    }
    //PVAR(TI.T_out);
    //print("rmin = ", rmin, "\trmax = ", rmax);
    double h = (rmax-rmin)/(N_bins-1);
    TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>(N_bins),
                                            Vector(N_bins,[rmin,rmax,N_bins](size_t i){
                                                double a = i/(N_bins-1.0);
                                                return rmin*(1-a) + a*rmax;
                                            }));
    TI.Trajectory.Grid[0] = 0;

    double x_ref = -1;
    quadric_interpolator qi;
    for(size_t i=1;i<TI.Trajectory.values.size();++i){
        double x0 = TI.Trajectory.values[i-1];
        size_t ix = phi.Grid.pos(x0);
        if(x_ref != phi.Grid[ix]){
            x_ref = phi.Grid[ix];
            qi = quadric_interpolator(phi,ix);
        }
        //print("qi: ", qi.a, ", ",qi.b);
        TI.Trajectory.Grid[i] = TI.Trajectory.Grid[i-1] + sqrt_integral(x0,h,qi,ep,lp);
    }
    TI.T_in = TI.Trajectory.Grid.grid().back();

    return TI;
}


template <class Generator>
/*MK generator of output nu'*/
inline auto NuOutTherm(Generator const & G,const vec3& Vcm,const vec3&Nu,
                                     double mp,double mk,double deltaE = 0){
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

    return MC::MCResult<vec3>(vNu1,1.0);
}

template <typename ScatterFuncType,typename Generator>
inline MC::MCResult<vec3> VoutTherm(double mk,double delta_mk,ScatteringType ST,Target T,
                              ScatterFuncType const &PhiFactor, const vec3& V_wimp, double Therm, Generator const & G){

    double mp = 0.938;
    double me = 0.52e-3;
    double factor = 1.0;


    vec3 V1 = RandomN(G)*sqrt(Therm*mp/2);//TODO: add thermal distribution of nuclei velocity

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
    double s = q/(alpha*me);
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
    */

    return MC::MCResult<vec3>(Nu1+Vcm,factor);
}


template <typename HistoType,typename Generator,typename ThermFuncType,typename VescFuncType,typename ScatterFuncType>
inline void TrajectoryIntegral(Generator const & G,
                        double e_nd,double l_nd,const TrajectoryInfo & TI,double VescMin,const ThermFuncType & ThermR,
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

        double v = sqrt(e_nd*VescMin*VescMin+vesc*vesc);
        if(std::isnan(v)){
            std::cout << "v is nan:" <<std::endl;
            PVAR(e_nd);
            PVAR(l_nd);
            PVAR(r);
            std::cout << "----------" <<std::endl;
            return;
        }
        vec3 Vin = RandomN(G)*v;
        auto Vout = VoutTherm(mk,delta_mk,ST,T,PhiFactor,Vin,Therm,G);

        double e_nd_1 = (Vout.Result*Vout.Result - vesc*vesc)/(VescMin*VescMin);
        double l_nd_1  = r*sqrt(Vout.Result.x*Vout.Result.x+Vout.Result.y*Vout.Result.y)/VescMin;


        double dens = Vout.RemainDensity/Nmk;
        if(!Out.putValue(dens,e_nd_1,l_nd_1)){
            EvaporationOut += dens;
        }
    }
}
#endif
