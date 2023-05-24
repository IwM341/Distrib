#ifndef TRAJ_H
#define TRAJ_H
#include <utils>
#include <random>
#include <cmath>
#include <tuple>
#include <type_traits>



double maxLnd(const BodyModel & BM, double End){
    if(End >= -0.5){
        return sqrt(1+End);
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

template <typename GridFuncType>
double maxLndf(const GridFuncType& phi, double End){
    if(End >= -0.5){
        return sqrt(1+End);
    }
    //std::vector<double> phi = -BM["phi"];

    //auto & R = BM["Radius"];

    //auto V = R*vmap(sqrt,2*(End-phi));

    size_t i1 = Function::find_more(phi.Values,-End);
    size_t i0 = 0;

    double L0 = 0;
    double L1 = 0;
    while(i0 + 2 < i1){
        size_t i_l = (i0+i1)/2;
        size_t i_r = i_l + 1;
        double L_l = phi.Grid[i_l]*sqrt(End+phi.Values[i_l]);
        double L_r = phi.Grid[i_r]*sqrt(End+phi.Values[i_r]);

        if(L_l > L_r){
            i1 = i_r;
            L1 = L_r;
        }
        else{
            i0 = i_l;
            L0 = L_l;
        }

    }
    return std::max(std::max(L0,L1),phi.Grid[i0+1]*sqrt(End+phi.Values[i0+1]));
}


template <typename ...Args>
double maxLndf(const Function::GridFunction<Args...>& phi, double End){
    if(End >= -0.5){
        return sqrt(1+End);
    }
    //std::vector<double> phi = -BM["phi"];

    //auto & R = BM["Radius"];

    //auto V = R*vmap(sqrt,2*(End-phi));

    size_t i1 = Function::find_more(phi.values,-End);
    size_t i0 = 0;

    double L0 = 0;
    double L1 = 0;
    while(i0 + 2 < i1){
        size_t i_l = (i0+i1)/2;
        size_t i_r = i_l + 1;
        double L_l = phi.Grid[i_l]*sqrt(End+phi.values[i_l]);
        double L_r = phi.Grid[i_r]*sqrt(End+phi.values[i_r]);

        if(L_l > L_r){
            i1 = i_r;
            L1 = L_r;
        }
        else{
            i0 = i_l;
            L0 = L_l;
        }

    }
    return std::max(std::max(L0,L1),phi.Grid[i0+1]*sqrt(End+phi.values[i0+1]));
}



struct TrajectoryInfo{
    //x = r, y = v_r, z = v_perp
    Function::GridFunction<double,Function::VectorGrid<double>,Function::LinearInterpolator> Trajectory;
    double T_in,T_out;
    double E,L;
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

    template <typename GridFuncType>
    inline quadric_interpolator(const GridFuncType & F, size_t ix){
        setAB(F.Grid[ix],F.Grid[ix+1],F.Values[ix],F.Values[ix+1]);
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
    //PVAR(sqrt_a_e_bl);
    //PVAR(y);
    //PVAR( dy );
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

template <typename GridFuncType>
inline double r_zero_precise(const GridFuncType & phi,size_t ix,double ep,double lp){
    quadric_interpolator qi(phi.Grid[ix],phi.Grid[ix+1],phi.Values[ix],phi.Values[ix+1]);
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
    TI.E = E_nd;
    TI.L = L_nd;

    double sqrt_e= sqrt(ep);
    double sqrt_1_e_l2 =sqrt(1 - ep - lp*lp);




    auto sqrt_function = [&phi,ep,lp](size_t i){
        return phi.values[i] - ep - lp*lp/(phi.Grid[i]*phi.Grid[i]+1e-100);
    };



    //PVAR(sqrt_function.toString());
    double rmin,rmax;
    if(ep > 0)
        _R(rmin,rmax) = quad_solve(ep,-1,L_nd*L_nd);
    else
        _R(rmin,rmax) = _T(L_nd*L_nd,1e10);
//    PVAR(rmin);
//    PVAR(rmax);
    if(rmin >= 1){
        TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,0}),std::vector<double>({1.,1.,}));
        TI.T_in = 0;
        TI.T_out = M_PI/(2*ep*sqrt_e);
//        print("all out");
        return TI;
    }
    else if(rmax >= 1 && rmin < 1){
//        PVAR(rmax);
//        PVAR(rmin);
        rmax = 1;
//        print("partly out");
        if(lp == 0){
            rmin = 0;
        }
        else{
            rmin = r_zero_precise(phi,find_zero_index(sqrt_function,0,phi.Grid.size()-1),ep,lp);
        }
        if(ep != 0)
           TI.T_out = (M_PI-2*asin((2*ep-1)/sqrt(1-4*ep*lp*lp)))/(4*ep*sqrt_e) + sqrt_1_e_l2/ep;
        else
           TI.T_out = 1e10;
//        if(std::isnan(rmin) or std::isnan(rmax)){
//            print("rmin: ",rmin,", rmax: ",rmax);
//            PVAR(ep);
//            PVAR(sqrt_function(phi.Grid.size()-1));
//            PVAR(phi.values[phi.Grid.size()-1]);
//        }
//        PVAR(TI.T_out );
    }
    else{
//        print("all in");
        TI.T_out = 0;
        if(lp == 0){
            if(phi[0] - ep >0){
                rmin = 0;
                rmax = r_zero_precise(phi,find_zero_index(sqrt_function,0,phi.num_of_element()-1),ep,lp);
                size_t i =(find_zero_index(sqrt_function,0,phi.num_of_element()-1));
                //COMPARE(sqrt_function(0),sqrt_function(1));
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

//                COMPARE(rt,phi.Grid[im]);
//                PVAR(qi.b);

                TI.T_in = M_PI/(2*sqrt(qi.b));
                TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,TI.T_in}),std::vector<double>({phi.Grid[im],phi.Grid[im],}));
                return TI;
            }
            else{
                rmin = r_zero_precise(phi,find_zero_index(sqrt_function,i0,im),ep,lp);
                rmax = r_zero_precise(phi,find_zero_index(sqrt_function,im,i1),ep,lp);
//                if(std::isnan(rmin) or std::isnan(rmax))
//                    printd(", ",i0,im,i1,sqrt_function(i0),sqrt_function(im),sqrt_function(i1));
            }
        }

    }
//    PVAR(TI.T_out);
//    print("rmin = ", rmin, "\trmax = ", rmax);
    double h = (rmax-rmin)/(N_bins-1);
    TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>(N_bins),
                                            Vector(N_bins,[rmin,rmax,N_bins](size_t i){
                                                double a = i/(N_bins-1.0);
                                                return rmin*(1-a) + a*rmax;
                                            }));
    TI.Trajectory.Grid[0] = 0;

    double x_ref = -1;
    //quadric_interpolator qi;
    quadric_interpolator qi_s;
    for(size_t i=1;i<TI.Trajectory.values.size();++i){
        double x0 = TI.Trajectory.values[i-1];

        size_t ix = phi.Grid.pos(x0);
        size_t ix_1 =  std::max(ix+1,phi.Grid.pos(TI.Trajectory.values[i]));
        double x_0 = phi.Grid[ix];
        double x_1 = phi.Grid[ix_1];

        //qi = quadric_interpolator(phi,ix);
        qi_s.setAB(x_0,x_1,phi[ix],phi[ix_1]);
//            COMPARE(qi.a,qi_s.a);
//        if(i==5){
//        printd(", ",TI.Trajectory.values[i-1],TI.Trajectory.values[i],ix,ix_1,x_0,x_1,phi[ix],phi[ix_1]);
//            COMPARE(qi.b,qi_s.b);
//        printd(", ",qi_s.a,qi_s.b);
//        print("intergeral from ",x0," to ",x0+h);
//        PVAR(sqrt_integral(x0,h,qi_s,ep,lp));
//        }
//        if(qi.b > 100){
//            PVAR(x0);
//            PVAR(phi.Grid.size());
//            PVAR(phi.Grid[1]);
//            printd(", ",phi[ix],phi[ix]-phi[ix+1]);
//            PVAR(x_ref);
//        }
        //print("qi: ", qi.a, ", ",qi.b);
        TI.Trajectory.Grid[i] = TI.Trajectory.Grid[i-1] + sqrt_integral(x0,h,qi_s,ep,lp);
    }
    TI.T_in = TI.Trajectory.Grid.grid().back();

    return TI;
}


struct TrajectoryPreInfo{
    double rmin,rmax;
    double T;

    inline operator double()const noexcept{return T;}
};

template <typename PhiFunctype>
TrajectoryPreInfo CalculatePeriod(const PhiFunctype & phi,double E_nd,double L_nd, size_t N_bins){
    const double ep = -E_nd;
    const double lp = L_nd;
    //TrajectoryInfo TI;
    //TI.E = E_nd;
    //TI.L = L_nd;

    double sqrt_e= sqrt(ep);
    double sqrt_1_e_l2 =sqrt(1 - ep - lp*lp);


    auto sqrt_function = [&phi,ep,lp](size_t i){
        return phi.Values[i] - ep - lp*lp/(phi.Grid[i]*phi.Grid[i]+1e-100);
    };


    double T_in = 0,T_out = 0;
    //PVAR(sqrt_function.toString());
    double rmin,rmax,r_max1;
    if(ep > 0)
        _R(rmin,rmax) = quad_solve(ep,-1,L_nd*L_nd);
    else
        _R(rmin,rmax) = _T(L_nd*L_nd,1e10);
//    PVAR(rmin);
//    PVAR(rmax);
    r_max1 = rmax;
    if(rmin >= 1){
        //TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,0}),std::vector<double>({1.,1.,}));
//        print("all out");
        return {rmin,rmax,M_PI/(2*ep*sqrt_e)};
    }
    else if(rmax >= 1 && rmin < 1){
        r_max1 = 1;
        if(lp == 0){
            rmin = 0;
        }
        else{
            rmin = r_zero_precise(phi,find_zero_index(sqrt_function,0,phi.Grid.size()-1),ep,lp);
        }
        if(ep != 0)
           T_out = (M_PI-2*asin((2*ep-1)/sqrt(1-4*ep*lp*lp)))/(4*ep*sqrt_e) + sqrt_1_e_l2/ep;
        else
           T_out = 1e10;
//        if(std::isnan(rmin) or std::isnan(rmax)){
//            print("rmin: ",rmin,", rmax: ",rmax);
//            PVAR(ep);
//            PVAR(sqrt_function(phi.Grid.size()-1));
//            PVAR(phi.values[phi.Grid.size()-1]);
//        }
//        PVAR(TI.T_out );
    }
    else{
//        print("all in");
        T_out = 0;
        if(lp == 0){
            if(phi[0] - ep >0){
                rmin = 0;
                rmax = r_zero_precise(phi,find_zero_index(sqrt_function,0,phi.Grid.size()-1),ep,lp);
                size_t i =(find_zero_index(sqrt_function,0,phi.Grid.size()-1));
                //COMPARE(sqrt_function(0),sqrt_function(1));
            }
            else{
                T_in = M_PI/(2*sqrt(quadric_interpolator(phi,0).b));
                //Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,TI.T_in}),std::vector<double>({0.,0.,}));
                return {rmin,rmax,T_in+T_out};
            }
        }
        else{
            size_t i0 = 0;
            size_t i1 = phi.Grid.size()-1;
            size_t im = argmax(sqrt_function,i0,i1);
            if(sqrt_function(im) <=0){
                double rt = 1;
                quadric_interpolator qi;
                if(im >= phi.Grid.size()-1){
                    qi= quadric_interpolator(phi,phi.Grid.size()-2);
                }
                else{
                    qi = quadric_interpolator(phi,im);
                    rt =  sqrt(sqrt(lp*lp/qi.b));
                }

//                COMPARE(rt,phi.Grid[im]);
//                PVAR(qi.b);

                T_in = M_PI/(2*sqrt(qi.b));
                //TI.Trajectory = decltype(TI.Trajectory)(std::vector<double>({0,TI.T_in}),std::vector<double>({phi.Grid[im],phi.Grid[im],}));
                return {rmin,rmax,T_in+T_out};
            }
            else{
                rmin = r_zero_precise(phi,find_zero_index(sqrt_function,i0,im),ep,lp);
                rmax = r_zero_precise(phi,find_zero_index(sqrt_function,im,i1),ep,lp);
//                if(std::isnan(rmin) or std::isnan(rmax))
//                    printd(", ",i0,im,i1,sqrt_function(i0),sqrt_function(im),sqrt_function(i1));
            }
        }
        r_max1 = rmax;
    }
//    PVAR(TI.T_out);
//    print("rmin = ", rmin, "\trmax = ", rmax);
    double h = (r_max1-rmin)/(N_bins-1);

    auto _x = [N_bins,rmin,r_max1](size_t i){
        double a = i/(N_bins-1.0);return rmin*(1-a) + a*r_max1;
    };

    double x_ref = -1;
    //quadric_interpolator qi;
    quadric_interpolator qi_s;
    for(size_t i=1;i<N_bins;++i){
        double x0 = _x(i-1);

        size_t ix = phi.Grid.pos(x0);
        size_t ix_1 =  std::max(ix+1,phi.Grid.pos(_x(i)));
        double x_0 = phi.Grid[ix];
        double x_1 = phi.Grid[ix_1];

        //qi = quadric_interpolator(phi,ix);
        qi_s.setAB(x_0,x_1,phi[ix],phi[ix_1]);
//            COMPARE(qi.a,qi_s.a);
//        if(i==5){
//        printd(", ",TI.Trajectory.values[i-1],TI.Trajectory.values[i],ix,ix_1,x_0,x_1,phi[ix],phi[ix_1]);
//            COMPARE(qi.b,qi_s.b);
//        printd(", ",qi_s.a,qi_s.b);
//        print("intergeral from ",x0," to ",x0+h);
//        PVAR(sqrt_integral(x0,h,qi_s,ep,lp));
//        }
//        if(qi.b > 100){
//            PVAR(x0);
//            PVAR(phi.Grid.size());
//            PVAR(phi.Grid[1]);
//            printd(", ",phi[ix],phi[ix]-phi[ix+1]);
//            PVAR(x_ref);
//        }
        //print("qi: ", qi.a, ", ",qi.b);
        T_in += sqrt_integral(x0,h,qi_s,ep,lp);
    }

    return {rmin,rmax,T_in+T_out};
}


#endif
