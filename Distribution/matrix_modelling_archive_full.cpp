
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>

#include "func/matrix_functions.hpp"
#include <array>
#include <regex>
#include <string>
#include <string>


#include "func/arg_parser.hpp"
#include <sstream>
#include <boost/filesystem.hpp>
#include <grid_objects.hpp>
#include <csv_io.hpp>
#include "func/cbals_functions.hpp"

typedef  double vtype;
#define vgemm cblas_dgemm




template <typename T>
struct ResultData{
    Function::GridFunction<T,Function::VectorGrid<T>,Function::LinearInterpolator> n_l;
    decltype(n_l) n_h;
    decltype(n_l) N_L,N_H;
    std::vector<T> d_l,d_h,D_L,D_H;
};

template <typename T>
ResultData<T> MatrixMultMethod(size_t NH,size_t NL,std::vector<T> &A_LL,std::vector<T> &A_LH_T,
                               std::vector<T> &A_HL_T,std::vector<T> &A_HH,
                               std::vector<T> &D_L_in,std::vector<T> &D_H_in,size_t log2_n,T tau){
    ResultData<T>  RD;
    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;

    std::vector<T>B_HL_T(A_HL_T.size());
    std::vector<T>B_LH_T(A_LH_T.size());
    std::vector<T>B_LL(A_LL.size());
    std::vector<T>B_HH(A_HH.size());

    RD.n_l = decltype(RD.n_l)(Vector(log2_n+1,[tau](size_t i)->T{return (i==0 ? 0 :tau*pow(2.0,i-1));}));
    RD.N_H = RD.N_L= RD.n_h = RD.n_l;

    RD.n_l[0] = vector_sum(D_L_in);
    RD.n_h[0] = vector_sum(D_H_in);
    RD.N_H[0] = 0;
    RD.N_L[0] = 0;

    RD.d_l = D_L_in;
    RD.d_h = D_H_in;
    auto LD_last = D_L_in;
    auto HD_last = D_H_in;

    RD.D_L = D_L_in*0;
    RD.D_H = D_H_in*0;
    print("multiply matricies");


    for(size_t i=0;i<=log2_n;++i){
        std::cout <<"step: "<< i <<std::endl;
        std::swap(RD.d_l,LD_last);
        std::swap(RD.d_h,HD_last);
        if(i>0){
            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NL,
                        1.0,A_LL.data(),ld_A_LL,A_LL.data(),ld_A_LL,0.0,B_LL.data(),ld_A_LL);
            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NH,
                                1.0,A_LH_T.data(),ld_A_LH,A_HL_T.data(),ld_A_HL,1.0,B_LL.data(),ld_A_LL);
            //DGEMM()

            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NL,
                        1.0,A_LL.data(),ld_A_LL,A_LH_T.data(),ld_A_LH,0.0,B_LH_T.data(),ld_A_LH);
            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NH,
                                1.0,A_LH_T.data(),ld_A_LH,A_HH.data(),ld_A_HH,1.0,B_LH_T.data(),ld_A_LH);


            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NL,
                        1.0,A_HL_T.data(),ld_A_HL,A_LL.data(),ld_A_LL,0.0,B_HL_T.data(),ld_A_HL);
            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NH,
                                1.0,A_HH.data(),ld_A_HH,A_HL_T.data(),ld_A_HL,1.0,B_HL_T.data(),ld_A_HL);

            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NL,
                        1.0,A_HL_T.data(),ld_A_HL,A_LH_T.data(),ld_A_LH,0.0,B_HH.data(),ld_A_HH);
            gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NH,
                                1.0,A_HH.data(),ld_A_HH,A_HH.data(),ld_A_HH,1.0,B_HH.data(),ld_A_HH);


            std::swap(A_LL,B_LL);
            std::swap(A_HL_T,B_HL_T);
            std::swap(A_LH_T,B_LH_T);
            std::swap(A_HH,B_HH);
        }
        //printPairMatrix(NL,NH,A_LL.data(),A_HL_T.data(),A_LH_T.data(),A_HH.data());
        if(D_L_in.size() && D_H_in.size() && RD.d_l.size() && RD.d_h.size()){
            gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,A_LL.data(),ld_A_LL,D_L_in.data(),1,0.0,RD.d_l.data(),1);
            gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,A_LH_T.data(),ld_A_LH,D_H_in.data(),1,1.0,RD.d_l.data(),1);

            gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,A_HL_T.data(),ld_A_HL,D_L_in.data(),1,0.0,RD.d_h.data(),1);
            gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,A_HH.data(),ld_A_HH,D_H_in.data(),1,1.0,RD.d_h.data(),1);
        }

        double dt = RD.n_l.Grid[i]-RD.n_l.Grid[i-1];
        RD.n_l[i] = vector_sum(RD.d_l);
        RD.n_h[i] = vector_sum(RD.d_h);
        RD.N_L[i] = RD.N_L[i-1] +0.5*(RD.n_l[i]+RD.n_l[i-1])*dt;
        RD.N_H[i] = RD.N_H[i-1] +0.5*(RD.n_h[i]+RD.n_h[i-1])*dt;
        RD.D_L = RD.D_L+ T(0.5*dt)*(RD.d_l+LD_last);
        RD.D_H = RD.D_H + T(0.5*dt)*(RD.d_h+HD_last);

    }
    return RD;
}
template <typename T>
ResultData<T> VectorMultMethod(size_t NH,size_t NL,std::vector<T> &A_LL,std::vector<T> &A_LH_T,
                               std::vector<T> &A_HL_T,std::vector<T> &A_HH,
                               std::vector<T> &D_L_in,std::vector<T> &D_H_in,size_t n,T tau){
    ResultData<T>  RD;
    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;

    RD.n_l = decltype(RD.n_l)(Vector(n+1,[tau](size_t i)->T{return i*tau;}));
    RD.N_H = RD.N_L= RD.n_h = RD.n_l;

    RD.n_l[0] = vector_sum(D_L_in);
    RD.n_h[0] = vector_sum(D_H_in);
    RD.N_H[0] = 0;
    RD.N_L[0] = 0;

    RD.d_l = D_L_in;
    RD.d_h = D_H_in;
    auto LD_last = D_L_in;
    auto HD_last = D_H_in;

    RD.D_L = D_L_in*0;
    RD.D_H = D_H_in*0;


    show_prog prog_lh(100);
    prog_lh.show(0);
    size_t total_progress_lh = n;
    size_t curr_progress_lh = 0;
    prog_lh.show(curr_progress_lh/((float)total_progress_lh));
    for(size_t i=1;i<=n;++i){
        std::swap(RD.d_l,LD_last);
        std::swap(RD.d_h,HD_last);
        //printPairMatrix(NL,NH,A_LL.data(),A_HL_T.data(),A_LH_T.data(),A_HH.data());
        if(D_L_in.size() && D_H_in.size() && RD.d_l.size() && RD.d_h.size()){
            gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,A_LL.data(),ld_A_LL,LD_last.data(),1,0.0,RD.d_l.data(),1);
            gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,A_LH_T.data(),ld_A_LH,HD_last.data(),1,1.0,RD.d_l.data(),1);

            gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,A_HL_T.data(),ld_A_HL,LD_last.data(),1,0.0,RD.d_h.data(),1);
            gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,A_HH.data(),ld_A_HH,HD_last.data(),1,1.0,RD.d_h.data(),1);
        }

        double dt = RD.n_l.Grid[i]-RD.n_l.Grid[i-1];
        RD.n_l[i] = vector_sum(RD.d_l);
        RD.n_h[i] = vector_sum(RD.d_h);
        RD.N_L[i] = RD.N_L[i-1] +0.5*(RD.n_l[i]+RD.n_l[i-1])*dt;
        RD.N_H[i] = RD.N_H[i-1] +0.5*(RD.n_h[i]+RD.n_h[i-1])*dt;
        RD.D_L = RD.D_L+ T(0.5*dt)*(RD.d_l+LD_last);
        RD.D_H = RD.D_H + T(0.5*dt)*(RD.d_h+HD_last);

        curr_progress_lh++;
        prog_lh.show(curr_progress_lh/((float)total_progress_lh));
    }
    prog_lh.end();

    return RD;
}

template <typename T>
struct ResultLH{
    grob::GridObject<grob::GridVector<T>,std::vector<T>> N_L;
    decltype(N_L) N_H;
    decltype(N_L) n_l,n_h;
    std::vector<T> C_L_final;
    std::vector<T> C_H_final;
    std::vector<T> D_L_final;
    std::vector<T> D_H_final;
    std::vector<std::vector<T>> D_L_tmp;
    std::vector<std::vector<T>> D_H_tmp;
    T S_LL,S_LH,S_HL,S_HH;
};

///
/// \brief
/// \param R function, R(c_l,c_h,c_l1,c_h1) writes to c_l1 c_h1
/// result of application R to c_l,c_h
/// \param Summator function, should give Summator(C_L,C_L_buffer,D_L_acc) ->
/// D_L_acc += interpolation(C_L,C_L_buffer)
/// \param DiffStepFunctype(D_L,D_H) -> {S_LL,S_LH,S_HL,S_HH}
template <typename T,typename OperatorFunctype,
          typename SummatorFunctype,
          typename DiffStepFunctype>
ResultLH<T> RunScheme(OperatorFunctype &&R,
                    SummatorFunctype &&Summator,
                    DiffStepFunctype && DStep,
                    std::vector<T> & C_L,std::vector<T> & C_H,
                    std::vector<T> & C_L_buffer,std::vector<T> & C_H_buffer,
                    T tau,T T_full,std::vector<T> view_moments = {}){
    size_t N_steps = T_full/tau+0.5;


    std::vector<T> D_L_acc(C_L.size(),0);
    std::vector<T> D_H_acc(C_H.size(),0);

    auto moments = Vector(N_steps+1,[=](size_t i){return i*tau;});
    std::vector<T> N_L(N_steps+1);
    N_L[0] = 0;
    std::vector<T> N_H(N_steps+1);
    N_H[0] = 0;
    std::vector<T> n_l(N_steps+1);
    n_l[0] = vector_sum(C_L);
    std::vector<T> n_h(N_steps+1);
    n_h[0] = vector_sum(C_H);

    progress_bar PB(N_steps,100);
    size_t view_count = 0;
    std::vector<std::vector<T>> N_L_view,N_H_view;
    for(size_t i=1;i<=N_steps;++i){
        R(C_L,C_H,C_L_buffer,C_H_buffer);
        Summator(C_L,C_L_buffer,tau,D_L_acc);
        Summator(C_H,C_H_buffer,tau,D_H_acc);
        std::swap(C_L,C_L_buffer);
        std::swap(C_H,C_H_buffer);

        n_l[i] = vector_sum(C_L);
        n_h[i] = vector_sum(C_H);

        T N_L_tmp = vector_sum(D_L_acc);
        T N_H_tmp = vector_sum(D_H_acc);



        N_L[i] = N_L_tmp;
        N_H[i] = N_H_tmp;
        PB++;
        if(view_count && i*tau <= view_moments[view_count] &&
                view_moments[view_count] < tau*(i+1)){
            N_L_view.push_back(D_L_acc/N_L_tmp);
            N_H_view.push_back(D_H_acc/N_H_tmp);
            view_count ++;
        }
    }
    auto [S_LL,S_LH,S_HL,S_HH] = DStep(D_L_acc,D_H_acc);
    D_L_acc /= N_L.back();
    D_H_acc /= N_H.back();
    C_L /= vector_sum(C_L);
    C_H /= vector_sum(C_H);
    return ResultLH<T>{
        decltype(ResultLH<T>::N_L)(moments,std::move(N_L)),
        decltype(ResultLH<T>::N_H)(moments,std::move(N_H)),
        decltype(ResultLH<T>::N_L)(moments,std::move(n_l)),
        decltype(ResultLH<T>::N_L)(moments,std::move(n_h)),
        std::move(C_L),
        std::move(C_H),
        std::move(D_L_acc),
        std::move(D_H_acc),
        std::move(N_L_view),
        std::move(N_H_view),
        S_LL,S_LH,S_HL,S_HH
    };
}



constexpr int NO_INC = 0;
template <typename T>
ResultLH<T> SchemeEuler(std::vector<T> &A_LL,std::vector<T> & A_LH,std::vector<T> & A_HL,std::vector<T> & A_HH,
                      std::vector<T> & C_L,std::vector<T> & C_H,
                      T tau,T T_full,size_t skip_pow = 0, std::vector<T> view_moments = {}){
    size_t NL = A_LL.size();
    size_t NH = A_HH.size();

    A_LL *= tau;
    A_LH *= tau;
    A_HL *= tau;
    A_HH *= tau;
    C_L *= tau;
    C_H *= tau;


    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;
    size_t N_steps = T_full/tau+0.5;
    std::vector<T> D_L_acc(NL,0);
    std::vector<T> C_L_tmp(NL,0);
    std::vector<T> D_H_acc(NH,0);
    std::vector<T> C_H_tmp(NH,0);

    std::vector<std::vector<T>> N_L_view;
    std::vector<std::vector<T>> N_H_view;

    auto moments = Vector(N_steps+1,[=](size_t i){return i*tau;});
    std::vector<T> N_L(N_steps+1);
    std::vector<T> N_H(N_steps+1);
    std::vector<T> n_l(N_steps+1);
    std::vector<T> n_h(N_steps+1);
    N_L[0] = 0;
    N_H[0] = 0;
    n_l[0] = vector_sum(C_L);
    n_h[0] = vector_sum(C_H);

    progress_bar PB(N_steps,100);
    //prog_lh.show(curr_progress_lh/((float)total_progress_lh));

    size_t view_cout = 0;
    if(NL*NH < 100){
        print("S matrix = ");
        printPairMatrixD(NL,NH,A_LL.data(),A_LH.data(),A_HL.data(),A_HH.data());
        print();
    }
    for(size_t i=1;i<=N_steps;++i){
        //gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,A_LL.data(),ld_A_LL,D_L0_tmp.data(),1,0.0,D_L1_tmp.data(),1);
        //cblas_sgbmv(CblasColMajor,CblasNoTrans,NL,NL,0,0,1,A_LL.data(),ld_A_LL,C_L.size(),NO_INC,1,D_L1_tmp.data(),NO_INC);
        //print(i);
        for(size_t _i=0;_i<NL;++_i){
            C_L_tmp[_i] = C_L[_i]+A_LL[_i]*C_L[_i];
        }
        for(size_t _i=0;_i<NH;++_i){
            C_H_tmp[_i] = C_H[_i]+A_HH[_i]*C_H[_i];
        }
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,A_LH.data(),ld_A_LH,C_H.data(),1,1.0,C_L_tmp.data(),1);
        D_L_acc += C_L;
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,A_HL.data(),ld_A_HL,C_L.data(),1,1.0,C_H_tmp.data(),1);
        //gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,A_HH.data(),ld_A_HH,D_H0_tmp.data(),1,1.0,D_H1_tmp.data(),1);
        D_H_acc += C_H;
        std::swap(C_L,C_L_tmp);
        std::swap(C_H,C_H_tmp);

        n_l[i] = vector_sum(C_L);
        n_h[i] = vector_sum(C_H);

        T N_L_tmp = vector_sum(D_L_acc);
        T N_H_tmp = vector_sum(D_H_acc);

        PB++;

        N_L[i] = N_L_tmp;
        N_H[i] = N_H_tmp;

        if(view_cout && i*tau <= view_moments[view_cout] &&
                view_moments[view_cout] < tau*(i+1)){
            N_L_view.push_back(D_L_acc/N_L_tmp);
            N_H_view.push_back(D_H_acc/N_H_tmp);
            view_cout ++;
        }

    }
    D_L_acc /= N_L.back();
    D_H_acc /= N_H.back();

    auto D_L1_tmp = Vector(NL,[&](size_t i){
            return (1-A_LL[i])*D_L_acc[i];}
            );

    auto D_H1_tmp = Vector(NH,[&](size_t i){
            return (1-A_HH[i])*D_H_acc[i];}
            );
    T S_LL = vector_sum(D_L1_tmp)/tau;
    T S_HH = vector_sum(D_H1_tmp)/tau;
    //cblas_dgemv()
    gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,A_LH.data(),ld_A_LH,D_H_acc.data(),1,0.0,D_L1_tmp.data(),1);
    T S_LH = vector_sum(D_L1_tmp)/tau;

    gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,A_HL.data(),ld_A_HL,D_L_acc.data(),1,0.0,D_H1_tmp.data(),1);
    T S_HL = vector_sum(D_H1_tmp)/tau;


    C_L/=vector_sum(C_L);
    C_H/=vector_sum(C_H);
    return ResultLH<T>{
        decltype(ResultLH<T>::N_L)(moments,std::move(N_L)),
        decltype(ResultLH<T>::N_H)(moments,std::move(N_H)),
        decltype(ResultLH<T>::N_L)(moments,std::move(n_l)),
        decltype(ResultLH<T>::N_L)(moments,std::move(n_h)),
        std::move(C_L),
        std::move(C_H),
        std::move(D_L_acc),
        std::move(D_H_acc),
        std::move(N_L_view),
        std::move(N_H_view),
        S_LL,S_LH,S_HL,S_HH
    };
}


template <typename T>
ResultLH<T> SchemeEulerMult(std::vector<T> &A_LL,std::vector<T> & A_LH,
                          std::vector<T> & A_HL,std::vector<T> & A_HH,
                          std::vector<T> & C_L,std::vector<T> & C_H,
                          T tau,T T_full,size_t skip_pow = 0,std::vector<T> view_moments = {}){

    size_t NL = A_LL.size();
    size_t NH = A_HH.size();
    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;

    auto tau_eff = tau/(1<<skip_pow);
    std::vector<T> R_LL(NL*NL,0);
    std::vector<T> R_LH = A_LH*tau_eff;
    A_LH.clear();
    std::vector<T> R_HL = A_HL*tau_eff;
    A_HL.clear();
    std::vector<T> R_HH(NH*NH,0);
    for(size_t i=0;i<NL;++i){
        R_LL[NL*i+i] = 1 + tau_eff*A_LL[i];
    }
    for(size_t i=0;i<NH;++i){
        R_HH[NH*i+i] = 1 + tau_eff*A_HH[i];
    }

    std::vector<T> D_L0_tmp = C_L;
    std::vector<T> D_L1_tmp(NL);
    std::vector<T> D_H0_tmp = C_H;
    std::vector<T> D_H1_tmp(NH);

    PowMatrix(NL,NH,R_LL,R_LH,R_HL,R_HH,skip_pow);
    print("R = ");
    printPairMatrix(NL,NH,R_LL.data(),R_LH.data(),R_HL.data(),R_HH.data());
    print("Run Scheme");
    return RunScheme(
                [&](auto & C_L,auto & C_H,auto & C_L_new,auto & C_H_new){
                gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,R_LL.data(),ld_A_LL,C_L.data(),1,0.0,C_L_new.data(),1);
                gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,R_LH.data(),ld_A_LH,C_H.data(),1,1.0,C_L_new.data(),1);
                gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,R_HL.data(),ld_A_HL,C_L.data(),1,0.0,C_H_new.data(),1);
                gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,R_HH.data(),ld_A_HH,C_H.data(),1,1.0,C_H_new.data(),1);
                },
        [](auto & D1,auto & D2,T tau,auto & Result){
            for(size_t i=0;i<D1.size();++i){
                Result[i] += 0.5*tau*(D1[i]+D2[i]);
            }
    },[&](auto &,auto &){
        T ret = 0;
        return std::make_tuple(ret,ret,ret,ret);
    },D_L0_tmp,D_H0_tmp,D_L1_tmp,D_H1_tmp,tau,T_full,view_moments);
}


template <typename T>
ResultLH<T> Scheme2Order(std::vector<T> &A_LL,std::vector<T> & A_LH,std::vector<T> & A_HL,std::vector<T> & A_HH,
                      std::vector<T> & C_L,std::vector<T> & C_H,
                      T tau,T T_full,size_t skip_pow = 0,std::vector<T> view_moments = {}){
    size_t NL = A_LL.size();
    size_t NH = A_HH.size();
    size_t N_steps = T_full/tau+0.5;

    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;
    auto tau_eff = tau/(1<<skip_pow);

    std::vector<std::vector<T>> N_L_view;
    std::vector<std::vector<T>> N_H_view;

    auto moments = Vector(N_steps+1,[=](size_t i){return i*tau;});
    std::vector<T> N_L(N_steps+1);
    std::vector<T> N_H(N_steps+1);
    N_L[0] = 0;
    N_H[0] = 0;



    //assuming col major
    auto MultMatrix_XD = [](size_t N,size_t M,
            auto const & X,auto const & Mat_Diag,auto & result){
        for(size_t j=0;j<M;++j){
            for(size_t i=0;i<N;++i){
                result[j*N+i]=X[j*N+i]*Mat_Diag[j];
            }
        }
    };
    auto MultMatrix_DX = [](size_t N,size_t M,
            auto const & Mat_Diag,auto const & X,auto & result){
        for(size_t j=0;j<M;++j){
            for(size_t i=0;i<N;++i){
                result[j*N+i]=Mat_Diag[i]*X[j*N+i];
            }
        }
    };
    auto SummMatrix_XD = [](size_t N,
            auto & X,auto const & Mat_Diag,T matCoeff){
        for(size_t j=0;j<N;++j){
            X[j*N+j]=X[j*N+j] + matCoeff*Mat_Diag[j];
        }
    };

    auto MatrixMult_XY = [](size_t N,size_t M,size_t K,
            auto const & X,auto const & Y,auto & result,
            T alpha = 1,T beta = 0){
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,N,K,M,
                    alpha,X.data(),N,Y.data(),M,beta,result.data(),N);
    };
    auto MatrixMult_XY_add = [](size_t N,size_t M,size_t K,
            auto const & X,auto const & Y,auto & result){
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,N,K,M,
                    1.0,X.data(),N,Y.data(),M,1.0,result.data(),N);
    };

    A_LL *= (0.5*tau_eff);
    A_LH *= (0.5*tau_eff);
    A_HL *= (0.5*tau_eff);
    A_HH *= (0.5*tau_eff);

    print("geting A^-1, D^-1, X and Y");
    std::vector<T> A_inv_mat = Vector(NL,[&](size_t i)->T{return 1/(1-A_LL[i]);});
    std::vector<T> D_inv_mat = Vector(NH,[&](size_t i)->T{return 1/(1-A_HH[i]);});

    std::vector<T> X(NL*NH);
    MultMatrix_XD(NL,NH,A_LH,D_inv_mat,X);
    X*=(-1);



    std::vector<T> Y(NH*NL);
    MultMatrix_XD(NH,NL,A_HL,A_inv_mat,Y);
    Y*=(-1);

    /**/
    //print("<A^i,X,Y,D^i>");
    //printPairMatrixD(NL,NH,A_inv_mat.data(),X.data(),Y.data(),D_inv_mat.data());
    /**/

    //************************//


    std::vector<T> XY_1(NL*NL,0);

    for(size_t i=0;i<NL;++i){
        XY_1[i*NL+i] = 1;
    }
    MatrixMult_XY(NL,NH,NL,X,Y,XY_1,-1,1);


    std::vector<T> YX_1(NH*NH,0);

    for(size_t i=0;i<NH;++i){
        YX_1[i*NH+i] = 1;
    }
    MatrixMult_XY(NH,NL,NH,Y,X,YX_1,-1,1);

    //print("print 1-XY");
    //printMatrix(NL,NL,XY_1.data());
    //print("print 1-YX");
    //printMatrix(NH,NH,YX_1.data());

    print("inverting 1-XY and 1-YX");
    std::vector<T> XY_1_inv = std::move(XY_1);
    std::vector<T> YX_1_inv = std::move(YX_1);

    std::vector<int32_t> ipiv_L(NL+1);
    auto ret_L = getrf<T>(CblasColMajor,NL,NL,XY_1_inv.data(),NL,ipiv_L.data());
    if(ret_L!=0){
        throw std::runtime_error("error decompoeing L*U matrix, ret_L = " + std::to_string(ret_L));
    }
    getri<T>(CblasColMajor,NL,XY_1_inv.data(),NL,ipiv_L.data());
    if(ret_L!=0){
        throw std::runtime_error("error inverting matrix, ret_L = " + std::to_string(ret_L));
    }

    std::vector<int32_t> ipiv_H(NH+1);
    auto ret_H = getrf<T>(CblasColMajor,NH,NH,YX_1_inv.data(),NH,ipiv_H.data());
    if(ret_H!=0){
        throw std::runtime_error("error decomposing L*U matrix, ret_L = " + std::to_string(ret_H));
    }
    getri<T>(CblasColMajor,NH,YX_1_inv.data(),NH,ipiv_H.data());
    if(ret_H!=0){
        throw std::runtime_error("error inverting matrix, ret_L = " + std::to_string(ret_H));
    }

    //print("print (1-XY)^i");
    //printMatrix(NL,NL,XY_1_inv.data());
    //print("print (1-YX)^i");
    //printMatrix(NH,NH,YX_1_inv.data());

    //*******************//
    //print("getting pre matricies");
    //MultMatrix_DX(NL,NL,A_inv_mat,XY_1_inv,XY_1_inv);
    //MultMatrix_DX(NH,NH,D_inv_mat,YX_1_inv,YX_1_inv);

    std::vector<T> R_LL(NL*NL);
    std::vector<T> R_LH(NL*NH);
    std::vector<T> R_HL(NL*NH);
    std::vector<T> R_HH(NH*NH);


    A_LL += 1;
    A_HH += 1;

    print("getting R_LL");
    std::vector<T> __LL(NL*NL);
    MatrixMult_XY(NL,NH,NL,X,A_HL,__LL,-1,0);
    SummMatrix_XD(NL,__LL,A_LL,1);
    MatrixMult_XY(NL,NL,NL,XY_1_inv,__LL,R_LL);
    MultMatrix_DX(NL,NL,A_inv_mat,R_LL,R_LL);
    __LL.resize(0);

    print("getting R_HL");
    for(size_t i=0;i<NL;++i){
        for(size_t j=0;j<NH;++j){
            A_HL[i*NH+j] -= Y[i*NH+j]*A_LL[i];
        }
    }
    MatrixMult_XY(NH,NH,NL,YX_1_inv,A_HL,R_HL);
    MultMatrix_DX(NH,NL,D_inv_mat,R_HL,R_HL);

    print("getting R_HH");
    std::vector<T> __HH(NH*NH);
    MatrixMult_XY(NH,NL,NH,Y,A_LH,__HH,-1,0);
    SummMatrix_XD(NH,__HH,A_HH,1);
    MatrixMult_XY(NH,NH,NH,YX_1_inv,__HH,R_HH);
    MultMatrix_DX(NH,NH,D_inv_mat,R_HH,R_HH);
    __HH.resize(0);
    print("getting R_LH");
    for(size_t i=0;i<NH;++i){
        for(size_t j=0;j<NL;++j){
            A_LH[i*NL+j] -= X[i*NL+j]*A_HH[i];
        }
    }
    MatrixMult_XY(NL,NL,NH,XY_1_inv,A_LH,R_LH);
    MultMatrix_DX(NL,NH,A_inv_mat,R_LH,R_LH);
    A_HH.clear();
    A_LL.clear();
    A_LH.clear();
    A_HL.clear();
    X.clear();
    Y.clear();
    XY_1_inv.clear();
    YX_1_inv.clear();
    /**************************************/
    //print("R matrix = ");
    //printPairMatrix(NL,NH,R_LL.data(),R_LH.data(),R_HL.data(),R_HH.data());

    std::vector<T> D_L0_tmp = C_L;
    std::vector<T> D_L1_tmp(NL);
    std::vector<T> D_H0_tmp = C_H;
    std::vector<T> D_H1_tmp(NH);
    std::vector<T> D_L_acc(NL,0);
    std::vector<T> D_H_acc(NH,0);

    PowMatrix(NL,NH,R_LL,R_LH,R_HL,R_HH,skip_pow);
    print("Run scheme");
    return RunScheme(
                [&](auto & C_L,auto & C_H,auto & C_L_new,auto & C_H_new){
                gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,R_LL.data(),ld_A_LL,C_L.data(),1,0.0,C_L_new.data(),1);
                gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,R_LH.data(),ld_A_LH,C_H.data(),1,1.0,C_L_new.data(),1);
                gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,R_HL.data(),ld_A_HL,C_L.data(),1,0.0,C_H_new.data(),1);
                gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,R_HH.data(),ld_A_HH,C_H.data(),1,1.0,C_H_new.data(),1);
                },
        [](auto & D1,auto & D2,T tau,auto & Result){
            for(size_t i=0;i<D1.size();++i){
                Result[i] += 0.5*tau*(D1[i]+D2[i]);
            }
    },[&](auto &,auto &){
        T ret = 0;
        return std::make_tuple(ret,ret,ret,ret);
    },D_L0_tmp,D_H0_tmp,D_L1_tmp,D_H1_tmp,tau,T_full,view_moments);
    /*
    progress_bar PB(N_steps,100);
    size_t view_cout = 0;
    for(size_t i=1;i<=N_steps;++i){
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,R_LL.data(),ld_A_LL,D_L0_tmp.data(),1,0.0,D_L1_tmp.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,R_LH.data(),ld_A_LH,D_H0_tmp.data(),1,1.0,D_L1_tmp.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,R_HL.data(),ld_A_HL,D_L0_tmp.data(),1,0.0,D_H1_tmp.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,R_HH.data(),ld_A_HH,D_H0_tmp.data(),1,1.0,D_H1_tmp.data(),1);

        D_L_acc += 0.5*(D_L0_tmp+D_L1_tmp);
        D_H_acc += 0.5*(D_H0_tmp+D_H1_tmp);

        std::swap(D_L0_tmp,D_L1_tmp);
        std::swap(D_H0_tmp,D_H1_tmp);

        T N_L_tmp = vector_sum(D_L_acc);
        T N_H_tmp = vector_sum(D_H_acc);

        PB++;

        N_L[i] = N_L_tmp;
        N_H[i] = N_H_tmp;

        if(view_cout && i*tau <= view_moments[view_cout] &&
                view_moments[view_cout] < tau*(i+1)){
            N_L_view.push_back(D_L_acc/N_L_tmp);
            N_H_view.push_back(D_H_acc/N_H_tmp);
            view_cout ++;
        }

    }
    D_L_acc /= N_L.back();
    D_H_acc /= N_H.back();
    D_L1_tmp *=0;
    D_L1_tmp += D_L_acc;
    gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,-1.0,R_LL.data(),ld_A_LL,D_L_acc.data(),1,1.0,D_L1_tmp.data(),1);
    T S_LL = vector_sum(D_L1_tmp)/tau;
    gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,R_LH.data(),ld_A_LH,D_H_acc.data(),1,0.0,D_L1_tmp.data(),1);
    T S_LH = vector_sum(D_L1_tmp)/tau;
    gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,R_HL.data(),ld_A_HL,D_L_acc.data(),1,0.0,D_H1_tmp.data(),1);
    T S_HL = vector_sum(D_H1_tmp)/tau;

    D_H1_tmp *=0;
    D_H1_tmp += D_H_acc;
    gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,-1.0,R_HH.data(),ld_A_HH,D_H_acc.data(),1,1.0,D_H1_tmp.data(),1);
    T S_HH = vector_sum(D_H1_tmp)/tau;

    D_L0_tmp/=vector_sum(D_L0_tmp);
    D_H0_tmp/=vector_sum(D_H0_tmp);
    return ResultLH<T>{
        decltype(Result<T>::N_L)(moments,std::move(N_L)),
        decltype(Result<T>::N_H)(std::move(moments),std::move(N_H)),
        std::move(D_L0_tmp),
        std::move(D_H0_tmp),
        std::move(D_L_acc),
        std::move(D_H_acc),
        std::move(N_L_view),
        std::move(N_H_view),
        S_LL,S_LH,S_HL,S_HH
    };
    */
}

#define pget get<std::string>
#define LOGIF(condition,expression) if(condition) PVAR(expression)
template <typename T>
void make_work(const boost::property_tree::ptree& CP){
    MatrixFormat MF = MatrixFormat::BINARY;
    if(CP.find("format") != CP.not_found()){
        auto fmt =CP.pget("format");
        if(fmt == "bin" || fmt == "binary" || fmt == "b")
            MF = MatrixFormat::BINARY;
        else if(fmt == "text" || fmt == "txt" || fmt == "dat"){
            MF = MatrixFormat::TEXT;
        }
    }
    bool log_std = ptree_condition(CP,"debug",false);
    //LOGIF(log_std,MF);
    //bool isVectorMethod = true;//CP.get<std::string>("method","vector") != "matrix";

    auto FPath = [&](const std::string & tree_path){
        return config_path_from(CP.pgets(tree_path),CP.pgets
                                ("config_path")).string();
    };

    auto Mat0_LH = loadMatrix<T>(FPath("lh_in"),MF);
    auto Mat0_HL = loadMatrix<T>(FPath("hl_in"),MF);

    size_t NL = std::get<1>(Mat0_HL);
    size_t NH = std::get<2>(Mat0_HL);

    std::vector<T> Ev_H;
    if(CP.find("eh_in") != CP.not_found()){
        Ev_H = std::get<0>(loadMatrix<T>(FPath("eh_in"),MF));
    }

    std::vector<T> Ev_L;
    if(CP.find("el_in")  != CP.not_found()){
        Ev_L = std::get<0>(loadMatrix<T>(FPath("el_in"),MF));
    }

    std::vector<T> H_Distrib(NH,1.0/NH); //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP.find("dh_in") != CP.not_found()){
        H_Distrib = std::get<0>(loadMatrix<T>(FPath("dh_in"),MF))*
                CP.get<T>("h_nuf",1);
    }

    //std::vector<T> EH_D_out(EH_Distrib.size());

    std::vector<T> L_Distrib(NL,1.0/NL);
    if(CP.find("dl_in") != CP.not_found()){
        L_Distrib = std::get<0>(loadMatrix<T>(FPath("dl_in"),MF))*
                CP.get<T>("l_nuf",1);;
    }

    //std::vector<T> EL_D_out(EL_Distrib.size());

    //auto EL_Distrib = EH_Distrib*0;



    auto & A_HL_T = std::get<0>(Mat0_HL);
    size_t ld_A_HL = NH;

    auto & A_LH_T = std::get<0>(Mat0_LH);
    size_t ld_A_LH = NL;

    std::vector<T> A_LL(NL,0);
    size_t ld_A_LL = NL;
    std::vector<T> A_HH(NH,0);
    size_t ld_A_HH = NH;

    for(size_t i=0;i<NL;++i){
        A_LL[i] = -std::accumulate(&A_HL_T[i*NH],&A_HL_T[i*NH]+NH,0.0);
        if(Ev_L.size())
            A_LL[i] -= Ev_L[i];
    }
    for(size_t i=0;i<NH;++i){
        A_HH[i] = -std::accumulate(&A_LH_T[i*NL],&A_LH_T[i*NL]+NL,0.0);
        if(Ev_H.size())
            A_HH[i] -= Ev_H[i];
    }

    T tau = CP.get<float>("tau",-1);
    int degree = CP.get<int>("deg",-1);
    T fullT = CP.get<float>("T",-1);


    T approp_tau = 1/std::max(-min(A_HH),-min(A_LL));

    if(tau < 0){
        if(degree < 0){
            tau = 0.1*approp_tau;
            if(fullT < 0)
                degree = 1;
            else
                degree = fullT/tau + 0.5;
        }
        else{
            if(fullT < 0){
                tau = 0.1*approp_tau;
                fullT = tau*degree;
            }
            else{
                degree = fullT/tau + 0.5;
            }
        }
    }else{
        if(degree > 0){
           fullT = tau*degree;
        }
        else{
            if(fullT < 0){
                degree = 1;
                fullT = tau*degree;
            }
            else{
                degree = fullT/tau + 0.5;
            }
        }
    }
    size_t sp = 0;
    if(CP.get<std::string>("skip_pow","") == "auto"){
        sp = std::max(0,int(log2(tau/approp_tau)+2.5));
    }else{
        sp = CP.get<int>("skip_pow",0);
    }

    if(tau/(1<<sp) > approp_tau){
        print("WARNING: tau too big");
        print("tau_max = ",approp_tau);
    }
    if(log_std){
        PVAR(tau);
        PVAR(tau/(1<<sp));
        PVAR(degree);
        PVAR(fullT);
    }

    /*
    A_LL*=tau;
    A_HH*=tau;
    A_HL_T*=tau;
    A_LH_T*=tau;
    Ev_L *=tau;
    Ev_H *=tau;

    for(size_t i=0;i<NL;++i){
        A_LL[i*NL+i] += 1;
    }
    for(size_t i=0;i<NH;++i){
        A_HH[i*NH+i] += 1;
    }
    */

    if(log_std and MF == MatrixFormat::TEXT){
        printPairMatrixD(NL,NH,A_LL.data(),A_LH_T.data(),A_HL_T.data(),A_HH.data());
        print();
    }


    //size_t iter_num = isVectorMethod ? degree : 1.5+log2(degree);
    /*
    auto RD = isVectorMethod ? VectorMultMethod<T>(NH,NL,A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,
                                                iter_num,tau) :
                               MatrixMultMethod<T>(NH,NL,A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,iter_num,tau);
    */

    Result<T> RD = CP.pgets("method","euler") == "euler" ?
                SchemeEulerMult<T>(A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,tau,fullT,
                                   sp):
                Scheme2Order<T>(A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,tau,fullT,
                                sp);


    if(log_std and MF == MatrixFormat::TEXT){
        print("D_vecs");
        printPairVector(NL,NH,RD.D_L_final.data(),RD.D_H_final.data());
    }
    //T S_LH_agg,S_HL_agg, E_H_agg, E_L_agg;

    std::ofstream pit_params(config_path_from(CP.pgets("pout","mat_pout.json.txt"),CP.pgets("config_path")));
    boost::property_tree::ptree P;
    P.put("ll_av",RD.S_LL);
    P.put("lh_av",RD.S_LH);
    P.put("hl_av",RD.S_HL);
    P.put("hh_av",RD.S_HH);

    boost::property_tree::write_json(pit_params,P);
    if(log_std and MF == MatrixFormat::TEXT){
        print("D_distrib_vecs");
        printPairVector(NL,NH,RD.D_L_final.data(),RD.D_H_final.data());
    }

    saveMatrix(RD.C_L_final.data(),1,NL,FPath("CL_out"),MF);
    saveMatrix(RD.C_H_final.data(),1,NH,FPath("CH_out"),MF);
    saveMatrix(RD.D_L_final.data(),1,NL,FPath("DL_out"),MF);
    saveMatrix(RD.D_H_final.data(),1,NH,FPath("DH_out"),MF);

    if(CP.find("nl_out") != CP.not_found()){
        std::ofstream nl_f(FPath("nl_out"));
        grob::as_csv(RD.n_l).save(nl_f);
    }
    if(CP.find("nh_out") != CP.not_found()){
        std::ofstream nh_f(FPath("nh_out"));
        grob::as_csv(RD.n_h).save(nh_f);
    }
    if(CP.find("NL_out") != CP.not_found()){
        std::ofstream NL_f(FPath("NL_out"));
        grob::as_csv(RD.N_L).save(NL_f);
    }
    if(CP.find("NH_out")  != CP.not_found()){
        std::ofstream NH_f(FPath("NH_out"));
        grob::as_csv(RD.N_H).save(NH_f);
    }
}

int main(int argc, char ** argv){
    boost::property_tree::ptree CP;
    auto ret_key = parse_command_line(argc,argv,CP);
    if(!ret_key.empty()){
        std::cout << "error : " << ret_key<<std::endl;
        return 0;
    }

    std::map<std::string,std::string> default_params{{"dl_out","dl_out.mat"},{"dh_out","dl_out.mat"},
                                                    {"DL_out","DL_out.mat"},{"DH_out","DH_out.mat"}};
    std::vector<std::string> required_params{"lh_in","hl_in"};
    for(const auto & param : required_params){
        if(CP.find(param)==CP.not_found()){
            std::cout << "required param " << param << std::endl;
            return 0;
        }
    }
    for(const auto & [param,value] : default_params){
        if(CP.find(param)==CP.not_found()){
            CP.put(param,value);
        }
    }
//    for(size_t i=0;i<CP.size();++i){
//        std::cout << key_names[i] << " : " << CP[i] << std::endl;
//    }
    if(CP.pgets("debug","false") == "true"){
        boost::property_tree::write_json(std::cout,CP);
    }
    if(CP.find("type")!=CP.not_found() && CP.pget("type") == "float")
        make_work<float>(CP);
    else
        make_work<double>(CP);

    return 0;
}
