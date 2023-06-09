#ifndef CBALS_FUNCTIONS_HPP
#define CBALS_FUNCTIONS_HPP
#include <cblas.h>
#include <lapacke.h>
#include <string>
#include <iostream>
#include <debug/debugdef.hpp>
#include "complex/complex_ex.hpp"
#include <grid_objects.hpp>
template <typename T>
constexpr auto gemm = std::conditional<std::is_same<T,float>::value,
                                    std::integral_constant<decltype((cblas_sgemm)),cblas_sgemm>,
                                    std::integral_constant<decltype((cblas_dgemm)),cblas_dgemm>>::type::value;
template <typename T>
constexpr auto gemv = std::conditional<std::is_same<T,float>::value,
                                    std::integral_constant<decltype((cblas_sgemv)),cblas_sgemv>,
                                    std::integral_constant<decltype((cblas_dgemv)),cblas_dgemv>>::type::value;

template <typename T>
constexpr auto gbmv = std::conditional<std::is_same<T,float>::value,
                                    std::integral_constant<decltype((cblas_sgbmv)),cblas_sgbmv>,
                                    std::integral_constant<decltype((cblas_dgbmv)),cblas_dgbmv>>::type::value;

template <typename T>
constexpr auto getri =  std::conditional<std::is_same<T,float>::value,
                        std::integral_constant<decltype((LAPACKE_sgetri)),LAPACKE_sgetri>,
                        std::integral_constant<decltype((LAPACKE_dgetri)),LAPACKE_dgetri>>::type::value;
template <typename T>
constexpr auto getrf =  std::conditional<std::is_same<T,float>::value,
                        std::integral_constant<decltype((LAPACKE_sgetrf)),LAPACKE_sgetrf>,
                        std::integral_constant<decltype((LAPACKE_dgetrf)),LAPACKE_dgetrf>>::type::value;




template <typename T>
void printMatrix(size_t N,size_t M,const T * data){
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<M;++j){
            std::cout << data[j*N+i] << "\t";
        }std::cout << std::endl;
    }std::cout << std::endl;
}
template <typename T>
void  printPairMatrix(size_t NL,size_t NH,const T * LL,const T*LH,const  T*HL,const T*HH){
    for(size_t il=0;il<NL;++il){
        for(size_t jl=0;jl<NL;++jl){
            std::cout << LL[il+jl*NL] << "\t";
        }
        for(size_t jh=0;jh<NH;++jh){
            std::cout << LH[il+jh*NL] << "\t";
        }std::cout << std::endl;
    }
    for(size_t ih=0;ih<NH;++ih){
        for(size_t jl=0;jl<NL;++jl){
            std::cout << HL[ih+jl*NH] << "\t";
        }
        for(size_t jh=0;jh<NH;++jh){
            std::cout << HH[ih+jh*NH] << "\t";
        }std::cout << std::endl;
    }
}
template <typename T>
void  printPairMatrixD(size_t NL,size_t NH,const T * LL,const T*LH,const  T*HL,const T*HH){
    for(size_t il=0;il<NL;++il){
        for(size_t jl=0;jl<NL;++jl){
            std::cout << (il == jl ? LL[il] : 0)<< "\t";
        }
        for(size_t jh=0;jh<NH;++jh){
            std::cout << LH[il+jh*NL] << "\t";
        }std::cout << std::endl;
    }
    for(size_t ih=0;ih<NH;++ih){
        for(size_t jl=0;jl<NL;++jl){
            std::cout << HL[ih+jl*NH] << "\t";
        }
        for(size_t jh=0;jh<NH;++jh){
            std::cout <<(ih == jh ?  HH[ih] :0) << "\t";
        }std::cout << std::endl;
    }
}

template <typename T>
void  printPairVector(size_t NL,size_t NH,const T * L,const T*H){
    for(size_t il=0;il<NL;++il){
        std::cout << L[il] << "\t";
    }
    for(size_t ih=0;ih<NH;++ih){
        std::cout << H[ih] << "\t";
    }
    print();
}

template <typename T>
void PowMatrix(size_t NL,size_t NH,std::vector<T> &A_LL,std::vector<T> & A_LH,std::vector<T> & A_HL,std::vector<T> & A_HH,
               size_t log_2_pow){

    if(log_2_pow == 0)
        return;
    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;

    std::vector<T> B_LL(A_LL.size());
    std::vector<T> B_LH(A_LH.size());
    std::vector<T> B_HL(A_HL.size());
    std::vector<T> B_HH(A_HH.size());
    print("exponentiat into pow 2^" + std::to_string(log_2_pow) +
          " = "+ std::to_string(1<<log_2_pow));
    progress_bar PB(log_2_pow,100);
    while(log_2_pow){
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NL,
                    1.0,A_LL.data(),ld_A_LL,A_LL.data(),ld_A_LL,0.0,B_LL.data(),ld_A_LL);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NH,
                            1.0,A_LH.data(),ld_A_LH,A_HL.data(),ld_A_HL,1.0,B_LL.data(),ld_A_LL);
        //DGEMM()

        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NL,
                    1.0,A_LL.data(),ld_A_LL,A_LH.data(),ld_A_LH,0.0,B_LH.data(),ld_A_LH);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NH,
                            1.0,A_LH.data(),ld_A_LH,A_HH.data(),ld_A_HH,1.0,B_LH.data(),ld_A_LH);


        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NL,
                    1.0,A_HL.data(),ld_A_HL,A_LL.data(),ld_A_LL,0.0,B_HL.data(),ld_A_HL);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NH,
                            1.0,A_HH.data(),ld_A_HH,A_HL.data(),ld_A_HL,1.0,B_HL.data(),ld_A_HL);

        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NL,
                    1.0,A_HL.data(),ld_A_HL,A_LH.data(),ld_A_LH,0.0,B_HH.data(),ld_A_HH);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NH,
                            1.0,A_HH.data(),ld_A_HH,A_HH.data(),ld_A_HH,1.0,B_HH.data(),ld_A_HH);


        std::swap(A_LL,B_LL);
        std::swap(A_HL,B_HL);
        std::swap(A_LH,B_LH);
        std::swap(A_HH,B_HH);
        PB++;
        log_2_pow--;
    }
}



template  <typename T>
struct BlockMatrix{
    std::vector<T> LL;
    std::vector<T> LH;
    std::vector<T> HL;
    std::vector<T> HH;
};
template <typename T>
void MultRMatrix(size_t NL,size_t NH,
                 std::vector<T> const &A_LL,std::vector<T> const & A_LH,std::vector<T> const & A_HL,std::vector<T> const & A_HH,
                 std::vector<T> const &B_LL,std::vector<T> const & B_LH,std::vector<T> const & B_HL,std::vector<T> const & B_HH,
                 std::vector<T> &C_LL,std::vector<T> & C_LH,std::vector<T> & C_HL,std::vector<T> & C_HH){
    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;
    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NL,
                1.0,A_LL.data(),ld_A_LL,B_LL.data(),ld_A_LL,0.0,C_LL.data(),ld_A_LL);
    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NH,
                        1.0,A_LH.data(),ld_A_LH,B_HL.data(),ld_A_HL,1.0,C_LL.data(),ld_A_LL);
    //DGEMM()

    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NL,
                1.0,A_LL.data(),ld_A_LL,B_LH.data(),ld_A_LH,0.0,C_LH.data(),ld_A_LH);
    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NH,
                        1.0,A_LH.data(),ld_A_LH,B_HH.data(),ld_A_HH,1.0,C_LH.data(),ld_A_LH);


    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NL,
                1.0,A_HL.data(),ld_A_HL,B_LL.data(),ld_A_LL,0.0,C_HL.data(),ld_A_HL);
    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NH,
                        1.0,A_HH.data(),ld_A_HH,B_HL.data(),ld_A_HL,1.0,C_HL.data(),ld_A_HL);

    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NL,
                1.0,A_HL.data(),ld_A_HL,B_LH.data(),ld_A_LH,0.0,C_HH.data(),ld_A_HH);
    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NH,
                        1.0,A_HH.data(),ld_A_HH,B_HH.data(),ld_A_HH,1.0,C_HH.data(),ld_A_HH);
}
template <typename T>
void MultRMatrix(size_t N,
                 std::vector<T> const &A,
                 std::vector<T> const &B,
                 std::vector<T> &C){
    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,
                1.0,A.data(),N,B.data(),N,0.0,C.data(),N);

}
template <typename T>
BlockMatrix<T> PowMatrixSumm(size_t NL,size_t NH,std::vector<T> &A_LL,std::vector<T> & A_LH,std::vector<T> & A_HL,std::vector<T> & A_HH,
               size_t log_2_pow,T tau){

    BlockMatrix<T> Summ{
        A_LL*(tau/2),
        A_LH*(tau/2),
        A_HL*(tau/2),
        A_HH*(tau/2)
    };
    for(size_t i=0;i<NL;++i){
        Summ.LL[NL*i+i] += tau/2;
    }
    for(size_t j=0;j<NH;++j){
        Summ.HH[NH*j+j] += tau/2;
    }
    if(log_2_pow == 0)
        return Summ;

    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;

    std::vector<T> B_LL(A_LL.size());
    std::vector<T> B_LH(A_LH.size());
    std::vector<T> B_HL(A_HL.size());
    std::vector<T> B_HH(A_HH.size());

    BlockMatrix<T> Summ_tmp{std::vector<T>(A_LL.size()),
                            std::vector<T>(A_LH.size()),
                            std::vector<T>(A_HL.size()),
                            std::vector<T>(A_HH.size())};

    print("exponentiat into pow 2^" + std::to_string(log_2_pow) +
          " = "+ std::to_string(1<<log_2_pow));
    progress_bar PB(log_2_pow,100);
    while(log_2_pow){

        MultRMatrix(NL,NH,A_LL,A_LH,A_HL,A_HH,
                    Summ.LL,Summ.LH,Summ.HL,Summ.HH,
                    Summ_tmp.LL,Summ_tmp.LH,Summ_tmp.HL,Summ_tmp.HH);
        Summ.LL += Summ_tmp.LL;
        Summ.LH += Summ_tmp.LH;
        Summ.HL += Summ_tmp.HL;
        Summ.HH += Summ_tmp.HH;

        MultRMatrix(NL,NH,A_LL,A_LH,A_HL,A_HH,
                    A_LL,A_LH,A_HL,A_HH,
                    B_LL,B_LH,B_HL,B_HH);

        std::swap(A_LL,B_LL);
        std::swap(A_HL,B_HL);
        std::swap(A_LH,B_LH);
        std::swap(A_HH,B_HH);
        PB++;
        log_2_pow--;
    }
    return Summ;
}

template <typename T>
void PowMatrix(size_t N,std::vector<T> &A,
               size_t log_2_pow){


    if(log_2_pow == 0)
        return;


    std::vector<T> B(A.size());



    print("exponentiat into pow 2^" + std::to_string(log_2_pow) +
          " = "+ std::to_string(1<<log_2_pow));
    progress_bar PB(log_2_pow,100);
    while(log_2_pow){


        MultRMatrix(N,A,
                    A,
                    B);

        std::swap(A,B);
        PB++;
        log_2_pow--;
    }

}

template <typename T>
std::vector<T> PowMatrixSumm(size_t N,std::vector<T> &A,
               size_t log_2_pow,T tau){

    std::vector<T> Summ{
        A*(tau/2)
    };
    for(size_t i=0;i<N;++i){
        Summ[N*i+i] += tau/2;
    }
    if(log_2_pow == 0)
        return Summ;


    std::vector<T> B(A.size());


    std::vector<T> Summ_tmp(N*N);

    print("exponentiat into pow 2^" + std::to_string(log_2_pow) +
          " = "+ std::to_string(1<<log_2_pow));
    progress_bar PB(log_2_pow,100);
    while(log_2_pow){

        MultRMatrix(N,A,
                    Summ,
                    Summ_tmp);

        Summ += Summ_tmp;

        MultRMatrix(N,A,
                    A,
                    B);

        std::swap(A,B);
        PB++;
        log_2_pow--;
    }
    return Summ;
}

template <typename T>
std::vector<T> RMatrixEuler(size_t N,
                            std::vector<T> &A,
                            T tau){

    std::vector<T> R = A*tau;

    for(size_t i=0;i<N;++i){
        R[N*i+i] = 1;
    }
    return R;
}

template <typename T>
std::vector<T> RMatrix2Order(size_t N,
                            std::vector<T> &A,
                            T tau){

    std::vector<T> R_1 = A*(-tau/2);

    for(size_t i=0;i<N;++i){
        R_1[N*i+i] += 1;
    }
    std::vector<int32_t> ipiv_L(N+1);
    auto ret_L = getrf<T>(CblasColMajor,N,N,R_1.data(),N,ipiv_L.data());
    if(ret_L!=0){
        throw std::runtime_error("error decompoeing L*U matrix, ret_L = " + std::to_string(ret_L));
    }
    getri<T>(CblasColMajor,N,R_1.data(),N,ipiv_L.data());
    if(ret_L!=0){
        throw std::runtime_error("error inverting matrix, ret_L = " + std::to_string(ret_L));
    }

    std::vector<T> R_2 = R_1;

    gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,
                tau/2.0,R_1,N,A.data(),N,0.0,R_2.data(),N);
    return R_2;
}
template <typename T>
BlockMatrix<T> RMatrixEuler(size_t NL,size_t NH,
                            std::vector<T> &A_LL,std::vector<T> & A_LH,
                            std::vector<T> & A_HL,std::vector<T> & A_HH,
                            T tau){

    std::vector<T> R_LL(NL*NL,0);
    std::vector<T> R_LH = A_LH*tau;
    A_LH.clear();
    std::vector<T> R_HL = A_HL*tau;
    A_HL.clear();
    std::vector<T> R_HH(NH*NH,0);
    for(size_t i=0;i<NL;++i){
        R_LL[NL*i+i] = 1 + tau*A_LL[i];
    }
    for(size_t i=0;i<NH;++i){
        R_HH[NH*i+i] = 1 + tau*A_HH[i];
    }
    return BlockMatrix<T>{std::move(R_LL),std::move(R_LH),
                        std::move(R_HL),std::move(R_HH)};

}

template <typename T>
BlockMatrix<T> RMatrix2Order(size_t NL,size_t NH,
                            std::vector<T> &A_LL,std::vector<T> & A_LH,
                            std::vector<T> & A_HL,std::vector<T> & A_HH,
                            T tau){


    size_t ld_A_HL = NH;
    size_t ld_A_LH = NL;
    size_t ld_A_HH = NH;
    size_t ld_A_LL = NL;


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

    A_LL *= (0.5*tau);
    A_LH *= (0.5*tau);
    A_HL *= (0.5*tau);
    A_HH *= (0.5*tau);

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

    return BlockMatrix<T>{std::move(R_LL),std::move(R_LH),
                          std::move(R_HL),std::move(R_HH)};
}

template <typename T>
struct ResultLH_F{
    grob::GridUniform<T> T_grid;
    std::vector<T> c_l_t,c_h_t;
    std::vector<T> n_l_t,n_h_t;
    std::vector<T>  a_t;
    std::vector<T> C_L_F,C_H_F;
    std::vector<T> D_L_F,D_H_F;
};
template <typename T>
struct Result_F{
    grob::GridUniform<T> T_grid;
    std::vector<T> c_t;
    std::vector<T> n_t;
    std::vector<T>  a_t;
    std::vector<T> C_F;
    std::vector<T> D_F;
};



template <typename ResultContainer_t,
          typename ElementFunction_t,
          typename...Containers_t>
void vector_apply(ResultContainer_t & Result,
                  ElementFunction_t && Function,
                  Containers_t const &...Containers){
    for(size_t i=0;i<Result.size();++i){
        Result[i] = Function(Containers[i]...);
    }
}

enum class AnnCalcPolicy{
    IGNORE,
    ONLY_OUT,
    FULL
};

template <typename T>
ResultLH_F<T> Evolution_LH(size_t NL,size_t NH,std::vector<T> &A_LL,std::vector<T> & A_LH,
                        std::vector<T> & A_HL,std::vector<T> & A_HH,
                           std::string scheme,
                           std::vector<T> & ANN_HL,
                           AnnCalcPolicy AP,
                           std::vector<T> &C_L,std::vector<T> & C_H,
                           std::vector<T> &C_L_I,std::vector<T> & C_H_I,
                           std::vector<T> &D_L_I,std::vector<T> & D_H_I,
                        T tau,T T_full,size_t skip_pow =0){
    T tau_eff = tau / (1 << skip_pow);
    BlockMatrix<T> R = (scheme == "order2" ? RMatrix2Order(NL,NH,A_LL,A_LH,A_HL,A_HH,tau_eff) :
                                   RMatrixEuler(NL,NH,A_LL,A_LH,A_HL,A_HH,tau_eff));
    BlockMatrix<T> SM = PowMatrixSumm(NL,NH,R.LL,R.LH,R.HL,R.HH,skip_pow,tau_eff);
    auto gemv_LH = [NL,NH](std::vector<T> & _LL,std::vector<T> & _LH,
                        std::vector<T> & _HL, std::vector<T> & _HH,
                        std::vector<T> & _VL, std::vector<T> & _VH,
                        std::vector<T> & _V1L, std::vector<T> & _V1H,
                        T alpha = 1,T betha = 0){
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,alpha,_LL.data(),NL,_VL.data(),1,betha,_V1L.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,alpha,_LH.data(),NL,_VH.data(),1,1.0,_V1L.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,alpha,_HL.data(),NH,_VL.data(),1,betha,_V1H.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,alpha,_HH.data(),NH,_VH.data(),1,1.0,_V1H.data(),1);
    };


    auto get_ann_matrix = [NL,NH,&ANN_HL](std::vector<T> & _D_L, std::vector<T> & _D_H,
            std::vector<T> & _ANN_L, std::vector<T> & _ANN_H){
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,ANN_HL.data(),NH,_D_L.data(),1,0.0,_ANN_H.data(),1);
        gemv<T>(CblasColMajor,CblasTrans,NL,NH,1.0,ANN_HL.data(),NH,_D_H.data(),1,0.0,_ANN_L.data(),1);
    };

    auto v_exp = [](std::vector<T> & _A,T _t,std::vector<T> & Result){
        if(Result.size() != _A.size()){
            Result.resize(_A.size());
        }
        for(size_t i=0;i<_A.size();++i){
            Result[i] = exp(_A[i]*_t);
        }
    };

    auto diag_v_mult = [](std::vector<T> const  &X,std::vector<T>const   & Y,
                        std::vector<T> & Z){
        for(size_t i=0;i<X.size();++i){
            Z[i] = X[i]*Y[i];
        }
    };

    std::vector<T> D_L = (D_L_I.size() ? D_L_I : std::vector<T>(NL,0));
    std::vector<T> D_H =   (D_H_I.size() ? D_H_I : std::vector<T>(NH,0));
    std::vector<T> C_L_tmp_1 = (C_L_I.size() ? C_L_I : C_L);
    std::vector<T> C_H_tmp_1 = (C_H_I.size() ? C_H_I : C_H);

    std::vector<T> C_L_tmp_2 = C_L_tmp_1;
    std::vector<T> C_H_tmp_2 = C_H_tmp_1;

    std::vector<T> C_L_eff(NL);
    std::vector<T> C_H_eff(NH);


    size_t N_steps = T_full/tau+0.5;


    bool calc_ann = ANN_HL.size()>0;
    std::vector<T> ANN_L(calc_ann ? NL : 0);
    std::vector<T> ANN_H(calc_ann ? NH : 0);

    std::vector<T> R_ANN_L(calc_ann ? NL : 0);
    std::vector<T> R_ANN_H(calc_ann ? NH : 0);

    std::vector<T> n_l_t(N_steps+1);
    n_l_t[0] = vector_sum(D_L);
    std::vector<T> n_h_t(N_steps+1);
    n_h_t[0] = vector_sum(D_H);;

    std::vector<T> c_l_t(N_steps+1);
    c_l_t[0] = vector_sum(C_L_tmp_1);
    std::vector<T> c_h_t(N_steps+1);
    c_h_t[0] = vector_sum(C_H_tmp_1);

    std::vector<T> ann_t(N_steps+1,0);
    get_ann_matrix(D_L,D_H,ANN_L,ANN_H);
    ann_t[0] = vector_sum(ANN_L*D_L);

    progress_bar PB(N_steps,100);

    for(size_t i=1;i<=N_steps;++i){
        if(calc_ann  && (AP != AnnCalcPolicy::IGNORE) ){
            get_ann_matrix(D_L,D_H,ANN_L,ANN_H);
            ann_t[i] = 0.5*vector_sum(ANN_L*D_L);
            if(AP == AnnCalcPolicy::FULL){
                v_exp(ANN_L,-tau/2,R_ANN_L);
                v_exp(ANN_H,-tau/2,R_ANN_H);

                diag_v_mult(R_ANN_L,C_L_tmp_1,C_L_tmp_1);
                diag_v_mult(R_ANN_H,C_H_tmp_1,C_H_tmp_1);
            }
        }
        gemv_LH(SM.LL,SM.LH,SM.HL,SM.HH,
                C_L_tmp_1,C_H_tmp_1,C_L_eff,C_H_eff);
        D_L += C_L_eff;
        D_H += C_H_eff;
        gemv_LH(R.LL,R.LH,R.HL,R.HH,
                C_L_tmp_1,C_H_tmp_1,
                C_L_tmp_2,C_H_tmp_2);
        std::swap(C_L_tmp_1,C_L_tmp_2);
        std::swap(C_H_tmp_1,C_H_tmp_2);
        if(calc_ann && (AP == AnnCalcPolicy::FULL)){
            get_ann_matrix(D_L,D_H,ANN_L,ANN_H);
            ann_t[i] += 0.5*vector_sum(ANN_L*D_L);
            v_exp(ANN_L,-tau/2,R_ANN_L);
            v_exp(ANN_H,-tau/2,R_ANN_H);

            diag_v_mult(R_ANN_L,C_L_tmp_1,C_L_tmp_1);
            diag_v_mult(R_ANN_H,C_H_tmp_1,C_H_tmp_1);
        }

        n_l_t[i] = vector_sum(D_L);
        n_h_t[i] = vector_sum(D_H);;
        c_l_t[i] = vector_sum(C_L_tmp_1);
        c_h_t[i] = vector_sum(C_H_tmp_1);
        PB++;
    }
    return ResultLH_F<T>{
        Vector(N_steps+1,[=](size_t i){return i*tau;}),
        std::move(c_l_t),std::move(c_h_t),
        std::move(n_l_t),std::move(n_h_t),
        std::move(ann_t),
        std::move(C_L_tmp_1),std::move(C_H_tmp_1),
        std::move(D_L),std::move(D_H)
    };
}



template <typename T>
Result_F<T> Evolution(size_t N,size_t NH,std::vector<T> &A,
                           std::string scheme,
                           std::vector<T> & ANN,
                           AnnCalcPolicy AP,
                           std::vector<T> &C,
                           std::vector<T> &C_I,
                           std::vector<T> &D_I,
                        T tau,T T_full,size_t skip_pow =0){

    T tau_eff = tau / (1 << skip_pow);
    std::vector<T> R = (scheme == "order2" ? RMatrix2Order(N,A,tau_eff) :
                                   RMatrixEuler(N,A,tau_eff));

    std::vector<T> SM = PowMatrixSumm(N,R,skip_pow,tau_eff);
    auto gemv_ = [N](std::vector<T> & _A,
                        std::vector<T> & _V,
                        std::vector<T> & _V1,
                        T alpha = 1,T betha = 0){
        gemv<T>(CblasColMajor,CblasNoTrans,N,N,alpha,_A.data(),N,_V.data(),1,betha,_V1.data(),1);
    };


    auto get_ann_matrix = [N,&ANN](std::vector<T> & _D,
            std::vector<T> & _ANN){
        gemv<T>(CblasColMajor,CblasNoTrans,N,N,1.0,ANN.data(),N,_D.data(),1,0.0,_ANN.data(),1);
    };

    auto v_exp = [](std::vector<T> & _A,T _t,std::vector<T> & Result){
        if(Result.size() != _A.size()){
            Result.resize(_A.size());
        }
        for(size_t i=0;i<_A.size();++i){
            Result[i] = exp(_A[i]*_t);
        }
    };

    auto diag_v_mult = [](std::vector<T> const  &X,std::vector<T>const   & Y,
                        std::vector<T> & Z){
        for(size_t i=0;i<X.size();++i){
            Z[i] = X[i]*Y[i];
        }
    };

    std::vector<T> D_ = (D_I.size() ? D_I : std::vector<T>(N,0));

    std::vector<T> C_tmp_1 = (C_I.size() ? C_I : C);

    std::vector<T> C_tmp_2 = C_tmp_1;

    std::vector<T> C_eff(N);

    gemv_LH(SM,
            C,C_eff);

    size_t N_steps = T_full/tau+0.5;


    bool calc_ann = ANN.size()>0;
    std::vector<T> ANN_(calc_ann ? N : 0);

    std::vector<T> R_ANN(calc_ann ? N : 0);

    std::vector<T> n_t(N_steps+1);
    n_t[0] = vector_sum(D_);

    std::vector<T> c_t(N_steps+1);
    c_t[0] = vector_sum(C_tmp_1);

    std::vector<T> ann_t(N_steps+1);
    get_ann_matrix(D_,ANN_);
    ann_t[0] = vector_sum(ANN_*D_);

    progress_bar PB(N_steps,100);

    for(size_t i=1;i<=N_steps;++i){
        if(calc_ann && (AP != AnnCalcPolicy::IGNORE)){
            get_ann_matrix(D_,ANN_);
            ann_t[i] = 0.5*vector_sum(ANN_*D_);
            if(AP == AnnCalcPolicy::FULL){
                v_exp(ANN_,-tau/2,R_ANN);
                v_exp(ANN_,-tau/2,R_ANN);
                //diag_v_mult(R_ANN,D_,D_);

                diag_v_mult(R_ANN,C_tmp_1,C_tmp_1);
            }
        }
        gemv_(SM,C_tmp_1,
                C_eff);
        D_ += C_eff;
        gemv_(R,C_tmp_1,
                C_tmp_2);
        std::swap(C_tmp_1,C_tmp_2);
        if(calc_ann && (AP == AnnCalcPolicy::FULL)){
            get_ann_matrix(D_,ANN_);
            ann_t[i] = 0.5*vector_sum(ANN_*D_);
            v_exp(ANN_,-tau/2,R_ANN);
            v_exp(ANN_,-tau/2,R_ANN);
            //diag_v_mult(R_ANN,D_,D_);
            diag_v_mult(R_ANN,C_tmp_1,C_tmp_1);
        }

        n_t[i] = vector_sum(D_);
        c_t[i] = vector_sum(C_tmp_1);
        PB++;
    }
    return Result_F<T>{
        Vector(N_steps+1,[=](size_t i){return i*tau;}),
        std::move(c_t),
        std::move(n_t),
        std::move(ann_t),
        std::move(C_tmp_1),
        std::move(D_)
    };
}


#endif//CBALS_FUNCTIONS_HPP
