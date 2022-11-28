#ifndef MATRIX_FUNCTIONS_H
#define MATRIX_FUNCTIONS_H
#include "../functions.hpp"


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






#endif
