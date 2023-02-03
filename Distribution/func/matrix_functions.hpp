#ifndef MATRIX_FUNCTIONS_H
#define MATRIX_FUNCTIONS_H
#include <utils>
#include <random>
#include <cmath>
#include <tuple>
#include <type_traits>


typedef std::vector<double> dvector;
typedef std::vector<std::vector<double>> dmatrix;

double vnorm_0(const dvector&V){

    double nv = 0;
    for(size_t i=0;i<V.size();++i){
        if(nv < std::abs(V[i]))
            nv = std::abs(V[i]);
    }
    return nv;
}

double vnorm_1(const dvector&V){
    double nv = 0;
    for(size_t i=0;i<V.size();++i){
        nv += std::abs(V[i]);
    }
    return nv;
}

double mnorm_0(const dmatrix&V){
    double nm = 0;
    for(size_t i=0;i<V.size();++i){
        nm = std::max(nm,vnorm_1(V[i]));
    }
    return nm;
}

double mnorm_1(const dmatrix&V){
    double nm = 0;
    for(size_t i=0;i<V.size();++i){
        nm += vnorm_0(V[i]);
    }
    return nm;
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




template <typename T >
void SaveMatrixString(const T *X,size_t N,size_t M,std::ostream &S){
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<M;++j)
            S << X[i*M+j] << "\t";
        S << std::endl;
    }
}
template <typename T >
void SaveMatrixString(const T *X,size_t N,size_t M,const std::string&fname){
    std::ofstream _of(fname);
    SaveMatrixString<T>(X,N,M,_of);
}


template <typename T>
auto LoadMatrixString(std::istream &S){

    size_t M = 0;
    size_t N = 1;
    std::vector<T> retMatrix;
    std::string first_line;
    std::getline(S,first_line);
    std::stringstream Line1(first_line);
    double x;
    while (Line1 >> x) {
        retMatrix.push_back(x);
        M++;
    }
    while(S>>x){
        retMatrix.push_back(x);
    }
    N = retMatrix.size()/M;
    return _T(retMatrix,N,M);
}
template <typename T>
auto LoadMatrixString(const std::string&fname){
    std::ifstream S;
    S.open(fname);
    return LoadMatrixString<T>(S);
}



template <typename T>
void SaveMatrixBinary(const T*X,size_t N,size_t M,std::ostream &S){
    S.write((char*)&N,sizeof(N));
    S.write((char*)&M,sizeof(M));
    S.write((char*)X,N*M*sizeof(double));
}

template <typename T>
void SaveMatrixBinary(const T*X,size_t N,size_t M,const std::string&fname){
    std::ofstream _of;
    _of.open(fname,std::ios::binary);
    SaveMatrixBinary<T>(X,N,M,_of);
}

template <typename T>
auto LoadMatrixBinary(std::istream &S){
    size_t N;
    size_t M;
    S.read((char*)&N,sizeof(N));
    S.read((char*)&M,sizeof(M));
    std::vector<T> X(N*M);
    S.read((char*)X.data(),N*M*sizeof(double));
    return _T(X,N,M);
}
template <typename T>
auto LoadMatrixBinary(const std::string&fname){
    std::ifstream _in;
    _in.open(fname,std::ios::binary);
    return LoadMatrixBinary<T>(_in);
}

bool is_binary(const std::string & fname){
    std::regex end(".*\\.(.*)");
    std::smatch sm;
    regex_match(fname, sm, end);
    if(sm.size() > 1){
        std::string fend = sm[1];
        return fend == "bmat" || fend == "bvec" || fend == "bmatrix" || fend == "bvector";
    }
    return 0;

}

enum class MatrixFormat{
    BINARY,TEXT
};
template <typename T>
void saveMatrix(const T*Mat,size_t N,size_t M,const std::string&fname,MatrixFormat MF = MatrixFormat::TEXT){
    if(MF == MatrixFormat::BINARY)
        SaveMatrixBinary<T>(Mat,N,M,fname);
    else
        SaveMatrixString<T>(Mat,N,M,fname);
}
template <typename T>
auto loadMatrix(const std::string&fname,MatrixFormat MF = MatrixFormat::TEXT){
    if(MF == MatrixFormat::BINARY)
        return LoadMatrixBinary<T>(fname);
    else
        return LoadMatrixString<T>(fname);
}

template <typename T>
void saveMatrixSmart(const T*Mat,size_t N,size_t M,const std::string&fname){
    if(is_binary(fname))
        SaveMatrixBinary<T>(Mat,N,M,fname);
    else
        SaveMatrixString<T>(Mat,N,M,fname);
}
template <typename T,typename string_type>
void SaveMatrix(const T*Mat,size_t N,size_t M, string_type&&fname,MatrixFormat MF){
    if(MF == MatrixFormat::BINARY)
        SaveMatrixBinary<T>(Mat,N,M,fname);
    else
        SaveMatrixString<T>(Mat,N,M,fname);
}

template <typename T>
auto loadMatrixSmart(const std::string&fname){
    if(is_binary(fname))
        return LoadMatrixBinary<T>(fname);
    else
        return LoadMatrixString<T>(fname);
}

template <typename T,typename string_type>
auto LoadMatrix(string_type&&fname,MatrixFormat MF){
    if(MF == MatrixFormat::BINARY)
        return LoadMatrixBinary<T>(fname);
    else
        return LoadMatrixString<T>(fname);
}





void PrintMatrix(const std::vector<std::vector<double>> &Mat,std::ostream &os = std::cout){
    for(size_t i=0;i<Mat.size();++i){
        for(size_t j=0;j<Mat[i].size();++j){
            os << Mat[i][j] <<"\t";
        }
        os<<std::endl;
    }
}


auto standart_mult = [](const auto &x,const auto &y){return x*y;};
template <typename T = double,typename MultFunctType = decltype(standart_mult)>
auto smart_pow(const T &x,size_t N,const MultFunctType &mult = standart_mult){
    T res = x;
    T factor = x;
    N-=1;
    while(N > 0){
        if(N%2)
            res = mult(res,factor);
        factor = mult(factor,factor);
        N = N/2;
    }
    return res;
}

#endif
