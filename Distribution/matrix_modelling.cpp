#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
#include "func/el_histo.hpp"
#include "func/matrix_functions.hpp"
#include <array>
#include <regex>
#include <string>
#include <string>

#include <cblas.h>

typedef  double vtype;
#define vgemm cblas_dgemm


template <typename T>
constexpr auto gemm = std::conditional<std::is_same<T,float>::value,
                                    std::integral_constant<decltype((cblas_sgemm)),cblas_sgemm>,
                                    std::integral_constant<decltype((cblas_dgemm)),cblas_dgemm>>::type::value;
template <typename T>
constexpr auto gemv = std::conditional<std::is_same<T,float>::value,
                                    std::integral_constant<decltype((cblas_sgemv)),cblas_sgemv>,
                                    std::integral_constant<decltype((cblas_dgemv)),cblas_dgemv>>::type::value;

struct index{
    size_t N,M;
    bool traspose;
    index():N(0),M(0),traspose(0){}
    index(size_t N,size_t M):N(N),M(M),traspose(0){}
    size_t operator()(size_t n,size_t m){
        if(traspose){
            return m*N+n;
        }
        else{
            return n*M+m;
        }
    }
};

/*
struct config_params{
    std::string LH_in;
    std::string HL_in;
    std::string EH_in;
    std::string EL_in;
    std::string DL_in;
    std::string DH_in;

    std::string DL_out;
    std::string DH_out;

    std::string HN_out;
    std::string LN_out;
};
*/

typedef std::array<std::string,14> config_params;
const std::array<bool,14> is_optional_param = {0,0,1,1,1,1,0,0,1,1,1,1,1,1};

enum ErrorCode{
    PARSE_OK = 0,
    PARSE_TOO_FEW_PARAMETRS = 16,
    PARSE_DUBLICATION_OF_PARAMETRS = 8,
    PARSE_UNDEFINED_KEY = 4,
    PARSE_EXPECT_VALUE = 2,
    PARSE_NO_SUCH_FILE = 1
};

enum KeyParam{
    KEY_LH_IN,
    KEY_HL_IN,
    KEY_EL_IN,
    KEY_EH_IN,

    KEY_DL_IN,
    KEY_DH_IN,
    KEY_DL_OUT,
    KEY_DH_OUT,
    KEY_NL_OUT,
    KEY_NH_OUT,
    KEY_TYPE,
    KEY_TAU,
    KEY_FAC2,
    KEY_ERROR
};
const std::array<std::string,14> key_names = {"lh_in","hl_in","el_in","eh_in",
                                             "dl_in","dh_in","dl_out","dh_out",
                                             "nl_out","nh_out","type","tau","fac2","error"};
char my_toupper(char ch)
{
    return static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
}
std::string str_toupper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                // static_cast<int(*)(int)>(std::toupper)         // wrong
                // [](int c){ return std::toupper(c); }           // wrong
                // [](char c){ return std::toupper(c); }          // wrong
                   [](unsigned char c){ return std::toupper(c); } // correct
                  );
    return s;
}
KeyParam KeyStr(const std::string&par){
    auto S = str_toupper(par);
    if(S == "-LH_IN")
        return KEY_LH_IN;
    else if(S == "-HL_IN")
        return KEY_HL_IN;
    else if(S == "-EL_IN")
        return KEY_EL_IN;
    else if(S == "-EH_IN")
        return KEY_EH_IN;
    else if(S == "-DL_IN")
        return KEY_DL_IN;
    else if(S == "-DH_IN")
        return KEY_DH_IN;
    else if(S == "-DL_OUT")
        return KEY_DL_OUT;
    else if(S == "-DH_OUT")
        return KEY_DH_OUT;
    else if(S == "-NL_OUT")
        return KEY_NL_OUT;
    else if(S == "-NH_OUT")
        return KEY_NH_OUT;
    else if(S == "-TYPE")
        return KEY_TYPE;
    else if(S == "-TAU")
        return KEY_TAU;
    else if(S == "-FAC2")
        return KEY_FAC2;
    else
        return KEY_ERROR;
}

int parse_config_file(const std::string &filename,config_params & cp){
    std::regex jsr("\"(.*)\"\\s*:\\s*\"(.*)\"");
    std::smatch sm;
    std::ifstream input_file(filename);
    if(!input_file.is_open())
        return PARSE_NO_SUCH_FILE;
    auto S = std::string((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
    while(regex_search(S, sm, jsr))
    {
        auto key = KeyStr("-"+sm[1].str());
        if(key == KEY_ERROR)
            return PARSE_UNDEFINED_KEY;
        cp[key] = sm[2];
        S = sm.suffix();
    }
    return PARSE_OK;
}
int parse_cmd(int argc,char ** argv,config_params & CP){
    for(size_t i=1;i<argc-1;i+=2){
        if(std::string(str_toupper(argv[i])) == "-CONFIG"){
            auto RetCode = parse_config_file(argv[i+1],CP);
            if(RetCode!=0){
                return RetCode;
            }
        }   else{

            auto key = KeyStr(argv[i]);
            if(key == KEY_ERROR){
                return PARSE_UNDEFINED_KEY;
            }
            else{
                CP[key] = argv[i+1];
            }
        }
    }
    for(size_t i=0;i<CP.size();++i){
        if(!is_optional_param[i] && !CP[i].size()){
            return PARSE_TOO_FEW_PARAMETRS + i;
        }
    }
    return 0;
}


template <typename T>
void  printPairMatrix(size_t NL,size_t NH,const T * LL,const T*HL,const  T*LH,const T*HH){
    for(size_t il=0;il<NL;++il){
        for(size_t jl=0;jl<NL;++jl){
            std::cout << LL[il+jl*NL] << "\t";
        }
        for(size_t jh=0;jh<NH;++jh){
            std::cout << HL[il+jh*NL] << "\t";
        }std::cout << std::endl;
    }
    for(size_t ih=0;ih<NH;++ih){
        for(size_t jl=0;jl<NL;++jl){
            std::cout << LH[ih+jl*NH] << "\t";
        }
        for(size_t jh=0;jh<NH;++jh){
            std::cout << HH[ih+jh*NH] << "\t";
        }std::cout << std::endl;
    }
}

template <typename T>
void make_work(const config_params & CP){
    auto Mat0_LH = loadMatrixSmart<T>(CP[KEY_LH_IN]);
    auto Mat0_HL = loadMatrixSmart<T>(CP[KEY_HL_IN]);

    size_t NH = std::get<1>(Mat0_HL);
    size_t NL = std::get<2>(Mat0_HL);

    std::vector<T> Ev_H;
    if(CP[KEY_EH_IN].size()){
        Ev_H = std::get<0>(loadMatrixSmart<T>(CP[KEY_EH_IN]));
    }

    std::vector<T> Ev_L;
    if(CP[KEY_EL_IN].size())
        Ev_L = std::get<0>(loadMatrixSmart<T>(CP[KEY_EL_IN]));

    std::vector<T> EH_Distrib(NH,1.0); //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP[KEY_DH_IN].size()){
        EH_Distrib = std::get<0>(loadMatrixSmart<T>(CP[KEY_DH_IN]));
    }

    //std::vector<T> EH_D_out(EH_Distrib.size());

    std::vector<T> EL_Distrib(NL,1.0);
    if(CP[KEY_DL_IN].size()){
        EL_Distrib = std::get<0>(loadMatrixSmart<T>(CP[KEY_DL_IN]));
    }

    //std::vector<T> EL_D_out(EL_Distrib.size());

    //auto EL_Distrib = EH_Distrib*0;



    auto & A_HL_T = std::get<0>(Mat0_HL);
    size_t ld_A_HL = NL;

    auto & A_LH_T = std::get<0>(Mat0_LH);
    size_t ld_A_LH = NH;

    std::vector<T> A_LL(NL*NL,0);
    size_t ld_A_LL = NL;
    for(size_t i=0;i<NL;++i){
        A_LL[i*NL+i] = -std::accumulate(&A_LH_T[i*NH],&A_LH_T[i*NH]+NH,0.0);
        if(Ev_L.size())
            A_LL[i*NL+i] -= Ev_L[i];
    }
    std::vector<T> A_HH(NH*NH,0);
    size_t ld_A_HH = NH;
    for(size_t i=0;i<NH;++i){
        A_HH[i*NH+i] = -std::accumulate(&A_HL_T[i*NL],&A_HL_T[i*NL]+NL,0.0);
        if(Ev_H.size())
            A_HH[i*NH+i] -= Ev_H[i];
    }

    T tau;
    try {
        tau = std::stod(CP[KEY_TAU]);
    }  catch (std::exception & e) {
        tau = 0.05/std::max(-min(A_HH),-min(A_LL));
    }


    A_LL*=tau;
    A_HH*=tau;
    A_HL_T*=tau;
    A_LH_T*=tau;

    for(size_t i=0;i<NL;++i){
        A_LL[i*NL+i] += 1;
    }
    for(size_t i=0;i<NH;++i){
        A_HH[i*NH+i] += 1;
    }

    std::vector<T>B_HL_T(A_HL_T.size());
    std::vector<T>B_LH_T(A_LH_T.size());
    std::vector<T>B_LL(A_LL.size());
    std::vector<T>B_HH(A_HH.size());
    //print();

    size_t it_number = 10;

    try {
        it_number = std::stoi(CP[KEY_FAC2]);
    }  catch (std::exception & e) {
        it_number = 10;
    }

    Function::GridFunction<T,Function::VectorGrid<T>,Function::LinearInterpolator> LN_t(
                Vector(it_number+1,[](size_t i)->T{return pow(2.0,i);}));
    auto HN_t = LN_t;
    LN_t[0] = vector_sum(EL_Distrib);
    HN_t[0] = vector_sum(EH_Distrib);

    auto LD_tmp = EL_Distrib;
    auto HD_tmp = EH_Distrib;

    /**/
    //printPairMatrix(NL,NH,A_LL.data(),A_HL_T.data(),A_LH_T.data(),A_HH.data());
    for(size_t i=1;i<=it_number;++i){
        std::cout << it_number <<std::endl;
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NL,
                    1.0,A_LL.data(),ld_A_LL,A_LL.data(),ld_A_LL,0.0,B_LL.data(),ld_A_LL);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NL,NH,
                            1.0,A_HL_T.data(),ld_A_HL,A_LH_T.data(),ld_A_LH,1.0,B_LL.data(),ld_A_LL);


        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NL,
                    1.0,A_LL.data(),ld_A_LL,A_HL_T.data(),ld_A_HL,0.0,B_HL_T.data(),ld_A_HL);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NL,NH,NH,
                            1.0,A_HL_T.data(),ld_A_HL,A_HH.data(),ld_A_HH,1.0,B_HL_T.data(),ld_A_HL);


        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NL,
                    1.0,A_LH_T.data(),ld_A_LH,A_LL.data(),ld_A_LL,0.0,B_LH_T.data(),ld_A_LH);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NL,NH,
                            1.0,A_HH.data(),ld_A_HH,A_LH_T.data(),ld_A_LH,1.0,B_LH_T.data(),ld_A_LH);

        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NL,
                    1.0,A_LH_T.data(),ld_A_LH,A_HL_T.data(),ld_A_HL,0.0,B_HH.data(),ld_A_HH);
        gemm<T>(CblasColMajor,CblasNoTrans,CblasNoTrans,NH,NH,NH,
                            1.0,A_HH.data(),ld_A_HH,A_HH.data(),ld_A_HH,1.0,B_HH.data(),ld_A_HH);

        std::swap(A_LL,B_LL);
        std::swap(A_HL_T,B_HL_T);
        std::swap(A_LH_T,B_LH_T);
        std::swap(A_HH,B_HH);

        //printPairMatrix(NL,NH,A_LL.data(),A_HL_T.data(),A_LH_T.data(),A_HH.data());
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,A_LL.data(),ld_A_LL,EL_Distrib.data(),1,0.0,LD_tmp.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,A_HL_T.data(),ld_A_HL,EH_Distrib.data(),1,1.0,LD_tmp.data(),1);

        gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,A_LH_T.data(),ld_A_LH,EL_Distrib.data(),1,0.0,HD_tmp.data(),1);
        gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,A_HH.data(),ld_A_HH,EH_Distrib.data(),1,1.0,HD_tmp.data(),1);

        LN_t[i] = vector_sum(LD_tmp);
        HN_t[i] = vector_sum(HD_tmp);
    }
    LD_tmp *= (1.0/vector_sum(LD_tmp));
    HD_tmp *= (1.0/vector_sum(HD_tmp));
    saveMatrixSmart(LD_tmp.data(),1,NL,CP[KEY_DL_OUT]);
    saveMatrixSmart(HD_tmp.data(),1,NH,CP[KEY_DH_OUT]);

    if(CP[KEY_NL_OUT].size()){
        std::ofstream LN_f(CP[KEY_NL_OUT]);
        LN_f <<LN_t.toString();
    }
    if(CP[KEY_NH_OUT].size()){
        std::ofstream HN_f(CP[KEY_NH_OUT]);
        HN_f <<HN_t.toString();
    }
}

int main(int argc, char ** argv){
    config_params  CP;
    auto ret_key = parse_cmd(argc,argv,CP);
    if(ret_key !=0){
        std::cout << "error : " << ret_key<<std::endl;
        if(ret_key >= PARSE_TOO_FEW_PARAMETRS){
            std::cout << "need parametr " << key_names[ret_key-PARSE_TOO_FEW_PARAMETRS]<<std::endl;
        }
        return 0;
    }

//    for(size_t i=0;i<CP.size();++i){
//        std::cout << key_names[i] << " : " << CP[i] << std::endl;
//    }
    if(CP[KEY_TYPE] == "float")
        make_work<float>(CP);
    else
        make_work<double>(CP);



    /*
    PrintMatrix(B_LL,std::cout);
    PrintMatrix(B_HH,std::cout);
    PVAR(MT_LH_out);
    PVAR(MT_HL_out);
    PVAR(Ev_H);
    PVAR(Ev_L);


    double Nu_Max = std::max(mnorm_1(MT_LH),mnorm_1(MT_HL));
    double tau = 0.01;
    if(tau*Nu_Max > 0.1)
        tau = 0.1/Nu_Max;

    dmatrix LLmat,HLmat;
    dmatrix LHmat,HHmat;

    auto E1 = E_Matrix<double>(Ev_L.size());
    _R(LLmat,HLmat,LHmat,HHmat) = smart_pow(_T(E1-tau*MatrixDiagonal(MT_LH_out),tau*MT_HL,
                                               tau*MT_LH,E1-tau*MatrixDiagonal(MT_HL_out)),1000000,
                                            [](const auto &T1,const auto &T2){
        const auto & A11 = std::get<0>(T1);
        const auto & A12 = std::get<1>(T1);
        const auto & A21 = std::get<2>(T1);
        const auto & A22 = std::get<3>(T1);

        const auto & B11 = std::get<0>(T2);
        const auto & B12 = std::get<1>(T2);
        const auto & B21 = std::get<2>(T2);
        const auto & B22 = std::get<3>(T2);

        return  _T(MatrixMult(A11,B11)+MatrixMult(A12,B21),MatrixMult(A11,B12)+MatrixMult(A12,B22),
                   MatrixMult(A21,B11)+MatrixMult(A22,B21),MatrixMult(A21,B12)+MatrixMult(A22,B22));
    });

    std::ofstream LLf("D:/Desktop/del_files/distrib/LL.txt");
    PrintMatrix(LLmat,LLf);
    LLf.close();

    std::ofstream HLf("D:/Desktop/del_files/distrib/HL.txt");
    PrintMatrix(HLmat,HLf);
    HLf.close();

    std::ofstream LHf("D:/Desktop/del_files/distrib/LH.txt");
    PrintMatrix(LHmat,LHf);
    LHf.close();

    std::ofstream HHf("D:/Desktop/del_files/distrib/HH.txt");
    PrintMatrix(HHmat,HHf);
    HHf.close();

    PVAR(tau);

    double T = 1000;

    size_t N_pass = 1000;
    size_t Nt = T/tau;

    Function::UniformGrid<double> T_grid(0,tau*Nt,Nt/N_pass+1);
    Function::GridFunction<double,decltype(T_grid),Function::LinearInterpolator> N_L(T_grid);
    decltype(N_L) N_H(T_grid);


    show_for_progress(T_grid.size(),[&](size_t i)
    for(size_t i=0;i<T_grid.size();++i){
        N_L[i] = vector_sum(EL_Distrib);
        N_H[i] = vector_sum(EH_Distrib);
        for(size_t j=0;j<N_pass;++j)
            _R(EL_Distrib,EH_Distrib) = R_Func(EL_Distrib,EH_Distrib,tau);

        double Hmin = min(EH_Distrib);
        double Lmin = min(EL_Distrib);
        if(Hmin <0 or Lmin < 0){
            PVAR(Hmin);
            PVAR(Lmin);
            PVAR(T_grid[i]);
        }
    });
    auto Integrate = [](const std::vector<double> &V,double h){
        return h*(vector_sum(V)-0.5*(V[0]+V.back()));};


    double SH = Integrate(N_H.values,T_grid._h());
    double SL = Integrate(N_L.values,T_grid._h());
    PVAR(SH);
    PVAR(SL);

    Gnuplot gp;
    gp.plotd(N_H.toString(),"with lines title \"heavy\"");
    gp.plotd(N_L.toString(),"with lines title \"light\"");
    gp.show();

    std::cin.get();
    */
    return 0;
}
