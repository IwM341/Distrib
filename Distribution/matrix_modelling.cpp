#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>

#include "func/matrix_functions.hpp"
#include <array>
#include <regex>
#include <string>
#include <string>

#include <cblas.h>

#include "func/arg_parser.hpp"
#include <sstream>
#include <boost/filesystem.hpp>


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
void  printPairVector(size_t NL,size_t NH,const T * L,const T*H){
    for(size_t il=0;il<NL;++il){
        std::cout << L[il] << "\t";
    }
    for(size_t ih=0;ih<NH;++ih){
        std::cout << H[ih] << "\t";
    }
    print();
}


#define pget get<std::string>
#define LOGIF(condition,expression) if(condition) PVAR(expression)
template <typename T>
void make_work(const boost::property_tree::ptree& CP){
    MatrixFormat MF = BINARY;
    if(CP.find("format") != CP.not_found()){
        auto fmt =CP.pget("format");
        if(fmt == "bin" || fmt == "binary" || fmt == "b")
            MF = BINARY;
    }
    bool log_std = ptree_condition(CP,"debug",false);
    LOGIF(log_std,MF);

    auto FPath = [&](const std::string & tree_path){
        return config_path_from(CP.pgets(tree_path),CP.pgets
                                ("config_path")).string();
    };

    auto Mat0_LH = loadMatrix<T>(FPath("lh_in"),MF);
    auto Mat0_HL = loadMatrix<T>(FPath("hl_in"),MF);

    size_t NH = std::get<1>(Mat0_HL);
    size_t NL = std::get<2>(Mat0_HL);

    std::vector<T> Ev_H;
    if(CP.find("eh_in") != CP.not_found()){
        Ev_H = std::get<0>(loadMatrix<T>(FPath("eh_in"),MF));
    }

    std::vector<T> Ev_L;
    if(CP.find("el_in")  != CP.not_found()){
        Ev_L = std::get<0>(loadMatrix<T>(FPath("el_in"),MF));
    }

    std::vector<T> EH_Distrib(NH,1.0); //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP.find("dh_in") != CP.not_found()){
        EH_Distrib = std::get<0>(loadMatrix<T>(FPath("dh_in"),MF));
    }

    //std::vector<T> EH_D_out(EH_Distrib.size());

    std::vector<T> EL_Distrib(NL,1.0);
    if(CP.find("dl_in") != CP.not_found()){
        EL_Distrib = std::get<0>(loadMatrix<T>(FPath("dl_in"),MF));
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
        tau = CP.get<float>("tau");
    }  catch (std::exception & e) {
        tau = 0.05/std::max(-min(A_HH),-min(A_LL));
        PVAR(tau);
    }
    LOGIF(log_std,tau);

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
        it_number = CP.get<int>("fac2");
    }  catch (std::exception & e) {
        it_number = 10;
    }

    Function::GridFunction<T,Function::VectorGrid<T>,Function::LinearInterpolator> LN_t(
                Vector(it_number+1,[tau](size_t i)->T{return (i==0 ? 0 :tau*pow(2.0,i));}));
    auto HN_t = LN_t;
    LN_t[0] = vector_sum(EL_Distrib);
    HN_t[0] = vector_sum(EH_Distrib);

    auto LD_tmp = EL_Distrib;
    auto HD_tmp = EH_Distrib;

    /**/
    if(log_std)
        printPairMatrix(NL,NH,A_LL.data(),A_HL_T.data(),A_LH_T.data(),A_HH.data());

    for(size_t i=1;i<=it_number;++i){
        std::cout << i <<std::endl;
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
        if(EL_Distrib.size() && EH_Distrib.size() && LD_tmp.size() && LD_tmp.size()){
            gemv<T>(CblasColMajor,CblasNoTrans,NL,NL,1.0,A_LL.data(),ld_A_LL,EL_Distrib.data(),1,0.0,LD_tmp.data(),1);
            gemv<T>(CblasColMajor,CblasNoTrans,NL,NH,1.0,A_HL_T.data(),ld_A_HL,EH_Distrib.data(),1,1.0,LD_tmp.data(),1);

            gemv<T>(CblasColMajor,CblasNoTrans,NH,NL,1.0,A_LH_T.data(),ld_A_LH,EL_Distrib.data(),1,0.0,HD_tmp.data(),1);
            gemv<T>(CblasColMajor,CblasNoTrans,NH,NH,1.0,A_HH.data(),ld_A_HH,EH_Distrib.data(),1,1.0,HD_tmp.data(),1);
        }
        LN_t[i] = vector_sum(LD_tmp);
        HN_t[i] = vector_sum(HD_tmp);
    }
    if(log_std){
        printPairVector(NL,NH,LD_tmp.data(),HD_tmp.data());
    }
    auto LD_summ = vector_sum(LD_tmp);
    auto HD_summ = vector_sum(LD_tmp);
    if(LD_summ == 0.0){
        std::cout << "WARNING, LD summ =0" <<std::endl;
    }
    if(HD_summ == 0.0){
        std::cout << "WARNING, HD summ =0" <<std::endl;
    }
    LD_tmp *= (1.0/LD_summ);
    HD_tmp *= (1.0/HD_summ);

    if(log_std){
        printPairVector(NL,NH,LD_tmp.data(),HD_tmp.data());
    }

    saveMatrix(LD_tmp.data(),1,NL,FPath("dl_out"),MF);
    saveMatrix(HD_tmp.data(),1,NH,FPath("dh_out"),MF);

    if(CP.find("nl_out") != CP.not_found()){
        std::ofstream LN_f(FPath("nl_out"));
        LN_f <<LN_t.toString();
    }
    if(CP.find("nh_out")  != CP.not_found()){
        std::ofstream HN_f(FPath("nh_out"));
        HN_f <<HN_t.toString();
    }
}

int main(int argc, char ** argv){
    boost::property_tree::ptree CP;
    auto ret_key = parse_command_line(argc,argv,CP);
    if(!ret_key.empty()){
        std::cout << "error : " << ret_key<<std::endl;
        return 0;
    }

    std::map<std::string,std::string> default_params{{"dl_out","dl_out.mat"},{"dh_out","dl_out.mat"}};
    std::vector<std::string> required_params{"lh_in","hl_in"};
    for(const auto & param : required_params){
        if(CP.find(param)==CP.not_found()){
            std::cout << "required param " << param << std::endl;
            return 0;
        }
    }
//    for(size_t i=0;i<CP.size();++i){
//        std::cout << key_names[i] << " : " << CP[i] << std::endl;
//    }
    if(CP.find("type")!=CP.not_found() && CP.pget("type") == "float")
        make_work<float>(CP);
    else
        make_work<double>(CP);

    return 0;
}
