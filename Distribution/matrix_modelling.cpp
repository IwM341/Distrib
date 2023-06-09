
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

template <typename Arg>
void print_values(std::ostream & os,Arg const & value){
    os << value <<std::endl;
}

template <typename Arg,typename...Args>
void print_values(std::ostream & os,Arg const & value,Args const&...args){
    os << value << "\t";
    print_values(os,args...);
}

template <typename Vector,typename...Vectors>
void print_vectors_tab(std::ostream & os,Vector const & V,Vectors const &...Vs){
    for(size_t i=0;i<V.size();++i){
        print_values(os,V[i],Vs[i]...);
    }
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

    std::vector<T> H_D_I; //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP.find("dh_init") != CP.not_found()){
        H_D_I = std::get<0>(loadMatrix<T>(FPath("dh_init"),MF));
    }
    std::vector<T> L_D_I;
    if(CP.find("dh_init") != CP.not_found()){
        L_D_I = std::get<0>(loadMatrix<T>(FPath("dl_init"),MF));
    }
    std::vector<T> H_C_I; //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP.find("ch_init") != CP.not_found()){
        H_C_I = std::get<0>(loadMatrix<T>(FPath("ch_init"),MF));
    }
    std::vector<T> L_C_I;
    if(CP.find("ch_init") != CP.not_found()){
        L_C_I = std::get<0>(loadMatrix<T>(FPath("cl_init"),MF));
    }


    std::vector<T> C_H_Distrib  = std::get<0>(loadMatrix<T>(FPath("ch_in"),MF))*
            CP.get<T>("h_nuf",1);

    std::vector<T> C_L_Distrib = std::get<0>(loadMatrix<T>(FPath("cl_in"),MF))*
            CP.get<T>("l_nuf",1);

    std::vector<T> ANN_HL;
    if(CP.find("ann_hl") != CP.not_found()){
        ANN_HL = std::get<0>(loadMatrix<T>(FPath("ann_hl"),MF))*
                CP.get<float>("ann_coeff",1);
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
    std::string method = CP.pgets("method","euler");
    AnnCalcPolicy AP = ptree_contain(CP,"ann_hl") ?
                            (
                                ptree_contain(CP,"ann_ignore") ?
                                    AnnCalcPolicy::IGNORE:
                                    (
                                        ptree_contain(CP,"ann_full") ?
                                            AnnCalcPolicy::FULL :
                                            AnnCalcPolicy::ONLY_OUT
                                    )
                            )
                        : AnnCalcPolicy::IGNORE;
    ResultLH_F<T> RD = Evolution_LH<T>( NL, NH,A_LL, A_LH_T,
                                       A_HL_T,A_HH,
                                          method,
                                          ANN_HL,
                                            AP,
                                          C_L_Distrib,C_H_Distrib,
                                          L_C_I,H_C_I,
                                          L_D_I,H_D_I,
                                          tau, fullT, sp);


    if(log_std and MF == MatrixFormat::TEXT){
        print("D_vecs");
        printPairVector(NL,NH,RD.D_L_F.data(),RD.D_H_F.data());
    }
    //T S_LH_agg,S_HL_agg, E_H_agg, E_L_agg;


    if(log_std and MF == MatrixFormat::TEXT){
        print("C_distrib_vecs");
        printPairVector(NL,NH,RD.C_L_F.data(),RD.C_H_F.data());
    }

    saveMatrix(RD.C_L_F.data(),1,NL,FPath("cl_out"),MF);
    saveMatrix(RD.C_H_F.data(),1,NH,FPath("ch_out"),MF);
    saveMatrix(RD.D_L_F.data(),1,NL,FPath("dl_out"),MF);
    saveMatrix(RD.D_H_F.data(),1,NH,FPath("dh_out"),MF);


    std::ofstream all_f(FPath("all_out"));
    print_values(all_f,"t[Ts]","C_l(t)","C_h(t)","N_l(t)","N_h(t)","A(t)");
    print_vectors_tab(all_f,RD.T_grid,
                      RD.c_l_t,RD.c_h_t,
                      RD.n_l_t,RD.n_h_t,
                      RD.a_t);
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
