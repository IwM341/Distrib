
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
#include <grob/grid_objects.hpp>
#include <grob/csv_io.hpp>
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
    auto prog_path = CP.pgets("config_path");
    auto FPath = [&](const std::string & tree_path){
        return config_path_from(CP.pgets(tree_path),prog_path).string();
    };

    auto Mat0_SC = loadMatrix<T>(FPath("sc_in"),MF);

    size_t NN = std::get<1>(Mat0_SC);

    std::vector<T> Ev;
    if(CP.find("e_in")  != CP.not_found()){
        Ev = std::get<0>(loadMatrix<T>(FPath("e_in"),MF));
    }

    std::vector<T> D_I; //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP.find("d_init") != CP.not_found()){
        D_I = std::get<0>(loadMatrix<T>(FPath("d_init"),MF));
    }

    std::vector<T> C_I; //LoadVectorBinary<vtype>("D:/Desktop/del_files/matrix/H_capt.matrix");
    if(CP.find("c_init") != CP.not_found()){
        C_I = std::get<0>(loadMatrix<T>(FPath("c_init"),MF));
    }

    std::vector<T> C_Distrib  = std::get<0>(loadMatrix<T>(FPath("c_in"),MF));

    AnnCalcPolicy AP = ptree_contain(CP,"ann_ignore") ?
                                    AnnCalcPolicy::IGNORE:
                                    (
                                        ptree_contain(CP,"ann_full") ?
                                            AnnCalcPolicy::FULL :
                                            AnnCalcPolicy::ONLY_OUT
                                    );
    const int factor = (AP == AnnCalcPolicy::FULL ? 1 : 0);
    std::vector<T> ANN;
    if(CP.find("ann") != CP.not_found()){
        ANN = std::get<0>(loadMatrix<T>(FPath("ann"),MF));
    }

    auto & A_SC_T = std::get<0>(Mat0_SC);
    size_t ld_A_SC = NN;

    for(size_t i=0;i<NN;++i){
        A_SC_T[i+NN*i] = -std::accumulate(&A_SC_T[i*NN],&A_SC_T[i*NN]+NN,0.0);
    }
    for(size_t i=0;i<Ev.size();++i){
        A_SC_T[i+NN*i] -= Ev[i];
    }
    T max_sc = 0.0;
    for(size_t i=0;i<NN;++i){
        max_sc = std::max(A_SC_T[i+NN*i],max_sc);
    }

    T tau = CP.get<float>("tau",-1);
    int degree = CP.get<int>("deg",-1);
    T fullT = CP.get<float>("T",-1);


    T approp_tau = 1/max_sc;

    if(tau < 0){
        if(degree < 0){
            tau = approp_tau;
            if(fullT < 0)
                degree = 1;
            else
                degree = fullT/tau + 0.5;
        }
        else{
            if(fullT < 0){
                tau = approp_tau;
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
    std::string method = CP.pgets("method","euler");
    if(CP.get<std::string>("skip_pow","") == "auto"){
        //int pre_factor = (method == "euler" ? 10 : 1);
        sp = std::max(0,int(log2(10*tau/approp_tau)+1.5));
    }else{
        sp = CP.get<int>("skip_pow",0);
    }
    T tau_eff = (1<<sp);
    if(tau_eff> approp_tau){
        print("WARNING: tau too big");
        print("tau_max = ",approp_tau);
    }
    if(log_std){
        PVAR(tau);
        PVAR(tau/(1<<sp));
        PVAR(degree);
        PVAR(fullT);
    }


    if(log_std and MF == MatrixFormat::TEXT){
        printMatrix(NN,NN,A_SC_T.data());
        print();
    }


    //size_t iter_num = isVectorMethod ? degree : 1.5+log2(degree);
    /*
    auto RD = isVectorMethod ? VectorMultMethod<T>(NH,NL,A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,
                                                iter_num,tau) :
                               MatrixMultMethod<T>(NH,NL,A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,iter_num,tau);
    */
    std::vector<T> R;
    std::vector<T> SM;
    if(ptree_contain(CP,"load_state")){
        boost::filesystem::path backup_path = FPath("load_matrix");
        try{
            R = std::get<0>(loadMatrix<T>((backup_path / "R.bmat").string(),MF));
            SM = std::get<0>(loadMatrix<T>((backup_path / "I.bmat").string(),MF));
            goto Success;
        } catch(std::exception & e){
            std::cout << e.what()<< std::endl;
        }
        goto CreateMatrix;

    } else {
        CreateMatrix:

        R = (method == "order2" ? RMatrix2Order(NN,A_SC_T,tau_eff) :
                                   RMatrixEuler(NN,A_SC_T,tau_eff));
        SM = PowMatrixSumm(NN,R,sp,tau_eff);

        try{
            boost::filesystem::path backup_path = config_path_from(
                        CP.pgets("save_matrix","rs_matrix_backup")
                        ,prog_path);

            if (!boost::filesystem::exists(backup_path))
                boost::filesystem::create_directory(backup_path);


            saveMatrix<T>(R.data(),NN,NN,(backup_path / "R.bmat").string(),MF);

            saveMatrix<T>(SM.data(),NN,NN,(backup_path / "I.bmat").string(),MF);
        } catch (std::exception & e){
            print("problems saving matricies");
            print(e.what());

        }
    }
    Success:

    Result_F<T> RD = Evolution_E<T>( NN,R,SM,Ev,ANN,AP,
                                   C_Distrib,C_I,D_I,
                                    tau, fullT, sp);


    if(log_std and MF == MatrixFormat::TEXT){
        print("D_vec");
        print(RD.D_F);
    }
    //T S_LH_agg,S_HL_agg, E_H_agg, E_L_agg;


    if(log_std and MF == MatrixFormat::TEXT){
        print("C_distrib_vec");
        print(RD.C_F);
    }

    if(CP.get<bool>("norm",false)){
        auto normz = [](auto & x){x /= vector_sum(x);};
        normz(RD.C_F);
        normz(RD.D_F);
    }
    saveMatrix(RD.C_F.data(),1,NN,FPath("c_out"),MF);
    saveMatrix(RD.D_F.data(),1,NN,FPath("d_out"),MF);


    std::ofstream all_f(FPath("all_out"));
    print_values(all_f,"t[Ts]","C(t)","N(t)",
                 "E(t)","A(t)");
    print_vectors_tab(all_f,RD.T_grid,
                      RD.c_t,
                      RD.n_t,
                      RD.e_t,
                      RD.a_t);
}

int main(int argc, char ** argv){
    boost::property_tree::ptree CP;
    auto ret_key = parse_command_line_v1(argc,argv,CP);
    if(!ret_key.empty()){
        std::cout << "error : " << ret_key<<std::endl;
        return 0;
    }
    if(ptree_contain(CP,"help")){
        printd('\n',"params: ",
               "sc_in : [path to scatter matrix]",
               "c_in : [path to capture vector]",
               "e_in : [optional, path to evaporation vector]",
               "load_matrix : [optional, path to folder, where R matrix stores",
               "save_matrix : [optional, path to folder, where I matrix stores",
               "c_init : [optional, path to H full initial distri]",
               "d_init : [optional, path to H initial capture distrib]",
               "ann : [optional, path to annihilation matrix]",
               "a_cc : [optional, coeff matrix annihilation]",
               "ann_ignore : [optinal, flag or bool, if true then annihilation ignored]",
               "ann_full : [optinal, flag or bool, if true then annihilation is considered in equations]",
               "method : [optinal, euler - default or order2 - scheme in diff solver]",
               "T : [time to evolute]",
               "tau : [optional, step in equations]",
               "skip_pow : [optional, default is \"auto\" effective tau is tau/(2^sp)]");
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
