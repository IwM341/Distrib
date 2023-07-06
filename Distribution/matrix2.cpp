
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

    auto Mat0_LH = loadMatrix<T>(FPath("lh_in"),MF);
    auto Mat0_HL = loadMatrix<T>(FPath("hl_in"),MF);

    size_t NL = std::get<1>(Mat0_HL);
    size_t NH = std::get<2>(Mat0_HL);

    std::vector<T> LL(NL*NL,0),HH(NH*NH,0);
    try{
        LL = std::get<0>(loadMatrix<T>(FPath("ll_in"),MF));
    }catch(std::exception & load_LL_exception){
        PVAR(load_LL_exception.what());
    }
    try{
        HH = std::get<0>(loadMatrix<T>(FPath("hh_in"),MF));
    }catch(std::exception & load_HH_exception){
        PVAR(load_HH_exception.what());
    }



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

    AnnCalcPolicy AP = ptree_contain(CP,"ann_ignore") ?
                                    AnnCalcPolicy::IGNORE:
                                    (
                                        ptree_contain(CP,"ann_full") ?
                                            AnnCalcPolicy::FULL :
                                            AnnCalcPolicy::ONLY_OUT
                                    );
    const int factor = (AP == AnnCalcPolicy::FULL ? 1 : 0);
    std::vector<T> ANN_HL;
    if(CP.find("ann_hl") != CP.not_found()){
        ANN_HL = std::get<0>(loadMatrix<T>(FPath("ann_hl"),MF))*
                CP.get<float>("a_hl",1);
    }
    std::vector<T> ANN_LL;
    if(CP.find("ann_ll") != CP.not_found()){
        ANN_LL = std::get<0>(loadMatrix<T>(FPath("ann_ll"),MF))*
                CP.get<float>("a_ll",1);
    }
    std::vector<T> ANN_HH;
    if(CP.find("ann_hh") != CP.not_found()){
        ANN_HH = std::get<0>(loadMatrix<T>(FPath("ann_hh"),MF))*
                CP.get<float>("a_hh",1);
    }
    std::vector<T> ANN_L;
    if(CP.find("ann_el") != CP.not_found()){
        ANN_L = std::get<0>(loadMatrix<T>(FPath("ann_el"),MF))*
                CP.get<float>("a_el",1);
    }
    std::vector<T> ANN_H;
    if(CP.find("ann_eh") != CP.not_found()){
        ANN_H = std::get<0>(loadMatrix<T>(FPath("ann_eh"),MF))*
                CP.get<float>("a_eh",1);
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
        LL[i+NL*i] += -std::accumulate(&A_HL_T[i*NH],&A_HL_T[i*NH]+NH,0.0) -
                std::accumulate(&LL[i*NL],&LL[i*NL]+NL,0.0);
    }
    for(size_t i=0;i<Ev_L.size();++i){
        LL[i+NL*i] -= Ev_L[i];
    }
    for(size_t i=0;i<ANN_L.size();++i){
        LL[i+NL*i] -= ANN_L[i]*factor;
    }


    for(size_t i=0;i<NH;++i){
        HH[i+NH*i] += -std::accumulate(&A_LH_T[i*NL],&A_LH_T[i*NL]+NL,0.0)-
                std::accumulate(&HH[i*NH],&HH[i*NH]+NH,0.0);
    }
    for(size_t i=0;i<Ev_H.size();++i){
        HH[i+NH*i] -= Ev_H[i];
    }
    for(size_t i=0;i<ANN_H.size();++i){
        HH[i+NH*i] -= ANN_H[i]*factor;
    }
    T max_sc = 0;
    for(size_t i=0;i<NL;++i){
        max_sc = std::max(max_sc,-LL[i+NL*i]);
    }
    for(size_t i=0;i<NH;++i){
        max_sc = std::max(max_sc,-HH[i+NH*i]);
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
                tau = fullT/degree;
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

    if(CP.get<std::string>("skip_pow","auto") == "auto"){
        //int pre_factor = (method == "euler" ? 10 : 1);
        sp = std::max(0,int(log2(10*tau/approp_tau)+1.5));
    }else{
        sp = CP.get<int>("skip_pow");
    }
    T tau_eff = tau/pow(2,sp);
    if(tau_eff> approp_tau){
        print("WARNING: tau too big");
        print("tau_max = ",approp_tau);
    }
    if(log_std){
        PVAR(tau);
        PVAR(sp);
        PVAR(tau/pow(2,sp));
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
        printPairMatrix(NL,NH,LL.data(),A_LH_T.data(),A_HL_T.data(),HH.data());
        print();
    }


    //size_t iter_num = isVectorMethod ? degree : 1.5+log2(degree);
    /*
    auto RD = isVectorMethod ? VectorMultMethod<T>(NH,NL,A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,
                                                iter_num,tau) :
                               MatrixMultMethod<T>(NH,NL,A_LL,A_LH_T,A_HL_T,A_HH,L_Distrib,H_Distrib,iter_num,tau);
    */
    BlockMatrix<T> R;
    BlockMatrix<T> SM;
    if(ptree_contain(CP,"load_matrix")){
        boost::filesystem::path backup_path = FPath("load_matrix");
        try{
            R.LL = std::get<0>(loadMatrix<T>((backup_path / "R_LL.bmat").string(),MF));
            R.LH = std::get<0>(loadMatrix<T>((backup_path / "R_LH.bmat").string(),MF));
            R.HL = std::get<0>(loadMatrix<T>((backup_path / "R_HL.bmat").string(),MF));
            R.HH = std::get<0>(loadMatrix<T>((backup_path / "R_HH.bmat").string(),MF));

            SM.LL = std::get<0>(loadMatrix<T>((backup_path / "I_LL.bmat").string(),MF));
            SM.LH = std::get<0>(loadMatrix<T>((backup_path / "I_LH.bmat").string(),MF));
            SM.HL = std::get<0>(loadMatrix<T>((backup_path / "I_HL.bmat").string(),MF));
            SM.HH = std::get<0>(loadMatrix<T>((backup_path / "I_HH.bmat").string(),MF));
            goto Success;
        } catch(std::exception & e){
            std::cout << e.what()<< std::endl;
        }
        goto CreateMatrix;

    } else {
        CreateMatrix:

        R =  RMatrix2Order2(NL,NH,LL,A_LH_T,A_HL_T,HH,tau_eff);
        SM = PowMatrixSumm(NL,NH,R.LL,R.LH,R.HL,R.HH,sp,tau_eff);

        try{
            boost::filesystem::path backup_path = config_path_from(
                        CP.pgets("save_matrix","rs_matrix_backup")
                        ,prog_path);

            if (!boost::filesystem::exists(backup_path))
                boost::filesystem::create_directory(backup_path);


            saveMatrix<T>(R.LL.data(),NL,NL,(backup_path / "R_LL.bmat").string(),MF);
            saveMatrix<T>(R.LH.data(),NL,NH,(backup_path / "R_LH.bmat").string(),MF);
            saveMatrix<T>(R.HL.data(),NH,NL,(backup_path / "R_HL.bmat").string(),MF);
            saveMatrix<T>(R.HH.data(),NH,NH,(backup_path / "R_HH.bmat").string(),MF);

            saveMatrix<T>(SM.LL.data(),NL,NL,(backup_path / "I_LL.bmat").string(),MF);
            saveMatrix<T>(SM.LH.data(),NL,NH,(backup_path / "I_LH.bmat").string(),MF);
            saveMatrix<T>(SM.HL.data(),NH,NL,(backup_path / "I_HL.bmat").string(),MF);
            saveMatrix<T>(SM.HH.data(),NH,NH,(backup_path / "I_HH.bmat").string(),MF);
        } catch (std::exception & e){
            print("problems saving matricies");
            print(e.what());

        }
    }
    Success:

    ResultLH_F<T> RD = Evolution_LH<T>( NL, NH,R,SM,
                                        Ev_L,Ev_H,ANN_L,ANN_H,
                                        ANN_LL,ANN_HL,ANN_HH,
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

    if(CP.get<bool>("norm",false)){
        auto normz = [](auto & x){x /= vector_sum(x);};
        normz(RD.C_L_F);
        normz(RD.C_H_F);
        normz(RD.D_L_F);
        normz(RD.D_H_F);
    }
    saveMatrix(RD.C_L_F.data(),1,NL,FPath("cl_out"),MF);
    saveMatrix(RD.C_H_F.data(),1,NH,FPath("ch_out"),MF);
    saveMatrix(RD.D_L_F.data(),1,NL,FPath("dl_out"),MF);
    saveMatrix(RD.D_H_F.data(),1,NH,FPath("dh_out"),MF);


    std::ofstream all_f(FPath("evolution_out"));
    print_values(all_f,"t[Ts]","C_l(t)","C_h(t)","N_l(t)","N_h(t)",
                 "E_l(t)","E_h(t)","A_LL(t)","A_HL(t)","A_HH(t)","A_L(t)","A_H(t)","A(t)");
    print_vectors_tab(all_f,RD.T_grid,
                      RD.c_l_t,RD.c_h_t,
                      RD.n_l_t,RD.n_h_t,
                      RD.e_l_t,RD.e_h_t,
                      RD.a_ll_t,RD.a_hl_t,RD.a_hh_t,
                      RD.a_l_t,RD.a_h_t,
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
               "ll_in : [optional, path to scatter matrix from H to L]",
               "hh_in : [optional, to scatter matrix from L to H]",
               "lh_in : [path to scatter matrix from H to L]",
               "hl_in : [path to scatter matrix from L to H]",
               "cl_in : [path to L capture vector]",
               "ch_in : [path to H capture vector]",
               "cl_out : [normilized output captured L distribution]",
               "ch_out : [normilized output captured H distribution]",
               "dl_out : [normilized output full L distribution]",
               "dh_out : [normilized output full L distribution]",
               "norm : [optional, if true, result normalizet to 1]",
               "evolution_out : [output dat file of evolution]",
               "l_nuf : [coeff before L capture (or fraction of H in halo)]",
               "h_nuf : [coeff before H capture (or fraction of L in halo)]",
               "el_in : [optional, path to L evaporation matrix]",
               "eh_in : [optional, path to H evaporation matrix]",
               "load_matrix : [optional, path to folder, where R matrix stores",
               "save_matrix : [optional, path to folder, where I matrix stores",
               "cl_init : [optional, path to L full initial distrib]",
               "ch_init : [optional, path to H full initial distri]",
               "dl_init : [optional, path to L initial capture distrib]",
               "dh_init : [optional, path to H initial capture distrib]",
               "ann_ll : [optional, path to LL annihilation matrix]",
               "ann_hl : [optional, path to HL annihilation matrix]",
               "ann_hh : [optional, path to HH annihilation matrix]",
               "ann_el : [optional, path to L annihilation vector]",
               "ann_eh : [optional, path to H annihilation vector]",
               "a_ll : [optional, coeff LL annihilation]",
               "a_hl : [optional, coeff HL annihilation]",
               "a_hh : [optional, coeff HH annihilation]",
               "a_el : [optional, coeff L vector annihilation]",
               "a_eh : [optional, coeff H vector annihilationr]",
               "a_eh : [optional, coeff H vector annihilationr]",
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
