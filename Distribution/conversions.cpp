#include <iostream>
#include "functions.hpp"
#include <cmath>
#include <fstream>

#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/arg_parser.hpp"
#include "func/dens_grid.hpp"
#include "func/matrix_functions.hpp"
#include "func/load_histo.hpp"
#include "csv_io.hpp"

//#include "boost/type_index.hpp"
/*
std::string process_space(std::string const &S){
    std::regex fname("\\s.(\"([\\w].)\")?(\\w.)\\s.");
}
*/


template <typename HistoType,typename FuncType>
double EL_dens(const HistoType &H,FuncType && LE,char method,double E,double L){

    double Lmax = LE(E);
    double l = L/Lmax;
    if(l>1+1e-6){
        return 0;
    }
    auto [Index,RectR] = H.Grid.FindElement(E,l);
    double Value = H[Index];
    if(method == 'f')
        return Value;
    double E0 = RectR.left();
    double E1 = RectR.right();

    double Lmax0 = LE(E0);
    double Lmax1 = LE(E1);;
    double l0 = RectR.inner().left();
    double l1 = RectR.inner().right();

    double L00 = l0*Lmax0;
    double L01 = l1*Lmax0;
    double L10 = l0*Lmax1;
    double L11 = l1*Lmax1;
    if(method == 'l')
        return Value/(0.5*(E1-E0)*(l1-l0));
    else
        return Value/(0.5*(E1-E0)*(L01+L11-L00-L10));
}

template <typename FuncType,typename HistoType>
auto ELH_toFunction(const HistoType &H,FuncType && LE,char method,size_t EN,size_t LN){
    double E0 = grob::Grid1Cast(H.Grid.first()).front();
    double E1 = grob::Grid1Cast(H.Grid.first()).back();
    double L0 = 0;
    double L1 = LE(E1);
    return grob::make_function_f(grob::mesh_grids(grob::GridUniform<double>( E0,E1,EN),
                                                  grob::GridUniform<double>( L0,L1,LN)),
                                 [&](double E,double L){
                                        return EL_dens(H,LE,method,E,L);
                                    });
}

int main(int argc,char **argv){
    boost::property_tree::ptree cmd_params;
    auto cmd_parse_log = parse_command_line_v1(argc,argv,cmd_params);
    if(!cmd_parse_log.empty()){
        std::cout << cmd_parse_log <<std::endl;
        return 1;
    }
    if(ptree_contain(cmd_params,"help")){
        print("convertions -histo [path_to_histo] -LE [path_to_LE_grid_function] -o output");
        print("or instread of -histo: -grid [path to grid] -values [path_to_values.bmat]");
        return 0;
    }
    //std::cout << cmd_params.get<std::string>(0) << std::endl;

    auto FPath = [&](const std::string & tree_path){
        return config_path_from(cmd_params.pgets(tree_path),cmd_params.pgets("config_path")).string();
    };


    bool is_txt = (cmd_params.get<std::string>("format","") == "text");
    typedef grob::MultiGridHisto<grob::GridVectorHisto<double>,std::vector<grob::GridVectorHisto<double>>> GType;
    typedef grob::Histogramm<GType,std::vector<double>> HType;

    const char * histo_path;
    std::vector<std::string> h_pathes = {"histo","H","Histo", "HISTO","elh",
                                         "ELH", "EL_histo","el_histo",
                                        "EL_Histo","EL_HISTO"};
    bool found_h_path = false;
    for(auto const& s : h_pathes){
        if(ptree_contain(cmd_params,s)){
            histo_path = s.c_str();
            found_h_path = true;
            break;
        }
    }
    stools::PtreeSerializator<boost::property_tree::ptree> S{};
    boost::property_tree::ptree P;
    boost::property_tree::ptree GrP;
    boost::property_tree::ptree VlP;
    if(found_h_path)
        boost::property_tree::read_json(cmd_params.pgets(histo_path),P);
    else{
        std::ifstream ifs_grid(cmd_params.pgets("grid"));
        if(!ifs_grid){
            print("no such file",cmd_params.pgets("grid"));
        }
        boost::property_tree::read_json(ifs_grid,GrP);
        //if(is_txt)
        //    boost::property_tree::read_json(cmd_params.pgets("values"),VlP);
    }
    HType Histo =  (found_h_path ? stools::DeSerialize<HType>(P,S) :
                                   HType(stools::DeSerialize<GType>(GrP,S),
                                         is_txt ? stools::DeSerialize<std::vector<double>>(VlP,S) :
                                                  std::get<0>(loadMatrix<double>(cmd_params.pgets("values"),
                                                                                 is_txt ? MatrixFormat::TEXT : MatrixFormat::BINARY))
                                                  ));

    const char * le_path = "";
    bool found_le_path = false;
    std::vector<std::string> le_ps = {"le","LE","l_e","L_E"};
    for(auto const& s : le_ps){
        if(ptree_contain(cmd_params,s)){
            le_path = s.c_str();
            found_le_path =true;
            break;
        }
    }
    if(!found_le_path){
        std::runtime_error("not found LE path");
    }
    boost::property_tree::ptree LEP;
    boost::property_tree::read_json(cmd_params.pgets(le_path),LEP);

    typedef grob::GridFunction<grob::linear_interpolator,
                        grob::GridUniform<double>,
                        std::vector<double>> LE_Type;
    auto LE_func = stools::DeSerialize<LE_Type> (LEP,S);


    const char * out_p = "";
    bool found_out_p = false;
    std::vector<std::string> out_ps = {"o","O","out","Out","OUT"};
    for(auto const& s : out_ps){
        if(ptree_contain(cmd_params,s)){
            out_p = s.c_str();
            found_out_p =true;
            break;
        }
    }
    if(!found_out_p){
        print("error: need -o flag");
        return 0;
    }

    size_t NE = (ptree_contain(cmd_params,"NE") ? cmd_params.get<int>("NE") : 100);
    size_t NL = (ptree_contain(cmd_params,"NL") ? cmd_params.get<int>("NL") : 100);

    char method = 'L';
    if(ptree_contain(cmd_params,"func") ||
       ptree_contain(cmd_params,"f") ||
       ptree_contain(cmd_params,"F"))
        method = 'f';
    else if (ptree_contain(cmd_params,"El") ||
             ptree_contain(cmd_params,"l")){
          method = 'l';
    }

    as_csv(ELH_toFunction(Histo,LE_func,method,NE,NL)).save(std::ofstream(cmd_params.pgets(out_p)),6,std::defaultfloat);

	return 0;
}
