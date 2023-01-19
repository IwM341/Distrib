#include <iostream>


#include "functions.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/arg_parser.hpp"
#include "func/dens_grid.hpp"
#include "func/matrix_functions.hpp"


int main(int argc,char**argv){
    boost::property_tree::ptree cmd_params;
    auto cmd_parse_log = parse_command_line_v1(argc,argv,cmd_params);
    if(!cmd_parse_log.empty()){
        std::cout << cmd_parse_log <<std::endl;
        return 0;
    }
    auto BM = BodyModel::fromFile(2e-3,cmd_params.pgets("body"));
    std::vector<Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator>>
            Funcs;
    std::vector<std::string> names;
    Function::UniformGrid<double> R(0,1,BM["Radius"].size());


    std::stringstream par_list(cmd_params.pgets("list"));
    bool is_relative = ptree_condition(cmd_params,"relative",false);
    std::string EL;
    while(par_list>>EL){
        if(std::find(EL.begin(),EL.end(),'\"') == EL.end() and !EL.empty()){
            if(ME.find(EL) != ME.end() and !is_relative)
                Funcs.push_back(decltype(Funcs)::value_type(R,BM[EL]*BM["RhoND"]));
            else
                Funcs.push_back(decltype(Funcs)::value_type(R,BM[EL]));
            names.push_back((EL));
        }
    }
    using namespace std::string_literals;

    auto gp_path = ptree_condition<std::string>(cmd_params,"gnuplot_path","gnuplot");

    Gnuplot gp(gp_path);
    if(ptree_condition(cmd_params,"logscale",false))
        gp.command("set logscale");
    for(size_t i=0;i<Funcs.size();++i){
        gp.plotd(Funcs[i].toString(),"with lines title \""s+names[i]+"\"");
    }
    gp.show();
    pipe_switcher ps;
    ps.add_pipe(&gp,&Gnuplot::command,"gp");
    ps.exec();
	return 0;
}
