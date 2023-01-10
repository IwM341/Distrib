#include <iostream>
#include "func/arg_parser.hpp"
#include <sstream>
#include <boost/filesystem.hpp>
#include <utils>
#include <sstream>
int main(int argc, char**argv){
    boost::property_tree::ptree par;

    auto log = parse_command_line(argc,argv,par);
    std::cout << log <<std::endl;

    for(const auto & [key,value]: par){
        std::cout << key << " : " << value.get_value<std::string>() << std::endl;
    }
    std::cout << (par.find("config") != par.not_found())<<std::endl;
    std::cout << "config path" << par.get_value<std::string>("config") <<std::endl;

    boost::filesystem::path file(par.get<std::string>("myfile"));
    boost::filesystem::path cnfpath(par.get<std::string>("config"));
    std::cout << config_path_from(file,cnfpath) <<std::endl;
    std::ifstream ifs(config_path_from(file,cnfpath));
    std::string s;
    ifs>>s;
    PVAR(s);

    PVAR(par.get<std::string>("config1"));
    PVAR(par.get<std::string>("config"));

    return 0;
}
