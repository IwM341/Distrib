#include <iostream>
#include "func/arg_parser.hpp"
#include <sstream>
#include <boost/filesystem.hpp>
#include <utils>

int main(int argc, char**argv){
    boost::property_tree::ptree par;

    auto log = parse_command_line(argc,argv,par);
    if(!log.empty()){
        std::cout << log <<std::endl;
    }

    std::stringstream myJsonEncodedData;
    boost::property_tree::write_json(myJsonEncodedData, par);
    std::cout << myJsonEncodedData.str() <<std::endl;

    auto val = par.find("lh_in") ;
    if(val != par.not_found()){
        std::cout << val->second.data()<<std::endl;
    }

    boost::filesystem::path p1("D:\\Tmp\\probs\\c_test");
    boost::filesystem::path p2("D:/Tmp/probs/c_test");
    boost::filesystem::path p3("D:/Tmp/probs/../probs/c_test");
    boost::filesystem::path p4("C:/Windows/System");

    auto ex_p = {p1,p2,p3,p4};

    for(const auto & p : ex_p)
        std::cout <<p.lexically_normal()<<std::endl;
    std::cout <<std::endl;

    std::cout << boost::filesystem::relative(p4,p1).lexically_normal()<<std::endl;

    for(const auto & p : ex_p)
        std::cout <<p.lexically_relative(boost::filesystem::current_path())<<std::endl;
    std::cout <<std::endl;

    for(const auto & p : ex_p)
        std::cout <<p.lexically_normal().lexically_relative(boost::filesystem::current_path())<<std::endl;
    std::cout <<std::endl;

    std::cout <<boost::filesystem::absolute(p3,p4).lexically_normal()<<std::endl;
    std::cout <<std::endl;

    PVAR(p1.is_absolute());
    PVAR(p4.is_absolute());
    PVAR(p1.lexically_relative(boost::filesystem::current_path()).is_absolute());

    return 0;
}
