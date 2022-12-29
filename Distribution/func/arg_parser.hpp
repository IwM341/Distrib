#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP
#include <boost/property_tree/json_parser.hpp>
#include <vector>
#include <string>
#include <regex>


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


std::string add_config_file(boost::property_tree::ptree &tree,const std::string_view &fname){
    if(tree.find("config") == tree.not_found()){
        tree.put("config",fname);
    }
    return "";
}

std::string parse_command_line(int argc,char ** argv,boost::property_tree::ptree &tree){
    using namespace std::string_literals;
    std::regex parametr("-+([a-z]\\w*)");
    std::regex parametr_seter("-+([a-z]\\w*)=(.*)");
    std::string last_par;
    std::smatch sm;

    std::string parse_log;

    for(size_t i=1;i<argc;++i){
        std::string arg = argv[i];
        if(std::regex_match(arg, sm, parametr)){
            arg = sm[1];
            if(!last_par.empty()){
                if(str_toupper(last_par) != "CONFIG"){
                    tree.put(last_par,true);
                }
                else{
                    parse_log += "expecting file after -config flag\n";
                    return parse_log;
                }

            }
            else
                last_par = arg;
        }
        else if(std::regex_match(arg, sm, parametr_seter)){
            std::string key = sm[1];
            std::string value = sm[2];
            if(str_toupper(key) != "CONFIG"){
                last_par = "";
                tree.put(key,value);
            }
            else{
                try{
                    std::ifstream cnf_file(value);
                    boost::property_tree::read_json(cnf_file,tree);
                }catch(std::exception & e){
                    parse_log += "cant parse config file "s + value + "\n";
                    return parse_log;
                }
            }
        }
        else{
            if(last_par == ""){
                parse_log += "got value"s + arg + " without flag";
                return parse_log;
            }
            if(str_toupper(last_par) != "CONFIG"){
                tree.put(last_par,arg);
            }else{
                try{
                    std::ifstream cnf_file(arg);
                    boost::property_tree::read_json(cnf_file,tree);
                }catch(std::exception & e){
                    parse_log += "cant parse config file "s + arg + "\n";
                    return parse_log;
                }
            }
            last_par.clear();
        }
    }
    return parse_log;
}



#endif//ARG_PARSER_HPP
