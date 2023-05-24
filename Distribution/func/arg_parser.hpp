#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include <regex>
#include <stack>

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

boost::filesystem::path config_path_from(const boost::filesystem::path & filename,
                            const boost::filesystem::path & config){
    if(filename.lexically_normal().is_absolute()){
        return filename.lexically_normal();
    }
    else{
        return boost::filesystem::absolute(filename,config).lexically_normal();
    }
}

/*!
 * \brief gives relative tp config path of filename path
 */
boost::filesystem::path config_path_to(const boost::filesystem::path & filename,
                            const boost::filesystem::path & config){
    auto sr = boost::filesystem::relative(filename,config);
    if(sr.empty())
        return boost::filesystem::absolute(filename.lexically_normal());
    else
        return sr.lexically_normal();
}


#define pgets get<std::string>
void merge_ptree(boost::property_tree::ptree &tree,const boost::property_tree::ptree &tree1){
    for(const auto &item:tree1)
        tree.push_back(item);
}
inline bool ptree_contain(const boost::property_tree::ptree &tree,const std::string & element){
    return tree.find(element) != tree.not_found();
}

inline std::string ptree_gets(const boost::property_tree::ptree &tree,const std::string & element){
    try {
        return tree.get<std::string>(element);
    } catch (...) {
        return "";
    }
}

template <typename T>
T ptree_condition(const boost::property_tree::ptree &tree,
                            const std::string & element,
                            const T & default_option){
    if(ptree_contain(tree,element))
        return tree.get<T>(element);
    else
        return default_option;
}

std::string add_config_file(boost::property_tree::ptree &tree,const std::string &fname){
    using namespace std::string_literals;
    std::ifstream cnf_file(fname);
    if(!cnf_file.is_open()){
        return "cant open config file "s + fname+ "\n";
    }
    try{
        boost::property_tree::ptree ftree;
        boost::property_tree::read_json(cnf_file,ftree);
        merge_ptree(tree,ftree);
    }catch(std::exception & e){
        return  "cant parse config file "s + fname+ e.what() + "\n";
    }
    tree.put("config_path",boost::filesystem::path(fname).parent_path().string());
    return "";
}

std::string parse_command_line(int argc,char ** argv,boost::property_tree::ptree &tree){
    using namespace std::string_literals;
    std::regex parametr("-+([a-z]\\w*)");
    std::regex parametr_seter("-+([a-z]\\w*)=?(.*)");
    std::string last_par;
    std::smatch sm;
    tree.put("config_path","");
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
            last_par = arg;
        }
        else if(std::regex_match(arg, sm, parametr_seter)){
            std::string key = sm[1];
            std::string value = sm[2];
            if(!last_par.empty()){
                if(str_toupper(last_par) != "CONFIG"){
                    tree.put(last_par,true);
                }
                else{
                    parse_log += "expecting file after -config flag\n";
                    return parse_log;
                }
            }
            if(str_toupper(key) != "CONFIG"){
                tree.put(key,value);
            }
            else{
                parse_log += add_config_file(tree,value);
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
                parse_log += add_config_file(tree,arg);
            }
            last_par.clear();
        }
    }
    if(!last_par.empty()){
        tree.put(last_par,true);
    }
    return parse_log;
}
std::string parse_command_line_v1(int argc,char ** argv,boost::property_tree::ptree &tree){
    using namespace std::string_literals;
    std::regex parametr("-+([a-zA-Z]\\w*)");
    std::regex parametr_seter("-+([a-zA-Z]\\w*)=?(.*)");

    std::smatch sm;

    enum param_type{PARAM,VALUE,SETTER};

    std::stack<std::tuple<param_type,std::string,std::string>> st_params;

    tree.put("config_path","./");
    std::string parse_log;

    for(size_t i=1;i<argc;++i){
        std::string arg = argv[i];
        if(std::regex_match(arg, sm, parametr)){
            arg = sm[1];
            st_params.push(std::make_tuple(PARAM,arg,""));
        }
        else if(std::regex_match(arg, sm, parametr_seter)){
            std::string key = sm[1];
            std::string value = sm[2];
            st_params.push(std::make_tuple(SETTER,key,value));
        }
        else{
            st_params.push(std::make_tuple(VALUE,"",arg));
        }
    }
    while(!st_params.empty()){
        std::string key,val;
        auto [type,key_p,val_p] = st_params.top();
        st_params.pop();
        if(type == SETTER){
            key = key_p;
            val = val_p;
        }
        else if(type == VALUE){
            val = val_p;
            if(st_params.empty()){
                parse_log += "unexpected argument: "s + val + "\n";
                break;
            }
            auto [type1,key_p1,val_p1] = st_params.top();
            st_params.pop();
            if(type1 != PARAM){
                parse_log += "unexpected argument: "s + val + "\n";
                break;
            }
            key = key_p1;
        }
        else{
            if(str_toupper(key_p) == "CONFIG"){
                parse_log += "-config flag requires filename";
                break;
            }
            key = key_p;
            val = "true";
        }
        if(str_toupper(key) != "CONFIG"){
            tree.put(key,val);
        }
        else{
            parse_log += add_config_file(tree,val);
        }
    }
    return parse_log;
}




#endif//ARG_PARSER_HPP
