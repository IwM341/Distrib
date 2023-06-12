#ifndef LOAD_HISTO_HPP
#define LOAD_HISTO_HPP

#include <cmath>
#include <fstream>

#include <time.h>
#include "el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "arg_parser.hpp"
#include "dens_grid.hpp"
#include "matrix_functions.hpp"


template <typename grid_type,typename value_type>
auto loadHisto(std::string const & grid_fname,
               std::string const & values_fname,
               boost::filesystem::path const &programm_path){
    stools::PtreeSerializator<boost::property_tree::ptree> S{};
    std::ifstream grid_file(
                config_path_from(grid_fname,
                                 programm_path).string()
                );
    if(!grid_file.is_open()){
        throw std::runtime_error("can't open histo filename " +
                                 config_path_from(grid_fname,programm_path).string());
    }
    std::ifstream values_file(
                config_path_from(values_fname,
                                 programm_path).string(),
                    std::ios::binary
                );
    if(!values_file.is_open()){
        throw std::runtime_error("can't open histo filename " +
                                 config_path_from(values_fname,programm_path).string());
    }
    boost::property_tree::ptree grid_fn;
    boost::property_tree::read_json(grid_file,grid_fn);
    return grob::make_histo(stools::DeSerialize<grid_type>(
                                        grid_fn
                                        ,S),
                            std::get<0>(LoadMatrixBinary<value_type>(values_file))
                            );
}

template <typename grid_type,typename value_type>
auto loadHisto(boost::property_tree::ptree const &h_info,
               boost::filesystem::path const &programm_path){
    //typedef grob::MultiGridHisto<grob::GridVectorHisto<double>,
    //            std::vector<grob::GridVectorHisto<double>>
    //        > grid_type;
    stools::PtreeSerializator<boost::property_tree::ptree> S{};
    if(ptree_contain(h_info,"grid") && ptree_contain(h_info,"values")){
        std::ifstream grid_file(
                    config_path_from(h_info.get<std::string>("grid"),
                                     programm_path).string()
                    );
        std::ifstream values_file(
                    config_path_from(h_info.get<std::string>("values"),
                                     programm_path).string(),
                        std::ios::binary
                    );
        boost::property_tree::ptree grid_fn;
        boost::property_tree::read_json(grid_file,grid_fn);
        return grob::make_histo(stools::DeSerialize<grid_type>(
                                            grid_fn
                                            ,S),
                                std::get<0>(LoadMatrixBinary<value_type>(values_file))
                                );
    }else if(ptree_contain(h_info,"histo")){
        boost::property_tree::ptree histo_fn;
        std::ifstream histo_file(
                    config_path_from(h_info.get<std::string>("histo"),
                                     programm_path).string()
                    );
        boost::property_tree::read_json(histo_file,histo_fn);
        return stools::DeSerialize<grob::Histogramm<grid_type,
                std::vector<value_type>>>(histo_fn,S);
    }
    else{
        throw std::runtime_error("expect histo info");
    }
}




#endif//LOAD_HISTO_HPP
