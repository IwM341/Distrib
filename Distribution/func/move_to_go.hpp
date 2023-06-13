#ifndef MOVE_TO_GO_HPP
#define MOVE_TO_GO_HPP

#include <utils>
#include "grob/grid.hpp"
#include "grob/grid_objects.hpp"
#include "grob/object_serialization.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
template <typename T>
grob::GridUniform<T> grobGridPointCast(Function::UniformGrid<T> const& Grid){
    return grob::GridUniform<T>(Grid._a(),Grid._b(),Grid.size());
}


template <typename T>
grob::GridUniformHisto<T> grobGridHistoCast(Function::UniformGrid<T> const& Grid){
    return grob::GridUniformHisto<T>(Grid._a(),Grid._b(),Grid.size());
}

template <typename T>
grob::GridVector<T> grobGridPointCast(Function::VectorGrid<T> const& Grid){
    return grob::GridVector<T>(Grid.begin(),Grid.end());
}
template <typename T>
grob::GridVectorHisto<T> grobGridHistoCast(Function::VectorGrid<T> const& Grid){
    return grob::GridVectorHisto<T>(Grid.begin(),Grid.end());
}



template <typename T,typename GridType,typename Interpol>
auto grobGetPointGrid(Function::GridFunction<T,GridType,Interpol> const&F){
    return grobGridPointCast(F.Grid);
}
template <typename T,typename GridType,typename Interpol,typename GridType1,typename...Other>
auto grobGetPointGrid(Function::GridFunction<T,GridType,Interpol,GridType1,Other...> const&F){
    return grob::make_grid_f(grobGridPointCast(F.Grid),[&F](size_t i){
        return grobGetPointGrid(F.values[i]);
    });
}


template <typename T,typename...Args>
auto grobFunctionCast(Function::GridFunction<T,Args...> const&F){
    return grob::make_function(grobGetPointGrid(F),F.AllValues());
}

template <typename T,typename GridType>
auto grobGetHistoGrid(Function::Histogramm<T,GridType> const&H){
    return grobGridHistoCast(H.Grid);
}

template <typename T,typename GridType,typename GridType1,typename...GridTypes>
auto grobGetHistoGrid(Function::Histogramm<T,GridType,GridType1,GridTypes...> const&H){
    return grob::make_grid_f(grobGridHistoCast(H.Grid),
                                   [&H](size_t i){return grobGetHistoGrid(H.values[i]);}
                             );
}

template <typename T,typename...Args>
auto grobHistoCast(Function::Histogramm<T,Args...> const&H){
    return grob::make_function(grobGetHiastoGrid(H),H.AllValues());
}

template <typename T,typename Stream>
void WriteObject(T const & go_item,Stream && os){
    stools::PtreeSerializator<boost::property_tree::ptree> S;
    auto P = stools::Serialize(go_item,S);
    boost::property_tree::write_json(os,P);
}

template <typename T,typename Stream>
auto ReadObject(T const & go_item,Stream && os){
    stools::PtreeSerializator<boost::property_tree::ptree> S;
    boost::property_tree::ptree P;
    boost::property_tree::read_json(os,P);
    return stools::DeSerialize<T>(P,S);

}



#endif//MOVE_TO_GO_HPP
