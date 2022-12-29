#ifndef EL_HISTO_H
#define EL_HISTO_H
#include <utils>
template <typename T,typename EGridType,typename LGridType>
class EL_Histo:
        public Function::Histogramm<T,EGridType,LGridType>
{
public:
    typedef Function::Histogramm<T,EGridType,LGridType> HBase;
    Function::GridFunction<double,EGridType,Function::LinearExtrapolator> Lmax;
    EL_Histo(){}
    EL_Histo(size_t){}
    template <typename LambdaLmaxType>
    EL_Histo(const EGridType & Grid,const LambdaLmaxType & FLmax):
        Lmax(Grid,Vector(Grid.size(),[&Grid,FLmax](size_t i){
           return   FLmax(Grid[i]);
        })),HBase(Grid){}

    template <typename...HArgs,typename LambdaLmaxType>
    EL_Histo(const Function::Histogramm<HArgs...>& Hist,const LambdaLmaxType & FLmax):
        Lmax(Hist.Grid,Vector(Hist.Grid.size(),[&Hist,FLmax](size_t i){
           return   FLmax(Hist.Grid[i]);
        })),HBase(Hist){}

    bool putValue(T V,double e_nd,double l_nd){
        return HBase::putValue(V,e_nd,l_nd/Lmax(e_nd));
    }

};


template <typename V,typename GridType,typename...Args>
auto MergeHisto(const Function::Histogramm<V,GridType,Args...> &H){
    return Function::Histogramm<V,GridType>(H.Grid,vmap([](auto x){return x.summ();},H.values));
}

#endif
