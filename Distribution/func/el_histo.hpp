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
    template <typename LambdaLmaxType>
    EL_Histo(const EGridType & Grid,const LambdaLmaxType & FLmax):
        Lmax(Grid,Vector(Grid.size(),[&Grid,FLmax](size_t i){
           return   FLmax(Grid[i]);
        })),HBase(Grid){}

    bool putValue(T V,double e_nd,double l_nd){
        return HBase::putValue(V,e_nd,l_nd/Lmax(e_nd));
    }

};


#endif
