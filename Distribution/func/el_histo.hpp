#ifndef EL_HISTO_H
#define EL_HISTO_H
#include <utils>
template <typename T,typename EGridType,typename LGridType>
class EL_Histo:
        public Function::Histogramm<T,EGridType,LGridType>
{
public:
    typedef Function::Histogramm<T,EGridType,LGridType> HBase;
    Function::GridFunction<T,EGridType,Function::LinearExtrapolator> Lmax;
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

    template <typename ofstream_type>
    void save_text(ofstream_type && stream){
        for(size_t i=0;i<Lmax.Grid.size();++i){
            stream << Lmax.Grid[i] << "\t";
        }
        stream <<std::endl;
        for(size_t i=0;i<Lmax.Grid.size();++i){
            stream << Lmax.values[i] << "\t";
        }
        stream <<std::endl;
        for(size_t i=0;i<Lmax.Grid.size()-1;++i){
            for(size_t j=0;j<HBase::values[i].Grid.size();++j){
                stream << HBase::values[i].Grid[j] << "\t";
            }
            stream <<std::endl;
        }
    }

    template <typename ifstream_type>
    static EL_Histo load_text(ifstream_type && ifs){
        EL_Histo ELH;
        double tmp_value;




        std::string line;
        std::list<double> val_list;

        val_list.clear();
        std::getline(ifs,line);
        std::istringstream S(line);

        while(S>>tmp_value)
            val_list.push_back(tmp_value);

        std::vector<double> Egrid(val_list.begin(),val_list.end());


        val_list.clear();
        std::getline(ifs,line);
        std::istringstream S1(line);
        while(S1>>tmp_value)
            val_list.push_back(tmp_value);
        std::vector<double> Lmax_vals(val_list.begin(),val_list.end());

        (HBase&)ELH= HBase(Function::VectorGrid(Egrid));
        ELH.Lmax = decltype(ELH.Lmax)(Function::VectorGrid(Egrid),Lmax_vals);

        for(size_t i=0;i<ELH.Lmax.Grid.size()-1;++i){
            val_list.clear();
            std::getline(ifs,line);
            std::istringstream S2(line);

            while(S2>>tmp_value)
                val_list.push_back(tmp_value);

            ELH.values[i] = typename std::remove_reference<decltype(ELH.values[i])>::type(
                        Function::VectorGrid(std::vector<T>(val_list.begin(),val_list.end()))
                        );
        }
        return ELH;
    }

    template <typename ofstream_type>
    void save_bin(ofstream_type &&ofs ){

        size_t Ne = HBase::Grid.size();
        ofs.write((char*)&Ne,sizeof(Ne));
        for(size_t i=0;i<Lmax.Grid.size();++i){
            ofs.write((char*)&Lmax.Grid[i],sizeof(double));
        }
        ofs.write((char*)Lmax.values.data(),Lmax.values.size()*sizeof(double));
        for(size_t i=0;i<Lmax.Grid.size()-1;++i){
            size_t Nl = HBase::values[i].Grid.size();
            ofs.write((char*)&Nl,sizeof(Nl));
            for(size_t j=0;j<HBase::values[i].Grid.size();++j){
                ofs.write((char*)&HBase::values[i].Grid[j],sizeof(double));
            }
        }
    }
    template <typename ifstream_type>
    static EL_Histo load_bin(ifstream_type &&ifs){
        EL_Histo ELH;
        T tmp_value;

        size_t Ne;
        ifs.read((char*)&Ne,sizeof(Ne));
        std::vector<double> EG(Ne);
        std::vector<double> LG(Ne);
        ifs.read((char*)EG.data(),EG.size()*sizeof(double));
        ifs.read((char*)LG.data(),LG.size()*sizeof(double));
        ELH.Lmax = decltype(ELH.Lmax)(Function::VectorGrid<double>(EG),LG);
        (HBase&)ELH= HBase(Function::VectorGrid<double>(EG));
        for(size_t i=0;i<Ne-1;++i){
            size_t Nl;
            ifs.read((char*)&Nl,sizeof(Nl));
            std::vector<double> Lg(Nl);
            ifs.read((char*)Lg.data(),Lg.size()*sizeof(double));
            ELH.values[i] = typename std::remove_reference<decltype(ELH.values[i])>::type
                    (Function::VectorGrid<double>(Lg));
        }
        return ELH;
    }
};


template <typename V,typename GridType,typename...Args>
auto MergeHisto(const Function::Histogramm<V,GridType,Args...> &H){
    return Function::Histogramm<V,GridType>(H.Grid,vmap([](auto x){return x.summ();},H.values));
}

#endif
