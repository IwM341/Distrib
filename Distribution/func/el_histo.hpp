#ifndef EL_HISTO_H
#define EL_HISTO_H
#include <utils>
#include <limits>

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
        double Lin = l_nd/Lmax(e_nd);
        /*
        if(e_nd <0 and Lin > 1.01){
            PVAR(e_nd);
            PVAR(Lin);
        }*/
        return HBase::putValue(V,e_nd,l_nd/Lmax(e_nd));
    }

    template <typename ofstream_type>
    void save(ofstream_type && stream,bool isBinary = false){
        if(isBinary)
            save_bin(stream);
        else
            save_text(stream);
    }
    template <typename ifstream_type>
    static EL_Histo load(ifstream_type &&ifs,bool isBinary = false){
        if(isBinary)
            return load_bin(ifs);
        else
            return load_text(ifs);
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

inline double e_intersect(double E0,double E1,double L0,double L1,double L){
    if(L0==L1){
        return 1.0/0.0;
    }
    else{
        return ((L-L1)*E0-E1*(L-L0))/(L0-L1);
    }
}
inline double e_intersect(double Eav,double Lav,double slope,double L){
    if(slope == 0){
        return 1.0/0.0;
    }
    else{
        return Eav+(L-Lav)/slope;
    }
}

inline double l_e(double E0,double E1,double L0,double L1,double E){
    return L0  + (L1-L0)/(E1-E0)*(E-E0);
}
inline double l_e(double E0,double L0,double scope,double E){
    return L0  + scope*(E-E0);
}
struct EL_column{
    double E;
    double L0;
    double L1;
    EL_column(double E,double L0,double L1):E(E),L0(L0),L1(L1){}
    EL_column():L1(0),L0(-1){}
    inline bool valid()const{return L0<=L1;}
    friend std::ostream & operator << (std::ostream & os,const EL_column&R){
        os << "(" << R.E << " [" << R.L0 << ", " << R.L1  << "])";
        return os;
    }
};


struct EL_rect{
    EL_column C0;
    EL_column C1;
    EL_rect(double E0,double E1,double l00,double l01,double l10,double l11):
        C0(E0,l00,l01),C1(E1,l10,l11){}
    EL_rect():C0(0,0,0),C1(-1,0,0){}
    inline double & E0(){return C0.E;}
    inline const double & E0()const{return C0.E;}
    inline double & E1(){return C1.E;}
    inline const double & E1()const{return C1.E;}

    inline double & L00(){return C0.L0;}
    inline const double & L00()const{return C0.L0;}
    inline double & L01(){return C0.L1;}
    inline const double & L01()const{return C0.L1;}

    inline double & L10(){return C1.L0;}
    inline const double & L10()const{return C1.L0;}
    inline double & L11(){return C1.L1;}
    inline const double & L11()const{return C1.L1;}

    bool valid() const {
        return (E1() >= E0()) and (L01() >= L00()) and (L11()>=L10());
    }
    inline std::vector<EL_column> SquareIntersection(double E0,double E1,double L0,double L1)const{
        double slope0 = (C1.L0-C0.L0)/(C1.E-C0.E);
        double slope1 = (C1.L1-C0.L1)/(C1.E-C0.E);
        double Eav = 0.5*(C0.E+C1.E);
        double L0av = 0.5*(C1.L0+C0.L0);
        double L1av = 0.5*(C1.L1+C0.L1);

        //std::array<double,6> crit_e;
        std::vector<EL_column> Out_Columns;

        double oE0 = C0.E;
        double lmin_0 = std::max(L00(),L0);
        double lmax_0 = std::min(L01(),L1);
        if(oE0 > E1){
            return Out_Columns;
        }
        if(oE0 <E0){
            oE0 = E0;
            lmin_0 = std::max(l_e(Eav,L0av,slope0,oE0),L0);
            lmax_0 = std::min(l_e(Eav,L1av,slope1,oE0),L1);
        }
        if(lmin_0 <= lmax_0){
            Out_Columns.push_back(EL_column(oE0,lmin_0,lmax_0));
        }

        double oE1= C1.E;
        double lmin_1 = std::max(L10(),L0);
        double lmax_1 = std::min(L11(),L1);


        if(oE1 < oE0){
            return Out_Columns;
        }
        if(oE1 > E1){
            oE1 = E1;
            lmin_1 = std::max(l_e(Eav,L0av,slope0,oE1),L0);
            lmax_1 = std::min(l_e(Eav,L1av,slope1,oE1),L1);
        }
        if(lmin_1 <= lmax_1){
            Out_Columns.push_back(EL_column(oE1,lmin_1,lmax_1));
        }


        double e00 = e_intersect(Eav,L0av,slope0,L0);
        double lmin00 = L0;
        double lmax00 = std::min(l_e(Eav,L1av,slope1,e00),L1);
        if(e00 >= oE0 and e00 <= oE1 and lmax00>=lmin00){
            Out_Columns.push_back(EL_column(e00,lmin00,lmax00));
        }

        double e01 = e_intersect(Eav,L0av,slope0,L1);
        if(e01 >= oE0 and e01 <= oE1){
            Out_Columns.push_back(EL_column(e01,L1,L1));
        }
        //lmax01 = L1
        //lmax01 = L1
        double e10 = e_intersect(Eav,L1av,slope1,L0);
        if(e10 >= oE0 and e10 <= oE1){
            Out_Columns.push_back(EL_column(e10,L0,L0));
        }
        //lmin10 = L0
        //lmax10 = L0
        double e11 = e_intersect(Eav,L1av,slope1,L1);
        double lmax11 = L1;
        double lmin11 = std::max(l_e(Eav,L0av,slope0,e11),L0);
        //lmax10 = L1
        if(e11 >= oE0 and e11 <= oE1 and lmax11>=lmin11){
            Out_Columns.push_back(EL_column(e11,lmin11,lmax11));
        }

        std::sort(Out_Columns.begin(),Out_Columns.end(),
                  [](const EL_column & C1,const EL_column & C2){return C1.E < C2.E;});
        return Out_Columns;

    }
    friend std::ostream & operator << (std::ostream & os,const EL_rect&R){
        os << "[" << R.E0() << ", " << R.E1() << "]" << std::endl;
        os << "[<" << R.L00() << " | " << R.L01() << ">, <" <<  R.L10() << " | " << R.L11() << ">]";
        return os;
    }
};

template <typename Measure,typename Container>
inline double EL_CascadeIntegration(const Measure&m,const Container&Arr) noexcept{
    if(Arr.size() < 2)
        return 0.0;
    double sum = 0;
    for(size_t i=0;i<Arr.size()-1;++i)
        sum += m(Arr[i],Arr[i+1]);
    return sum;
}



inline double EL_rect_measure(const EL_rect &R) noexcept{
    return (R.C0.L1+R.C1.L1-R.C0.L0-R.C1.L0)*(R.C1.E-R.C0.E)/2;
}
inline double EL_column_measure(const EL_column &C0,const EL_column &C1) noexcept{
    return (C0.L1+C1.L1-C0.L0-C1.L0)*(C1.E-C0.E)/2;
}


inline double EL_reg(double x) noexcept{
    return sqrt(-x*x*x);
}
double EL_bin_volume(double E0,double E1,double L00,double L01,double L10,double L11,
                  double T00,double T01,double T10,double T11) noexcept{
    if(std::isinf(T00) or std::isinf(T01)
            or std::isinf(T10) or std::isinf(T11)
            or E1 >=0){
        return std::numeric_limits<double>::infinity();
    }
    double F00 = T00*EL_reg(E0);
    double F01 = T01*EL_reg(E0);
    double F10 = T10*EL_reg(E1);
    double F11 = T11*EL_reg(E1);
    double sE0 = sqrt(-E0);
    double sE1 = sqrt(-E1);

    double Fa0 = (F00+F10)/2;
    double Fd0 = (F10-F00);
    double Fa1 = (F01+F11)/2;
    double Fd1 = (F11-F01);
    double La0 = (L00+L10)/2;
    double Ld0 = (L10-L00);
    double La1 = (L01+L11)/2;
    double Ld1 = (L11-L01);


    double int_1 = 2*(E1-E0)/(sE0*sE1*(sE0+sE1));
    //integratin of 1
    double fac_3 = sE0-sE1;
    double int_2 = fac_3*fac_3*fac_3/(sE0*sE1*(E1-E0));
    //integration of x-1/2
    //where x = (E-E0)/(E1-E0)
    double int_1_part = ((2*Fa0+Fa1)*La0+La1*(Fa0+2*Fa1))*(La1-La0)/3;

    double int_2_part = ( (La1-La0)*(Fd1*(2*La1+La0)+Fd0*(2*La0+La1)) +
                          (Fa0-Fa1)*(La1-La0)*(Ld0-Ld1)+(3*(Fa0+Fa1))*(-La0*Ld0+La1*Ld1))/3;
    return (int_1*int_1_part+int_2_part*int_2);
}

template <typename Phitype,typename ...Args>
double ELH_density(const EL_Histo<Args...> &H,const Phitype&phi,double E,double L){
    if(E < H.Grid._a() or E > H.Grid._b()){
        return 0;
    }
    double Lmax = H.Lmax(E);
    if(Lmax == .0){
        return .0;
    }
    double l = L/Lmax;
    if(l>1){
        return 0;
    }
    size_t i = H.Grid.pos(E);
    size_t j = H.values[i].Grid.pos(l);
    double E0 = H.Grid[i];
    double E1 = H.Grid[i+1];
    double Value = H.values[i].values[j];

    double Lmax0 = H.Lmax[i];
    double Lmax1 = H.Lmax[i+1];
    double l0 = H.values[i].Grid[j];
    double l1 = H.values[i].Grid[j+1];

    double L00 = l0*Lmax0;
    double L01 = l1*Lmax0;
    double L10 = l0*Lmax1;
    double L11 = l1*Lmax1;

    auto TI = CalculateTrajectory(phi,E0,L00,100);
    double T00 = TI.T_in+TI.T_out;
    TI = CalculateTrajectory(phi,E0,L01,100);
    double T01 = TI.T_in+TI.T_out;
    TI = CalculateTrajectory(phi,E1,L10,100);
    double T10 = TI.T_in+TI.T_out;
    TI = CalculateTrajectory(phi,E1,L11,100);
    double T11 = TI.T_in+TI.T_out;
    return Value/EL_bin_volume(E0,E1,L00,L01,L10,L11,T00,T01,T10,T11);
}

#endif
