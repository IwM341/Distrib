#include "functions.hpp"
#include <cmath>
#include <fstream>
#include "complex/complex_ex.hpp"
#include <time.h>
const double mp = 0.938;
const double me = 0.51e-3;
auto PhiFactor55p = [](double mk,double q)->double{
    return fast_pow(q*q/(mp*mk),2);
};
auto PhiFactor55e = [](double mk,double q)->double{
    return fast_pow(q*q/(me*mk),2);
};
auto PhiFactorSS = [](double mk,double q)->double{
    return 1.0;
};

#if defined(MASS)
int main(void){


    const auto& BM = BodyModel(2.03e-4,"D:/Important/articles/DMFramework/celstial_models/jupiter_model.dat");


    std::vector<double> Vmk = {0.5,1,2,10};
    //std::vector<double> Vmk = {0.4};
    size_t N = Vmk.size();

    std::vector<std::pair<double,double>> Fel(N);
    std::vector<std::pair<double,double>> FelE(N);
    std::vector<std::pair<double,double>> Fione(N);
    std::vector<std::pair<double,double>> Fion(N);
    std::vector<std::pair<double,double>> Fmgd(N);

    auto t_start = clock();
    #pragma omp parallel for
    for(size_t i=0;i<N;++i){
        Fel[i] = SupressFactor(Vmk[i],0,ELASTIC,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
        FelE[i] = SupressFactor(Vmk[i],0,ELASTIC,ELECTRON,PhiFactor55e,BM,"H",1000000,U0*1.1,U0);
        Fion[i] = SupressFactor(Vmk[i],0,IONIZATION,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
        Fione[i] = SupressFactor(Vmk[i],0,IONIZATION,ELECTRON,PhiFactor55e,BM,"H",1000000,U0*1.1,U0);
        Fmgd[i] = SupressFactor(Vmk[i],0,MIGDAL,PROTON,PhiFactor55p,BM,"H",1000000,U0*1.1,U0);
    }
    std::cout << "time = " << (clock()- t_start) <<std::endl;
    std::ofstream ofsEl("el55.dat");
    std::ofstream ofsIon("ion55.dat");
    std::ofstream ofsIonE("ione55.dat");
    std::ofstream ofsElE("elE55.dat");
    std::ofstream ofsMgd("mgd55.dat");

    ofsEl << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fel) << std::endl;
    ofsIon << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fion) << std::endl;
    ofsIonE << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fione) << std::endl;
    ofsElE << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,FelE) << std::endl;

    ofsMgd << "mk\tF\tdF\n"<< Function::GridFunction1(Vmk,Fmgd) << std::endl;
    return 0;

}
#elif defined(DELTAMASS)
int main(void){


    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    const double mk = 10.0;

    std::vector<double> Vdmk = Vector(10,[](size_t i){return i*1e-4;});
    //PVAR(Vdmk);
    //PVAR(BM["Vesc"]);
    //std::vector<double> Vmk = {0.4};
    size_t N = Vdmk.size();

    std::vector<std::pair<double,double>> Fel(N);
    std::vector<std::pair<double,double>> FelE(N);
    std::vector<std::pair<double,double>> Fione(N);
    std::vector<std::pair<double,double>> Fion(N);
    std::vector<std::pair<double,double>> Fmgd(N);

    auto t_start = clock();
    #pragma omp parallel for
    for(size_t i=0;i<N;++i){
        Fel[i] = SupressFactor(mk,Vdmk[i],ELASTIC,PROTON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        FelE[i] = SupressFactor(mk,Vdmk[i],ELASTIC,ELECTRON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        Fion[i] = SupressFactor(mk,Vdmk[i],IONIZATION,PROTON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        Fione[i] = SupressFactor(mk,Vdmk[i],IONIZATION,ELECTRON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
        Fmgd[i] = SupressFactor(mk,Vdmk[i],MIGDAL,PROTON,PhiFactorSS,BM,"H",1000000,U0*1.1,U0);
    }
    std::cout << "time = " << (clock()- t_start) <<std::endl;
    std::ofstream ofsEl("el55dm.dat");
    std::ofstream ofsIon("ion55dm.dat");
    std::ofstream ofsIonE("ione55dm.dat");
    std::ofstream ofsElE("elE55dm.dat");
    std::ofstream ofsMgd("mgd55dm.dat");

    ofsEl << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fel) << std::endl;
    ofsIon << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fion) << std::endl;
    ofsIonE << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fione) << std::endl;
    ofsElE << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,FelE) << std::endl;

    ofsMgd << "d_mk\tF\tdF\n"<< Function::GridFunction1(Vdmk,Fmgd) << std::endl;
    return 0;

}
#elif defined(TEST_BORDERS)
//TEST SUCCEED
//https://colab.research.google.com/drive/1qUEolPBBhXua3gXKSBi_QDfXSrb03oyT#scrollTo=4ld9zqreuCTh
int main(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");

    double Emin = -BM["phi"][0];
    double E = -0.1;

    //auto V = BM["Radius"]*vmap([](double x){return sqrt(2*x);},E + BM["phi"]);
    //PVAR(BM["phi"]);
    //double Lrm = *std::max_element(V.begin(),V.end());

    auto L_simple = [&BM](double _E){
        auto V = BM["Radius"]*vmap(sqrt,2*_E + 2*BM["phi"]);
        return *std::max_element(V.begin(),V.end());
    };

    //PVAR(maxLnd(BM,E));
    //PVAR(maxLnd(BM,-BM["phi"][2]));
    //PVAR(Lrm);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator>
            F(Function::UniformGrid<double>(Emin,0,101));
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator>
            F_t(Function::UniformGrid<double>(Emin,0,101));
    F.map([&BM,E](double _E){return maxLnd(BM,_E);});
    F_t.map(L_simple);


    std::ofstream F_file("D:\\Tmp\\Downloads\\LE_func.dat");
    std::ofstream F_true_file("D:\\Tmp\\Downloads\\LE_true_func.dat");
    F_file << F.toString() <<std::endl;
    F_true_file << F_t.toString() <<std::endl;

    return 0;
}
#elif defined(TEST_DISTRIBUTIONS)
//TEST SUCCEED
//https://colab.research.google.com/drive/12ahUMNi_Isr9gfgwIeDUMncZzqwP_6Rv?usp=sharing
int main(void){
    size_t N = 1000000;
    srand(time(0));
    auto G = [](){return 1.0-rand()/(RAND_MAX+1.0);};

    std::vector<MC::MCResult<double>> A1(N),A2(N);
    double sm1 = 0;
    for(size_t i=0;i<N;++i){
        auto res = Velocity(G,0,1,1);
        A1[i] = MC::MCResult<double>(res.Result.norm(),res.RemainDensity);
        sm1+= res.RemainDensity/N;
    }
    PVAR(sm1);

    double sm2 = 0;
    for(size_t i=0;i<N;++i){
        auto res = VelocityConstrained(G,0,0.5,2.5,1,1);
        A2[i] = MC::MCResult<double>(res.Result.norm(),res.RemainDensity);
        sm2 += res.RemainDensity/N;
    }
    PVAR(sm2);
    std::ofstream f1("distrib_simple.txt");
    std::ofstream f2("distrib_sophist.txt");

    for(size_t i=0;i<N;++i){
        f1 << A1[i].Result << "\t" << A1[i].RemainDensity << std::endl;
        f2 << A2[i].Result << "\t" << A2[i].RemainDensity << std::endl;
    }
}

#elif defined(TEST_TRAJECTORY)
#include "diff_solve.hpp"
struct SchemeDiff{
    template <typename T,template <typename> typename GridType>
    static inline Function::Scheme<3> interpolate(const GridType<T> &Grid,T x){
        size_t i = Grid.pos(x);
        if(i<1)
            i = 1;
        if(i > Grid.size()-2)
            i =  Grid.size()-2;

        auto h = Grid.at(i+1)-Grid.at(i);
        auto a1 = (x - Grid.at(i))/h;
        T w = (x-Grid.at(i))/(Grid.at(i+1)-Grid.at(i));
        auto W_left = (-0.5+a1)/h;
        auto W_right= (0.5+a1)/h;
        return Function::Scheme<3>({{i-1,i,i+1},{W_left,-W_left-W_right,W_right}});
    }
};


int main()
{
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");

    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator>
            PhiL(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,SchemeDiff>
            F_grid(R,phi);
    auto F = [&F_grid](double x){return (x<=1 ? F_grid(x): -1/(x*x));};

    PVAR(F_grid[0]);

    double r_edge = 0.4;
    double v_nd = 1.4;
    double E_nd = v_nd*v_nd - PhiC(r_edge);
    PVAR(E_nd);
    PVAR(v_nd*v_nd - PhiL(r_edge));
    double lp = r_edge*v_nd;
    PVAR(lp);

    size_t Nb = 101;
    auto TI = CalculateTrajectory(PhiC,E_nd,lp,Nb);


    size_t N_traj = 2001;
    double h = TI.T_in/(Nb-1)/2;



    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator>
            ExactTrajectory(Function::UniformGrid<double>(0,h*(N_traj-1),N_traj),
                            vmap([](vec4 x){return sqrt(x.t*x.t+x.x*x.x);},solve_diff(
                                    vec4(r_edge,0,0,v_nd),
                                [&F](vec4 x){
                                double r = sqrt(x.t*x.t+x.x*x.x);
                                double n1 = (r != 0 ? x.t/r : 1);
                                double n2 = (r != 0 ? x.x/r : 0);
                                double f = F(r);
                                return vec4(x.y,x.z,0.5*f*n1,0.5*f*n2);
                            }
                                ,N_traj,h
                                ))
            );

    //PVAR(ExactTrajectory.toString());
    //PVAR(TI.Trajectory.toString());

    PVAR(TI.T_in);
    PVAR(TI.T_out);
    Gnuplot gp;
    gp.plotd(PhiC.toString(),"title \"phi\" with lines");
    gp.plotd(Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(R,F).toString(),
             "title \"F\" with lines");
    gp.plot(std::to_string(E_nd),"title \"E\"");
    gp.plotd(Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(R,[&PhiC,lp](double r)
    {return -PhiC(r)+lp*lp/(r*r);}).toString(),"title \"EffPhi\" with lines");
    gp.command("set yrange [-6:6]");
    gp.show();

    Function::GridFunction<double,Function::VectorGrid<double>,Function::LinearInterpolator> Traj(
                    (TI.Trajectory.values),
                    (Function::VectorGrid<double>(TI.Trajectory.Grid).grid())
                );
    Gnuplot gp2;
    gp2.plotd(TI.Trajectory.toString(), "title \"calculate\" with lines");
    gp2.plotd(ExactTrajectory.toString(), "title \"by method\" with lines");
    gp2.show();

    /*
    std::function<void(const std::string &)> ActualGp;
    ActualGp = std::bind(&Gnuplot::command,&gp,std::placeholders::_1);
    std::string S;
    while(1){
        std::cout << "command>";
        std::getline(std::cin,S);
        if(S == "quit"){
            break;
        }
        else if(S == "change 2")
            ActualGp = std::bind(&Gnuplot::command,&gp2,std::placeholders::_1);
        else if(S == "change 1")
            ActualGp = std::bind(&Gnuplot::command,&gp,std::placeholders::_1);
        else
            ActualGp(S);
    }
    */
    pipe_switcher s;
    s.add_pipe(&gp,&Gnuplot::command,"pot");
    s.add_pipe(&gp2,&Gnuplot::command,"traj");
    pipe_menu menu;
    menu.add_dialog_action([&s](){std::cout << s.tmp_pipe() << ">";});
    menu.add_action("show",[&s,&gp,&gp2](const std::string &){
        if(s.tmp_pipe() == "pot")
            gp.show();
        else
            gp2.show();
    });
    menu.add_action("pipes",[&s](const std::string &){print(s.pipe_names());});
    menu.add_action(&s,&pipe_switcher::process_string,[](const std::string &){});
    menu.exec();
    return 0;
}
#elif defined(TEST_SF_HISTO)

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



int main(void){
    const auto& BM = BodyModel::fromFile(2.056e-3,"D:/Important/articles/DMFramework/celstial_models/solar_model.dat");
    double Emin = -BM["phi"][0];
    auto G = [](){return rand()/(RAND_MAX+0.);};
    auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);

    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Therm(R,BM["Temp"]*1e-14);
    Function::GridFunction<double,Function::UniformGrid<double>,Function::LinearInterpolator> Vesc(R,BM["Vesc"]);

    size_t NE = 30+1;
    size_t NLmax = 30 + 1;
    Function::UniformGrid<double> E_grid(Emin,0,NE);
    EL_Histo<double,Function::UniformGrid<double>,Function::UniformGrid<double>> ELH(E_grid,Function::FunctorM(BM,maxLnd));
    double Lmax = max(ELH.Lmax.values);
    for(size_t i=0;i<ELH.values.size();++i){
        ELH.values[i] = Function::Histogramm<double,Function::UniformGrid<double>>(
                    Function::UniformGrid<double>(0,1,2 + ELH.Lmax.values[i]/Lmax*NLmax));
    }
    PVAR(ELH.gridStr());

    Function::GridFunction<decltype(ELH),Function::UniformGrid<double>,
            Function::LinearInterpolator,Function::UniformGrid<double>,Function::LinearInterpolator>
            ScatterPreHisto(E_grid,Vector(E_grid.size(),[&ELH](size_t i){
                                size_t N = ELH.values[std::min(i,ELH.values.size()-1)].Grid.size();
                                return Function::GridFunction<decltype(ELH),Function::UniformGrid<double>,Function::LinearInterpolator>
                                (Function::UniformGrid<double>(0,std::max(ELH.Lmax.values[i],1e-10),N),
                                std::vector<decltype(ELH)>(N,ELH));
                            }));

    auto EvaporationPreHisto = Function::GridFunctionCreator2<Function::LinearInterpolator>::Create(E_grid,
        Vector(E_grid.size(),[&ELH](size_t i){
            size_t N = ELH.values[std::min(i,ELH.values.size()-1)].Grid.size();
            return Function::GridFunctionCreator1<Function::LinearInterpolator>::Create(
                Function::UniformGrid<double>(0,ELH.Lmax.values[i],N),std::vector<double>(N,0));
        }));



    double mk = 3;
    double delta_mk = 1e-4;



    for(size_t i=0;i<ScatterPreHisto.Grid.size();++i){
        for(size_t j = 0;j<ScatterPreHisto.values[i].Grid.size();++j){
            auto TI = CalculateTrajectory(PhiC,ScatterPreHisto.Grid[i],ScatterPreHisto.values[i].Grid[j],10);
            TrajectoryIntegral(G,ScatterPreHisto.Grid[i],ScatterPreHisto.values[i].Grid[j],TI,
                               BM.VescMin(),Therm,Vesc,ScatterPreHisto.values[i].values[j],EvaporationPreHisto.values[i].values[j],
                               mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,10000);
        }
    }

    SupressFactor(ELH,mk,delta_mk,ELASTIC,PROTON,PhiFactorSS,BM,"H",10000);

    PVAR(vmap([](const auto &x){return x.summ();},ScatterPreHisto.AllValues()));
    PVAR(ScatterPreHisto.gridStr());

    auto MatHisto = Function::GridExtractorType<std::vector<double>,decltype(ScatterPreHisto)>::type::sameGrid(ScatterPreHisto);
    for(size_t i=0;i<ScatterPreHisto.Grid.size();++i){
        for(size_t j=0;j<ScatterPreHisto.values[i].Grid.size();++j){
            MatHisto.values[i].values[j] = ScatterPreHisto.values[i].values[j].AllValues();
        }
    }
    auto HMat = Function::GridExtractorType<std::vector<double>,decltype(ELH)::HBase>::type::sameGrid(ELH);
    for(size_t i=0;i<HMat.Grid.size()-1;++i){
        double E_av = 0.5*(HMat.Grid[i]+HMat.Grid[i+1]);
        for(size_t j=0;j<HMat.values[i].Grid.size()-1;++j){
            double L_av = 0.5*(HMat.values[i].Grid[j]+HMat.values[i].Grid[j+1])*ELH.Lmax.values[i];
            HMat.values[i].values[j] = MatHisto(E_av,L_av);
        }
    }
    auto SmatT = HMat.AllValues();
    auto Rmat = Rmatrix(SmatT,10.0);


    Gnuplot distrib;
    distrib.show_cmd = "splot";
    distrib.plotd(ELH.toFunction().toString(),"with lines title \"distrib\"");


    auto ResDistrib = ELH;

    auto ResVector = dot(MatrixPow(Rmat,1000),ELH.AllValues());
    ResDistrib.loadIter(ResVector.begin(),ResVector.end());

    PVAR(ELH.summ());
    PVAR(vector_sum(ResVector));

    distrib.plotd(ResDistrib.toFunction().toString(),"with lines title \"final distrib\"");

    distrib.show();
    std::cin.get();

    //PVAR(SmatT);
    //PVAR(EvaporationPreHisto.toString());
    //PVAR(ScatterPreHisto.gridStr());
    return 0;
}
#endif


