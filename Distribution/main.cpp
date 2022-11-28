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

#define SWITCH1 1
#define ON(x,reason) x
#define OFF(x,reason) x##reason


//TEST SUCCEED
//https://colab.research.google.com/drive/1qUEolPBBhXua3gXKSBi_QDfXSrb03oyT#scrollTo=4ld9zqreuCTh
int OFF(main,test_borders)(void){
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


//TEST SUCCEED
//https://colab.research.google.com/drive/12ahUMNi_Isr9gfgwIeDUMncZzqwP_6Rv?usp=sharing
int OFF(main,test_distribution)(void){
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
    return 0;
}


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


int OFF(main,test_distrib)()
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







