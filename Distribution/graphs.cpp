#include <iostream>
#include "functions.hpp"
#include <cmath>
#include <fstream>

#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
#include "func/arg_parser.hpp"
#include "func/dens_grid.hpp"
#include "func/matrix_functions.hpp"
//#include "boost/type_index.hpp"
/*
std::string process_space(std::string const &S){
    std::regex fname("\\s.(\"([\\w].)\")?(\\w.)\\s.");
}
*/


template <typename T>
auto get_DLHL(const  boost::property_tree::ptree& tree){
    auto FPath = [&](const std::string & tree_path){
        return config_path_from(tree.pgets(tree_path),tree.pgets("config_path")).string();
    };
    MatrixFormat  MF = BINARY;
    if(ptree_gets(tree,"format") == "text")
        MF = TEXT;
    std::vector<double> D_L = vmap([](T x)->double{return x;},
            std::get<0>(
                LoadMatrix<T>(
                    FPath("dl_final"),MF)));
    std::vector<double> D_H = vmap([](T x)->double{return x;},
            std::get<0>(
                LoadMatrix<T>(
                    FPath("dh_final"),MF)));
    return _T(D_L,D_H);
}

template <typename ...Args>
auto ELH_toFunction(const EL_Histo<Args...> &H,size_t EN,size_t LN,double esgm,double lsgm){
    Function::UniformGrid<float> EG(H.Grid._a(),H.Grid._b(),EN);
    Function::UniformGrid<float> LG(0,H.Lmax.values.back(),LN);
    return Function::GridFunction<double,Function::UniformGrid<float>,Function::LinearInterpolator,
            Function::UniformGrid<float>,
            Function::LinearInterpolator>(EG,[&](float e){
                return Function::GridFunction<double,Function::UniformGrid<float>,
                        Function::LinearInterpolator>(LG,[&](float l){
                                double summ = 0;
                                for(size_t i=0;i<H.Grid.size()-1;++i){
                                    float _e = 0.5*(H.Grid[i]+H.Grid[i+1]);
                                    for(size_t j=0;j<H.values[i].Grid.size()-1;++j){
                                        float _l = 0.5*(H.values[i].Grid[j]+H.values[i].Grid[j+1])*H.Lmax(_e);
                                        summ += exp(-(e-_e)*(e-_e)/(2*esgm*esgm)-
                                                    (l-_l)*(l-_l)/(2*lsgm*lsgm)
                                                    )/(2*M_PI*esgm*lsgm);
                                    }
                                }
                                return summ;
                            }
                            );
            });

}

inline auto CrossSquare(double x0,double x1,double _x0,double _x1){
    return _T(std::max(x0,_x0),std::min(x1,_x1));
}
inline bool IsIntersectSquare(double x0,double x1,double _x0,double _x1){
    return std::max(x0,_x0) < std::min(x1,_x1);
}




template <typename ...Args,typename KernelMeasure4p>
auto ELH_toFunction_v2(const EL_Histo<Args...> &H,size_t EN,size_t LN,double se,double sl,const KernelMeasure4p & mu){
    double Emin = H.Grid._a();
    double Emax = H.Grid._b();
    double he =  (H.Grid._b()-H.Grid._a())/(EN-3);
    double hl =  H.Lmax.values.back()/(LN-3);
    Function::UniformGrid<float> EG(H.Grid._a()-he,H.Grid._b()+he,EN);
    Function::UniformGrid<float> LG(-hl,H.Lmax.values.back()+hl,LN);
    Function::GridFunction<double,Function::UniformGrid<float>,Function::LinearInterpolator,
            Function::UniformGrid<float>,
            Function::LinearInterpolator> F(EG,LG);
    F.map([&](float e,float l){
                /*return Function::GridFunction<double,Function::UniformGrid<float>,
                        Function::LinearInterpolator>(LG,[&](float l){*/
                                double summ = 0;

                                double emin = std::clamp(e-se,Emin,Emax);
                                double emax = std::clamp(e+se,Emin,Emax);
                                //if(e-se)
                                if(emin==emax){
                                    return 0.0;
                                }
                                double Lmin = std::clamp((l-sl)/H.Lmax(emax),0.0,1.0);
                                double Lmax = (l+sl);
                                if(Lmax==0.0){
                                    return 0.0;
                                }
                                Lmax = std::clamp(Lmax/H.Lmax(emin),0.0,1.0);
                                size_t i0 = H.Grid.pos(emin);
                                size_t i1 = std::min(H.Grid.pos(emax),H.Grid.size()-2);
                                for(size_t i=i0;i<=i1;++i){
                                    size_t j0 = H.values[i].Grid.pos(Lmin);
                                    double de0 = (H.Grid[i]-e)/se;
                                    double de1 = (H.Grid[i+1]-e)/se;
                                    double lmax0 = H.Lmax[i];
                                    double lmax1 = H.Lmax[i+1];
                                    size_t j1 = std::min(H.values[i].Grid.pos(Lmax),H.values[i].Grid.size()-2);
                                    for(size_t j=j0;j<=j1;++j){
                                        double dl00 = (H.values[i].Grid[j]*lmax0-l)/sl;
                                        double dl01 = (H.values[i].Grid[j+1]*lmax0-l)/sl;
                                        double dl10 = (H.values[i].Grid[j]*lmax1-l)/sl;
                                        double dl11 = (H.values[i].Grid[j+1]*lmax1-l)/sl;
                                        EL_rect R(de0,de1,dl00,dl01,dl10,dl11);
                                        summ += mu(R)*
                                                H.values[i].values[j]/EL_rect_measure(R)*sl*se;

                                        if(mu(R) < 0 ){
                                            std::cout << R <<std::endl;
                                            mu(R);
                                        }
                                    }
                                }
                                return summ;
                            }
                            );
    return F;
}

struct Measure_11{
    inline double operator()(const EL_rect&R)const noexcept{
        auto Intersections = R.SquareIntersection(-1.0,1.0,-1.0,1.0);
        return EL_CascadeIntegration(EL_column_measure,Intersections)/4;
    }
};

int main(int argc,char **argv){
    boost::property_tree::ptree cmd_params;
    auto cmd_parse_log = parse_command_line_v1(argc,argv,cmd_params);
    //std::cout << cmd_params.get<std::string>(0) << std::endl;

    auto Poll = [&](const std::string & property){
        if(!ptree_contain(cmd_params,property)){
            std::string fname;
            std::cout << "input" << property <<" filename: ";
            std::getline(std::cin,fname);
            cmd_params.put(property,fname);
        }
    };
    auto FPath = [&](const std::string & tree_path){
        return config_path_from(cmd_params.pgets(tree_path),cmd_params.pgets("config_path")).string();
    };
    auto HistoToFunc = [&](const auto &ELH){
        return ELH_toFunction_v2(ELH,100,100,0.1,0.1,Measure_11());
    };


    Poll("l_grid");
    Poll("h_grid");
    Poll("dl_final");
    Poll("dh_final");

    auto ELH_L =
            EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>>::load
            (std::ifstream(FPath("l_grid")),ptree_gets(cmd_params,"grid_format") == "binary");

    auto ELH_H =
            EL_Histo<double,Function::VectorGrid<double>,Function::VectorGrid<double>>::load
            (std::ifstream(FPath("h_grid")),ptree_gets(cmd_params,"grid_format") == "binary");

    auto get_dldh = (ptree_gets(cmd_params,"type") == "float" ? get_DLHL<float> : get_DLHL<double>);

    auto [DL,DH] = get_dldh(cmd_params);
    ELH_L.loadIter(DL.begin(),DL.end());
    ELH_H.loadIter(DH.begin(),DH.end());

    PVAR(HistoToFunc(ELH_L).toString());



    Gnuplot gp_3d3(GNUPLOT_PATH);
    gp_3d3.command("set pm3d map");
    gp_3d3.show_cmd = "splot";
    gp_3d3.plotd(HistoToFunc(ELH_H).toString(),"with pm3d title \"H distrid\"");
    //gp_3d.plotd(ELH_H.toFunction().toString(),"with pm3d title \"H distrid\"");
    gp_3d3.show();

    Gnuplot gp_3d4(GNUPLOT_PATH);
    gp_3d4.command("set pm3d map");
    gp_3d4.show_cmd = "splot";
    gp_3d4.plotd(HistoToFunc(ELH_L).toString(),"with pm3d title \"L distrid\"");
    //gp_3d.plotd(ELH_H.toFunction().toString(),"with pm3d title \"H distrid\"");
    gp_3d4.show();

    Gnuplot gp_traj;
    if(ptree_contain(cmd_params,"nl") and ptree_contain(cmd_params,"nh")){
        gp_traj.plotf(cmd_params.pgets("nl"),"with lines title \"nl\"");
        gp_traj.plotf(cmd_params.pgets("nh"),"with lines title \"nh\"");
        gp_traj.show();
    }
    wait();
	return 0;
}
