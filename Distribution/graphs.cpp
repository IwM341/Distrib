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
    MatrixFormat  MF = MatrixFormat::BINARY;
    if(ptree_gets(tree,"format") == "text")
        MF = MatrixFormat::TEXT;
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
                                                H.values[i].values[j];//EL_rect_measure(R)/(sl*se);
                                        if(std::isnan(summ)){
                                            PVAR(e);
                                            PVAR(l);
                                            PVAR(H.Grid[i]);
                                            PVAR(H.values[i].Grid[j]);
                                        }
                                    }
                                }
                                return summ/*/std::max(l,1e-6f)*/;
                            }
                            );
    return F;
}
template <typename PhiType,typename ...Args>
auto ELH_toFunction_v3(const EL_Histo<Args...> &H,const PhiType & phi,size_t EN,size_t LN){
    double Emin = H.Grid._a();
    double Emax = H.Grid._b();
    Function::UniformGrid<float> EG(H.Grid._a(),H.Grid._b(),EN);
    Function::UniformGrid<float> LG(0,H.Lmax.values.back(),LN);
    Function::GridFunction<double,Function::UniformGrid<float>,Function::LinearInterpolator,
            Function::UniformGrid<float>,
            Function::LinearInterpolator> F(EG,LG);
    F.map([&](float e,float l){
        return ELH_density(H,phi,e,l);
    });
    return F;
}

template <typename ...Args>
double EL_density1(const EL_Histo<Args...> &H,double E,double L){
    if(E < H.Grid._a() or E > H.Grid._b()){
        return .0;
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

    return Value/(0.5*(E1-E0)*(L01+L11-L00-L10));
}
template <typename ...Args>
double EL_density3(const EL_Histo<Args...> &H,double E,double L){
    if(E < H.Grid._a() or E > H.Grid._b()){
        return .0;
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

    double l0 = H.values[i].Grid[j];
    double l1 = H.values[i].Grid[j+1];

    return Value/((l1-l0)*(E1-E0));
}

template <typename ...Args>
double EL_density2(const EL_Histo<Args...> &H,double E,double L){
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

    if(E < H.Grid[i] or E > H.Grid[i+1] or
            l < H.values[i].Grid[j] or l > H.values[i].Grid[j+1]){
        PVAR(E);
        PVAR(l);
        PVAR(i);
        PVAR(j);
    }
    double Value = H.values[i].values[j];
    /*
    if(E > -0.2)
        std::cout << "i = " << i << ", j = " << j <<std::endl;
    */
    return Value;
}

template <typename ...Args>
auto ELH_toFunction_v4(const EL_Histo<Args...> &H,size_t EN,size_t LN){
    double Emin = H.Grid._a();
    double Emax = H.Grid._b();
    Function::UniformGrid<float> EG(H.Grid._a(),H.Grid._b(),EN);
    Function::UniformGrid<float> LG(0,H.Lmax.values.back(),LN);
    Function::GridFunction<double,Function::UniformGrid<float>,Function::LinearInterpolator,
            Function::UniformGrid<float>,
            Function::LinearInterpolator> F(EG,LG);
    F.map([&](float e,float l){
        return EL_density1(H,e,l);
    });
    return F;
}

template <typename ...Args>
auto ELH_toFunction_v5(size_t type,const EL_Histo<Args...> &H,size_t EN,size_t LN){
    double Emin = H.Grid._a();
    double Emax = H.Grid._b();
    Function::UniformGrid<float> EG(H.Grid._a(),H.Grid._b(),EN);
    Function::UniformGrid<float> LG(0,H.Lmax.values.back(),LN);
    Function::GridFunction<double,Function::UniformGrid<float>,Function::LinearInterpolator,
            Function::UniformGrid<float>,
            Function::LinearInterpolator> F(EG,LG);
    if(type == 1){
        F.map([&](float e,float l){
            return EL_density1(H,e,l);
        });
    }
    else{
        F.map([&](float e,float l){
            return EL_density2(H,e,l);
        });
    }
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
    if(!cmd_parse_log.empty()){
        std::cout << cmd_parse_log <<std::endl;
        return 1;
    }
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


    Poll("l_grid");
    Poll("h_grid");
    Poll("dl_final");
    Poll("dh_final");
    Poll("body");
    auto BM = BodyModel::fromFile(1,FPath("body"));
    const auto & phi = BM["phi"];
    Function::UniformGrid<double> R(0,1,phi.size());
    Function::GridFunction<double,Function::UniformGrid<double>,Function::CubicInterpolator>
            PhiC(R,phi);
    std::string method =  cmd_params.get("tofunc","EL");
    auto HistoToFunc = [&](const auto &ELH){
        if(method == "rv")
            return ELH_toFunction_v3(ELH,PhiC,100,100);
        else if(method == "val")
            return ELH_toFunction_v5(3,ELH,201,201);
        else
            return ELH_toFunction_v5(1,ELH,201,201);
        //return ELH_toFunction_v3(ELH,PhiC,100,100);
    };


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

    //PVAR(HistoToFunc(ELH_L).toString());
    auto DB = ELH_H.toFunction();
    auto gp_path = ptree_condition<std::string>(cmd_params,"gnuplot_path","gnuplot");


    auto ELH_L_normed1 = HistoToFunc(ELH_L);
    double maxH_L =0;
    for(size_t i=0;i<ELH_L_normed1.num_of_element();++i){
        maxH_L =  std::max(ELH_L_normed1[i],maxH_L);
    }
    ELH_L_normed1 /= maxH_L;

    double maxH_H =0;
    auto ELH_H_normed1 = HistoToFunc(ELH_H);
    for(size_t i=0;i<ELH_H_normed1.num_of_element();++i){
        maxH_H = std::max(ELH_H_normed1[i],maxH_H);
    }

    ELH_H_normed1 /= maxH_H;

    Gnuplot gp_3d1(gp_path);
    gp_3d1.command("set pm3d map");
    gp_3d1.command("set palette maxcolors 20");
    gp_3d1.show_cmd = "splot";
    gp_3d1.plotd((ELH_L).toFunction().toString(),"title \"L distrib\" with image");
    //gp_3d.plotd(ELH_H.toFunction().toString(),"with pm3d title \"H distrid\"");
    gp_3d1.show();


    Gnuplot gp_3d2(gp_path);
    gp_3d2.command("set pm3d map");
    gp_3d2.command("set palette maxcolors 20");
    gp_3d2.show_cmd = "splot";
    gp_3d2.plotd((ELH_H).toFunction().toString()," title \"H distrib\" with image");
    //gp_3d.plotd(ELH_H.toFunction().toString(),"with pm3d title \"H distrid\"");
    gp_3d2.show();


    Gnuplot gp_3d3(gp_path);
    gp_3d3.command("set pm3d map");
    gp_3d3.command("set palette maxcolors 20");
    gp_3d3.show_cmd = "splot";
    gp_3d3.plotd(ELH_H_normed1.toString()," title \"H distrib dens\" with image ");
    //gp_3d.plotd(ELH_H.toFunction().toString(),"with pm3d title \"H distrid\"");
    gp_3d3.show();

    Gnuplot gp_3d4(gp_path);
    gp_3d4.command("set pm3d map");
    gp_3d4.command("set palette maxcolors 20");
    gp_3d4.show_cmd = "splot";
    gp_3d4.plotd(ELH_L_normed1.toString()," title \"L distrib dens\" with image");
    //gp_3d.plotd(ELH_H.toFunction().toString(),"with pm3d title \"H distrid\"");
    gp_3d4.show();

    auto Histo_E_H = Function::Histogramm<double,Function::VectorGrid<double>>(ELH_H.Grid,
            Vector((ELH_H.Grid.size()-1),[&](size_t i){return ELH_H.values[i].summ();})
            );

    Gnuplot gp_e_h(gp_path);
    gp_e_h.show_cmd = "plot";
    gp_e_h.plotd(Histo_E_H.toFunction().toString(),"with boxes title \"H E distrib\"");
    gp_e_h.show();

    Gnuplot gp_traj;
    if(ptree_contain(cmd_params,"nl") and ptree_contain(cmd_params,"nh")){
        gp_traj.plotf(cmd_params.pgets("nl"),"with lines title \"nl\"");
        gp_traj.plotf(cmd_params.pgets("nh"),"with lines title \"nh\"");
        gp_traj.show();
    }
    pipe_switcher ps;
    ps.add_pipe(&gp_3d1,&decltype(gp_3d1)::command,"L1");
    ps.add_pipe(&gp_3d2,&decltype(gp_3d2)::command,"H1");
    ps.add_pipe(&gp_3d3,&decltype(gp_3d3)::command,"H2");
    ps.add_pipe(&gp_3d4,&decltype(gp_3d4)::command,"L2");
    ps.add_pipe(&gp_e_h,&decltype(gp_e_h)::command,"HE");
    ps.exec();
	return 0;
}
