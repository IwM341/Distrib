#include "functions.hpp"
#include <cmath>
#include <fstream>

#include <time.h>
#include "func/el_histo.hpp"
#include <utils>
#include <random>
#include <numeric>
//#include "func/arg_parser.hpp"
//#include "func/dens_grid.hpp"
//#include "func/matrix_functions.hpp"
#include "func/arg_parser.hpp"
#include <sstream>
#include "func/move_to_go.hpp"
#include "grob/templates.hpp"
#include "grob/grid_objects.hpp"
#include "grob/object_serialization.hpp"
#include "grob/serialization.hpp"
#include "func/matrix_functions.hpp"
#include "func/load_histo.hpp"
#include "grob/csv_io.hpp"
#include "../factors/factors.hpp"


template <typename L_E_Functype,typename Phi_FuncType,
          typename Generator,typename phi_ann_t,typename...HistoArgs>
auto annihilation_r_distrib(grob::Histogramm<HistoArgs...> const&Histo_L,
                            grob::Histogramm<HistoArgs...> const& Histo_H,
                            L_E_Functype const & L_E,
                            Phi_FuncType const &phi,
                            Generator && G,phi_ann_t Phi_Ann_Func,size_t N_mk = 1,
                           const size_t N_r = 1000,double R_cut = 10) noexcept{

        print("start calculations");
        progress_bar progress(Histo_L.size()*Histo_H.size(),100);


        auto Mi = Histo_L.size();
        auto rd_func = grob::make_function_f(grob::GridUniform<double>(0,R_cut,N_r+1),
                              [](size_t i){return 0;});

        for(size_t i=0;i<Mi;++i)
        {
            auto MI = Histo_L.Grid.FromLinear(i);

            auto EL_rect_L = Histo_L.Grid[MI];
            auto E_L = EL_rect_L.template getR<0>().center();
            auto l_L = EL_rect_L.template getR<1>().center();
            auto L_L =l_L*L_E(E_L);

            double L00_L = EL_rect_L.template get<1,0b00>()*L_E(EL_rect_L.left());
            double L01_L = EL_rect_L.template get<1,0b10>()*L_E(EL_rect_L.left());
            double L10_L = EL_rect_L.template get<1,0b01>()*L_E(EL_rect_L.right());
            double L11_L = EL_rect_L.template get<1,0b11>()*L_E(EL_rect_L.right());

            TrajectoryPreInfo T_L = CalculatePeriod(phi,E_L,L_L,10);


            double rho_L = Histo_L[MI]/
                    EL_volume(EL_rect_L.left(),EL_rect_L.right(),
                                  L00_L,L01_L,L10_L,L11_L,
                                  T_L.T);


            for(grob::MultiIndex<2> MJ = {0,0};
                MJ!=Histo_H.Grid.MultiSize();
                Histo_H.Grid.MultiIncrement(MJ))
            {
                auto EL_rect_H = Histo_H.Grid[MJ];
                auto E_H = EL_rect_H.template getR<0>().center();
                auto l_H = EL_rect_H.template getR<1>().center();
                auto L_H =l_H*L_E(E_H);



                double L00_H = EL_rect_H.template get<1,0b00>()*L_E(EL_rect_H.left());
                double L01_H = EL_rect_H.template get<1,0b10>()*L_E(EL_rect_H.left());
                double L10_H = EL_rect_H.template get<1,0b01>()*L_E(EL_rect_H.right());
                double L11_H = EL_rect_H.template get<1,0b11>()*L_E(EL_rect_H.right());

                TrajectoryPreInfo T_H = CalculatePeriod(phi,E_H,L_H,10);
                double rho_H = Histo_H[MJ]/
                        EL_volume(EL_rect_H.left(),EL_rect_H.right(),
                                      L00_H,L01_H,L10_H,L11_H,
                                      T_H.T);

                const size_t rd_sz = rd_func.Grid.size();
                #pragma omp parallel for shared(rd_func)
                for(size_t j=0;j<rd_sz;++j){
                    double r = rd_func.Grid[i];
                    double rge_dens = (4*M_PI*4*M_PI/3)*r*r;

                    double d3v_L= d_3_v_mes(EL_rect_L,L_E,phi,r);
                    double d3v_H= d_3_v_mes(EL_rect_H,L_E,phi,r);


                    double v_L = del_nan(sqrt(E_L+phi(r)));
                    double v_t_L = del_nan(L_L/r);
                    double v_r_L = del_nan(sqrt(v_L*v_L-v_t_L*v_t_L));

                    double v_H = del_nan(sqrt(E_H+phi(r)));
                    double v_t_H = del_nan(L_H/r);
                    double v_r_H = del_nan(sqrt(v_H*v_H-v_t_H*v_t_H));
                    double fac  = 0;
                    for(size_t k = 0;k<N_mk;++k){
                        double v_L_dot_v_H = v_t_L*v_t_H*cos(G()*M_PI) + v_r_L*v_r_H;

                        double Vsum = sqrt(v_L*v_L+v_H*v_H+2*v_L_dot_v_H);
                        double Vdiff = sqrt(v_L*v_L+v_H*v_H-2*v_L_dot_v_H);
                        fac += rge_dens*d3v_L*d3v_H*rho_L*rho_H/(N_mk*2)*(
                                    Phi_Ann_Func(Vsum) + Phi_Ann_Func(Vdiff));

                    }
                    rd_func.Values[i] += fac;
                }
                progress++;
            }
        }
        return rd_func;
}


int main(int argc, char ** argv){
    print("getting annihilation matrix");
    print("params: ",
          "l_d.grid [filename of L grid]",
          "l_d.values [filename of L values]",
          "h_d.grid [filename of H grid]",
          "h_d.values [filename of H values]",
          "body [filename of body model]",
          "Rcut [default is 10]",
          "out [output filename]",
          "Nmk [default is 1000]");

    boost::property_tree::ptree cmd_params;
    auto ret_key = parse_command_line(argc,argv,cmd_params);
    if(!ret_key.empty()){
        std::cout << "error : " << ret_key<<std::endl;
        return 0;
    }
    boost::filesystem::path programm_path = ptree_gets(cmd_params,"config_path");

    std::string bm_filename;
    double Vesc1;
    try {
        bm_filename = cmd_params.get<std::string>("body");
        Vesc1 = cmd_params.get<double>("Vesc");

        if(ptree_gets(cmd_params,"debug") != ""){
            std::cout << "body filename: " << bm_filename<<std::endl;
            std::cout << "Vesc " << Vesc1<<std::endl;
        }
    } catch (std::exception &e) {
        std::cout << "can't get body model" <<std::endl;
        return 0;
    }
    auto BM = BodyModel::fromFile(Vesc1,
                             config_path_from(bm_filename,programm_path).string());

    auto H_distrib = loadHisto<grob::MultiGridHisto<
                                    grob::GridVectorHisto<double>,
                                    std::vector<grob::GridVectorHisto<double>>
                                >,
                                double>
            (cmd_params.get_child("h_d"),programm_path);

    auto L_distrib = loadHisto<grob::MultiGridHisto<
                                    grob::GridVectorHisto<double>,
                                    std::vector<grob::GridVectorHisto<double>>
                                >,
                                double>
            (cmd_params.get_child("l_d"),programm_path);


    double VescMax_nd = BM.VescMax()/BM.VescMin();
    auto RetHisto = grob::make_histo<double>(
                grob::mesh_grids(
                        grob::GridUniformHisto<double>(0,VescMax_nd,100),
                        grob::GridUniformHisto<double>(0,VescMax_nd,100)
                    )
                );
    auto G = [](){return rand()/(RAND_MAX+1.0);};
    grob::GridUniform<double> r_grid(0,1,BM["Vesc"].size());
    auto VescR = grob::make_function(r_grid,BM["Vesc"]);
    auto phiR = grob::make_function(r_grid,BM["phi"]);
    double Emin = -phiR.Values[0];
    auto LE_func = grob::make_function_f(
                    grob::GridUniform<double>(Emin ,0,1000),
                [&phiR](double e){return maxLndf(phiR,e);}
                );
    size_t Nmk = cmd_params.get<int>("Nmk",1);
    double Rcut = cmd_params.get<double>("Rcut",10);
    size_t Nr = cmd_params.get<int>("Nr",1000);
    auto rd_func = annihilation_r_distrib(L_distrib,H_distrib,LE_func,phiR,G,Phi_Fac_Ann{},Nmk,Nr,Rcut);

    grob::as_csv(rd_func).save(std::ofstream(cmd_params.pgets("out")),6,std::defaultfloat);

    /*
    if(ptree_contain(cmd_params,"csv")){

    }*/



}
