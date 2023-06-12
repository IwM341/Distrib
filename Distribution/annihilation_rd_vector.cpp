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
#include "templates.hpp"
#include "grid_objects.hpp"
#include "object_serialization.hpp"
#include "serialization.hpp"
#include "func/matrix_functions.hpp"
#include "func/load_histo.hpp"
#include "csv_io.hpp"



/*!
 * \f$\cfrac{dN}{dKd\nu d\sigma}\f$
 * \return ANN vector of values
*/
template <typename L_E_Functype,
          typename Phi_FuncType,typename VescFunctype,
          typename Generator,typename phi_ann_t,typename ...HistoArgs>
auto annihilation_rd_vector(grob::Histogramm<HistoArgs...> const&Histo,
                                        L_E_Functype const & L_E,VescFunctype const &VescR,
                                        double Vdisp,double U0_,
                                        Phi_FuncType const &phi,Generator const& G,
                                        phi_ann_t Phi_Ann_Func,
                           const size_t N_mk = 10000,const size_t N_r = 1000,double R_cut = 20) {

        auto rd_func = grob::make_function_f(grob::GridUniform<double>(0,R_cut,N_r+1),
                          [](size_t i){return 0.0;});

        print("start calculations");
        size_t N =Histo.Grid.size();

        progress_bar prog_lh(N,100);
        for(size_t i=0;i<N;++i)
        {

            auto MI = Histo.Grid.FromLinear(i);

            auto EL_rect_L = Histo.Grid[MI];
            auto E_L = EL_rect_L.template getR<0>().center();
            auto l_L = EL_rect_L.template getR<1>().center();
            auto L_L =l_L*L_E(E_L);

            double L00_L = EL_rect_L.template get<1,0b00>()*L_E(EL_rect_L.left());
            double L01_L = EL_rect_L.template get<1,0b10>()*L_E(EL_rect_L.left());
            double L10_L = EL_rect_L.template get<1,0b01>()*L_E(EL_rect_L.right());
            double L11_L = EL_rect_L.template get<1,0b11>()*L_E(EL_rect_L.right());


            TrajectoryPreInfo T_L = CalculatePeriod(phi,E_L,L_L,10);
            double rho_L = 1/
                    EL_volume(EL_rect_L.left(),EL_rect_L.right(),
                                  L00_L,L01_L,L10_L,L11_L,
                                  T_L.T);


            double r_min = T_L.rmin;
            double r_max = std::min(T_L.rmax,R_cut);

            const size_t Nr = rd_func.Grid.size();
            #pragma omp parallel for
            for(size_t j=0;j<Nr;++j){
                for(size_t k = 0;k<N_mk;++k){

                    double r = r_min + G()*(r_max-r_min);
                    double rge_dens = (4*M_PI)*r*r*(r_max-r_min);

                    double d3v_L= d_3_v_mes(EL_rect_L,L_E,phi,r);

                    double v_L = sqrt(E_L+phi(r));

                    auto [v_H,vh_rd] = Velocity(G,VescR(r),Vdisp,U0_);

                    v_H.z -= v_L;
                    double Vdiff = v_H.norm();
                    double fac = rge_dens*d3v_L*vh_rd*rho_L/(N_mk)*Phi_Ann_Func(Vdiff);
                    rd_func[j] += fac;
                }
            }
        }
    return rd_func;
}


int main(int argc, char ** argv){
    print("getting annihilation matrix");
    printd("\n","params: ",
          "distrib.grid [filename of L grid]",
          "distrib.values [filename of L values]",
          "body [filename of body model]",
          "Vbody [velocity of celestial object]",
          "Vdisp [disp velocity in halo]",
          "Rcut [default is 10]",
          "out [output filename]",
          "Nr [default is 1000]",
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

    auto M_distrib = loadHisto<grob::MultiGridHisto<
                                    grob::GridVectorHisto<double>,
                                    std::vector<grob::GridVectorHisto<double>>
                                >,
                                double>
            (cmd_params.get_child("distrib"),programm_path);


    double VescMax_nd = BM.VescMax()/BM.VescMin();
    auto G = [](){return rand()/(RAND_MAX+0.0);};
    grob::GridUniform<double> r_grid(0,1,BM["Vesc"].size());
    auto VescR = grob::make_function(r_grid,BM["Vesc"]);
    auto phiR = grob::make_function(r_grid,BM["phi"]);
    double Emin = -phiR.Values[0];
    auto LE_func = grob::make_function_f(
                    grob::GridUniform<double>(Emin ,0,1000),
                [&phiR](double e){return maxLndf(phiR,e);}
                );
    double V_body = ptree_condition(cmd_params,"Vbody",0.73e-3);
    double V_disp = ptree_condition(cmd_params,"Vdisp",0.52e-3);
    auto rd_func = annihilation_rd_vector(M_distrib,LE_func,VescR,V_disp,V_body,phiR,G,[](double){return 1;},
                             cmd_params.get<int>("Nmk",10000),
                             cmd_params.get<double>("Rcut",10));

    grob::as_csv(rd_func).save(std::ofstream(cmd_params.pgets("out")),6,std::defaultfloat);

    /*
    if(ptree_contain(cmd_params,"csv")){

    }*/



}
