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
#include "../factors/factors.hpp"



template <typename Grid_L_t,typename Grid_H_t,typename L_E_Functype,
          typename Phi_FuncType,typename Generator>
struct scatter_counter{

    Grid_L_t const& Grid_L;
    Grid_H_t const& Grid_H;
    L_E_Functype const & L_E;
    Phi_FuncType const &phi;
    Generator const& G;


    scatter_counter(Grid_L_t const& Grid_L,
                    Grid_H_t const& Grid_H,
                             L_E_Functype const & L_E,
                    Phi_FuncType const & phi,Generator const& G)noexcept:
                    Grid_L(Grid_L),Grid_H(Grid_H),L_E(L_E),phi(phi),
                    G(G){}

    /*!
     * \f$\cfrac{dN}{dKd\nu d\sigma}\f$
     * \return ANN_HL vector of values with NH major dim
    */
    template <typename phi_ann_t>
    std::vector<double> annihilation_matrix(phi_ann_t Phi_Ann_Func,
                               const size_t N_mk = 10000,double R_cut = 10)const noexcept{
        auto r_gen = [this,R_cut](){
            double r = G()*R_cut;
            return MC::MCResult<double>(r,4*M_PI*r*r);
        };

            print("start calculations");
            size_t NH =Grid_H.size();
            size_t NL =Grid_L.size();
            progress_bar prog_lh(NL*NH,100);
            std::vector<double> ANN_HL(NL*NH,0);
            #pragma omp parallel for shared(ANN_HL)
            for(size_t i=0;i<NL;++i)
            {
                auto Histo = grob::make_histo_ref(Grid_H,
                                                  grob::make_slice(ANN_HL,i*NH,NH)
                                               );

                auto MI = Grid_L.FromLinear(i);

                auto EL_rect_L = Grid_L[MI];
                auto E_L = EL_rect_L.template getR<0>().center();
                auto l_L = EL_rect_L.template getR<1>().center();
                auto L_L =l_L*L_E(E_L);

                double L00_L = EL_rect_L.template get<1,0b00>()*L_E(EL_rect_L.left());
                double L01_L = EL_rect_L.template get<1,0b10>()*L_E(EL_rect_L.left());
                double L10_L = EL_rect_L.template get<1,0b01>()*L_E(EL_rect_L.right());
                double L11_L = EL_rect_L.template get<1,0b11>()*L_E(EL_rect_L.right());

//                double T00_L = CalculatePeriod(phi,EL_rect_L.left(),L00_L,10);
//                double T01_L = CalculatePeriod(phi,EL_rect_L.left(),L01_L,10);
//                double T10_L = CalculatePeriod(phi,EL_rect_L.right(),L10_L,10);
//                double T11_L = CalculatePeriod(phi,EL_rect_L.right(),L11_L,10);
                TrajectoryPreInfo T_L = CalculatePeriod(phi,E_L,L_L,10);
                double rho_L = 1/
                        EL_volume(EL_rect_L.left(),EL_rect_L.right(),
                                      L00_L,L01_L,L10_L,L11_L,
                                      T_L.T);


                for(grob::MultiIndex<2> MJ = {0,0};
                    MJ!=Grid_H.MultiSize();
                    Grid_H.MultiIncrement(MJ))
                {
                    size_t j = Grid_H.LinearIndex(MJ);
                    auto EL_rect_H = Grid_H[MJ];
                    auto E_H = EL_rect_H.template getR<0>().center();
                    auto l_H = EL_rect_H.template getR<1>().center();
                    auto L_H =l_H*L_E(E_H);



                    double L00_H = EL_rect_H.template get<1,0b00>()*L_E(EL_rect_H.left());
                    double L01_H = EL_rect_H.template get<1,0b10>()*L_E(EL_rect_H.left());
                    double L10_H = EL_rect_H.template get<1,0b01>()*L_E(EL_rect_H.right());
                    double L11_H = EL_rect_H.template get<1,0b11>()*L_E(EL_rect_H.right());

                    TrajectoryPreInfo T_H = CalculatePeriod(phi,E_H,L_H,10);
                    double rho_H = 1/
                            EL_volume(EL_rect_H.left(),EL_rect_H.right(),
                                          L00_H,L01_H,L10_H,L11_H,
                                          T_H.T);

                    double & HL_Value = ANN_HL[i*NH+j];

                    double r_min = std::max(T_L.rmin,T_H.rmin);
                    double r_max = std::min(T_L.rmax,T_H.rmax);
                    if(r_min < r_max){
                        for(size_t i=0;i<N_mk;++i){


                            double r = r_min + G()*(r_max-r_min);
                            double rge_dens = (4*M_PI*4*M_PI/3)*r*r*(r_max-r_min);

                            double d3v_L= d_3_v_mes(EL_rect_L,L_E,phi,r);
                            double d3v_H= d_3_v_mes(EL_rect_H,L_E,phi,r);


                            double v_L = del_nan(sqrt(E_L+phi(r)));
                            double v_t_L = del_nan(L_L/r);
                            double v_r_L = sqrt(v_L*v_L-v_t_L*v_t_L);

                            double v_H = del_nan(sqrt(E_H+phi(r)));
                            double v_t_H = del_nan(L_H/r);
                            double v_r_H = del_nan(sqrt(v_H*v_H-v_t_H*v_t_H));

                            double v_L_dot_v_H = v_t_L*v_t_H*cos(G()*M_PI) + v_r_L*v_r_H;
                            double Vsum = sqrt(v_L*v_L+v_H*v_H+2*v_L_dot_v_H);
                            double Vdiff = sqrt(v_L*v_L+v_H*v_H-2*v_L_dot_v_H);
                            double fac = rge_dens*d3v_L*d3v_H*rho_L*rho_H/(N_mk*2);
                            HL_Value += fac*Phi_Ann_Func(Vdiff);
                            HL_Value += fac*Phi_Ann_Func(Vsum);
                        }
                    }
                    prog_lh++;
                }
            }
        return ANN_HL;
    }
};

template <typename...Args>
auto sc_cnt(Args &&...args){
    return scatter_counter<std::decay_t<Args>...>(
                std::forward<Args>(args)...
                );
}


int main(int argc, char ** argv){
    print("getting annihilation matrix");
    print("params: ",
          "l_grid [filename of L grid]",
          "h_grid [filename of H grid]",
          "body [filename of body model]",
          "Rcut [default is 10]",
          "ann [out annihilation col major matrix A_HL]",
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

    typedef grob::MultiGridHisto<
                grob::GridVectorHisto<double>,
                std::vector<grob::GridVectorHisto<double>>
            > GridType;

    stools::PtreeSerializator<boost::property_tree::ptree> S{};
    boost::property_tree::ptree Grid_L_p,Grid_H_p;

    boost::property_tree::read_json(
                config_path_from(
                    cmd_params.pgets("l_grid"),
                    programm_path
                    ).string(),Grid_L_p
                );
    boost::property_tree::read_json(
                config_path_from(
                    cmd_params.pgets("h_grid"),
                    programm_path
                    ).string(),Grid_H_p
                );

    auto L_Grid = stools::DeSerialize<GridType>(Grid_L_p,S);
    auto H_Grid = stools::DeSerialize<GridType>(Grid_H_p,S);


    auto G = [](){return rand()/(RAND_MAX+0.0);};
    grob::GridUniform<double> r_grid(0,1,BM["Vesc"].size());
    auto VescR = grob::make_function(r_grid,BM["Vesc"]);
    auto phiR = grob::make_function(r_grid,BM["phi"]);
    double Emin = -phiR.Values[0];
    auto LE_func = grob::make_function_f(
                    grob::GridUniform<double>(Emin ,0,1000),
                [&phiR](double e){return maxLndf(phiR,e);}
                );
    auto sc = sc_cnt(L_Grid,H_Grid,LE_func,phiR,G);
    auto Ann_Matrix = sc.annihilation_matrix(Phi_Fac_Ann{},cmd_params.get<int>("Nmk",100000),
                             cmd_params.get<double>("Rcut",10));


    saveMatrix(Ann_Matrix.data(),L_Grid.size(),H_Grid.size(),
               config_path_from(cmd_params.pgets("ann","ann_HL.bmat"),programm_path).string(),MatrixFormat::BINARY);




}
