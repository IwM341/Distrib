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


template <typename Histo_L,typename Histo_H,typename L_E_Functype,typename Phi_FuncType,typename Generator>
struct scatter_counter;

template <typename...Histo_LArgs,typename...Histo_HArgs,typename L_E_Functype,
          typename Phi_FuncType,typename Generator>
struct scatter_counter<grob::Histogramm<Histo_LArgs...>,
                        grob::Histogramm<Histo_HArgs...>,
                        L_E_Functype,
                        Phi_FuncType,Generator>{

    grob::Histogramm<Histo_LArgs...> const& Histo_L;
    grob::Histogramm<Histo_HArgs...> const& Histo_H;
    L_E_Functype const & L_E;
    Phi_FuncType const &phi;
    Generator const& G;


    scatter_counter(grob::Histogramm<Histo_LArgs...> const& Histo_L,
                    grob::Histogramm<Histo_HArgs...> const& Histo_H,
                             L_E_Functype const & L_E,
                    Phi_FuncType const & phi,Generator const& G)noexcept:
                    Histo_L(Histo_L),Histo_H(Histo_H),L_E(L_E),phi(phi),
                    G(G){}

    /**
     * \f$\cfrac{dN}{dKd\nu d\sigma}\f$
    */
    template <typename HistoType>
    void dens_scatter_function(HistoType & Out,
                               const size_t N_mk = 10000,double R_cut = 10)const noexcept{
        auto r_gen = [this,R_cut](){
            double r = G()*R_cut;
            return MC::MCResult<double>(r,4*M_PI*r*r);
        };
            print("start calculations");
            show_prog prog_lh(100);
            prog_lh.show(0);
            size_t total_progress_lh = Histo_L.size()*Histo_H.size();
            size_t curr_progress_lh = 0;

            auto Mi = Histo_L.size();
            #pragma omp parallel for shared(Out)
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

                double T00_L = CalculatePeriod(phi,EL_rect_L.left(),L00_L,10);
                double T01_L = CalculatePeriod(phi,EL_rect_L.left(),L01_L,10);
                double T10_L = CalculatePeriod(phi,EL_rect_L.right(),L10_L,10);
                double T11_L = CalculatePeriod(phi,EL_rect_L.right(),L11_L,10);

                double rho_L = Histo_L[MI]/
                        EL_bin_volume(EL_rect_L.left(),EL_rect_L.right(),
                                      L00_L,L01_L,L10_L,L11_L,
                                      T00_L,T01_L,T10_L,T11_L);


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

                    double T00_H = CalculatePeriod(phi,EL_rect_H.left(),L00_H,10);
                    double T01_H = CalculatePeriod(phi,EL_rect_H.left(),L01_H,10);
                    double T10_H = CalculatePeriod(phi,EL_rect_H.right(),L10_H,10);
                    double T11_H = CalculatePeriod(phi,EL_rect_H.right(),L11_H,10);

                    double rho_H = Histo_H[MJ]/
                            EL_bin_volume(EL_rect_H.left(),EL_rect_H.right(),
                                          L00_H,L01_H,L10_H,L11_H,
                                          T00_H,T01_H,T10_H,T11_H);


                    TrajectoryPreInfo T_L = CalculatePeriod(phi,E_L,L_L,10);
                    TrajectoryPreInfo T_H = CalculatePeriod(phi,E_H,L_H,10);

                    for(size_t i=0;i<N_mk;++i){
                        double r_min = std::max(T_L.rmin,T_H.rmin);
                        double r_max = std::min(T_L.rmax,T_H.rmax);

                        double r = r_min + G()*(r_max-r_min);
                        double rge_dens = 4*M_PI*r*r*(r_max-r_min);

                        double d3v_L= d_3_v_mes(EL_rect_L,L_E,phi,r);
                        double d3v_H= d_3_v_mes(EL_rect_H,L_E,phi,r);


                        double v_L = sqrt(E_L+phi(r));
                        double v_t_L = L_L/r;
                        double v_r_L = sqrt(v_L*v_L-v_t_L*v_t_L);

                        double v_H = sqrt(E_H+phi(r));
                        double v_t_H = L_H/r;
                        double v_r_H = sqrt(v_H*v_H-v_t_H*v_t_H);

                        double v_L_dot_v_H = v_t_L*v_t_H*cos(G()*M_PI) + v_r_L*v_r_H;
                        double Vsum = sqrt(v_L*v_L+v_H*v_H+2*v_L_dot_v_H);
                        double Vdiff = sqrt(v_L*v_L+v_H*v_H-2*v_L_dot_v_H);
                        double fac = rge_dens*d3v_L*d3v_H*rho_L*rho_H/(N_mk*2);
                        #pragma omp critical
                        {
                            Out.put(fac,Vsum,Vdiff);
                            Out.put(fac,Vdiff,Vsum);
                        }
                    }
                    curr_progress_lh++;
                    prog_lh.show(curr_progress_lh/((float)total_progress_lh));
                }
            }


            prog_lh.end();

    }
};

template <typename...Args>
auto sc_cnt(Args &&...args){
    return scatter_counter<std::decay_t<Args>...>(
                std::forward<Args>(args)...
                );
}


int main(int argc, char ** argv){
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
            (cmd_params.get_child("H_D"),programm_path);

    auto L_distrib = loadHisto<grob::MultiGridHisto<
                                    grob::GridVectorHisto<double>,
                                    std::vector<grob::GridVectorHisto<double>>
                                >,
                                double>
            (cmd_params.get_child("L_D"),programm_path);


    double VescMax_nd = BM.VescMax()/BM.VescMin();
    auto RetHisto = grob::make_histo<double>(
                grob::mesh_grids(
                        grob::GridUniformHisto<double>(0,VescMax_nd,100),
                        grob::GridUniformHisto<double>(0,VescMax_nd,100)
                    )
                );
    auto G = [](){return rand()/(RAND_MAX+0.0);};
    grob::GridUniform<double> r_grid(0,1,BM["Vesc"].size());
    auto VescR = grob::make_function(r_grid,BM["Vesc"]);
    auto phiR = grob::make_function(r_grid,BM["phi"]);
    double Emin = -phiR.Values[0];
    auto LE_func = grob::make_function_f(
                    grob::GridUniform<double>(Emin ,0,1000),
                [&phiR](double e){return maxLndf(phiR,e);}
                );
    auto sc = sc_cnt(L_distrib,H_distrib,LE_func,phiR,G);
    sc.dens_scatter_function(RetHisto,
                             cmd_params.get<int>("Nmk",100000),
                             cmd_params.get<int>("Rcut",10));

    stools::PtreeSerializator<boost::property_tree::ptree> S{};
    auto H_to_Write = stools::Serialize(RetHisto,S);
    boost::property_tree::write_json(config_path_from(cmd_params.pgets("vsd","vsd.txt"),programm_path).string()
                                     ,H_to_Write);

    /*
    if(ptree_contain(cmd_params,"csv")){

    }*/



}
