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

/**
 * @brief  \f$\int{dx\sqrt{Ax+B}}\f$
 * @param A  A
 * @param B  B
 * @return
 */
double sqrt_integral(double A,double B,double x0,double x1){
    if(A == 0){
        return (x1-x0)*sqrt(B);
    }
    double F02 = A*x0+B;
    double F0 = (F02 > 0 ? sqrt(F02) : 0);
    double F12 = A*x1+B;
    double F1 = (F12 > 0 ? sqrt(F12) : 0);
    return (2.0/3)*(F1*F12-F0*F02)/(A);
}

/**
 * @brief d_3_v_mes \f$\int{d^3v = \cfrac2{\pi vdv dL^2}{r\sqrt{r^2v^2-L^2}}}\f$
 * @param R - Rect, contaning E and L
 * @param L_E - function to get L(E)
 * @param r - radius, place, where to calculate
 */
template <typename RectType,typename L_E_Functype,typename PhiFunctype>
inline double d_3_v_mes(RectType const &R,L_E_Functype const & L_E,
                 PhiFunctype const& phi,double r) noexcept{
    if(r == 0)
        return 0;
    double L0 = L_E(R.left());
    double L1 = L_E(R.right());
    auto [l0,l1] = R.template getR<1>();
    double dL2_dE = (L1*L1-L0*L0)*(R.right()-R.left());
    double r2 = r*r;
    return M_PI*(
                sqrt_integral(r2-l1*l1*dL2_dE,
                              r2*phi(r)-l1*l1*(L0*L0-dL2_dE*R.left()),
                              R.left(),R.right()
                              )-
                sqrt_integral(r2-l0*l0*dL2_dE,
                                              r2*phi(r)-l0*l0*(L0*L0-dL2_dE*R.left()),
                                              R.left(),R.right()
                                              )
                          )/r;



}



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

                double T00_L = CalculatePeriod(phi,EL_rect_L.left(),L00_L,10);
                double T01_L = CalculatePeriod(phi,EL_rect_L.left(),L01_L,10);
                double T10_L = CalculatePeriod(phi,EL_rect_L.right(),L10_L,10);
                double T11_L = CalculatePeriod(phi,EL_rect_L.right(),L11_L,10);

                double rho_L = 1/
                        EL_bin_volume(EL_rect_L.left(),EL_rect_L.right(),
                                      L00_L,L01_L,L10_L,L11_L,
                                      T00_L,T01_L,T10_L,T11_L);


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

                    double T00_H = CalculatePeriod(phi,EL_rect_H.left(),L00_H,10);
                    double T01_H = CalculatePeriod(phi,EL_rect_H.left(),L01_H,10);
                    double T10_H = CalculatePeriod(phi,EL_rect_H.right(),L10_H,10);
                    double T11_H = CalculatePeriod(phi,EL_rect_H.right(),L11_H,10);

                    double rho_H = 1/
                            EL_bin_volume(EL_rect_H.left(),EL_rect_H.right(),
                                          L00_H,L01_H,L10_H,L11_H,
                                          T00_H,T01_H,T10_H,T11_H);


                    TrajectoryPreInfo T_L = CalculatePeriod(phi,E_L,L_L,10);
                    TrajectoryPreInfo T_H = CalculatePeriod(phi,E_H,L_H,10);

                    double & HL_Value = ANN_HL[i*NH+j];

                    double r_min = std::max(T_L.rmin,T_H.rmin);
                    double r_max = std::min(T_L.rmax,T_H.rmax);
                    if(r_min < r_max){
                        for(size_t i=0;i<N_mk;++i){


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
    auto Ann_Matrix = sc.annihilation_matrix([](auto){return 1.0;},cmd_params.get<int>("Nmk",100000),
                             cmd_params.get<int>("Rcut",10));


    saveMatrix(Ann_Matrix.data(),L_Grid.size(),H_Grid.size(),
               config_path_from(cmd_params.pgets("ann","ann.bmat"),programm_path).string(),MatrixFormat::BINARY);




}
