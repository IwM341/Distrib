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



template <typename L_E_Functype,typename Phi_FuncType,typename...HistoArgs>
auto r_distrib(grob::Histogramm<HistoArgs...> const&mHisto,
                            L_E_Functype const & L_E,
                            Phi_FuncType const &phi,
                           const size_t N_r = 1000,
               double R_cut = 10) noexcept{

        print("start calculations");
        progress_bar progress(mHisto.size(),100);


        auto Mi = mHisto.size();
        auto rd_func = grob::make_function_f(grob::GridUniform<double>(0,R_cut,N_r+1),
                              [](size_t i){return 0.0;});

        for(size_t i=0;i<Mi;++i)
        {
            auto MI = mHisto.Grid.FromLinear(i);

            auto EL_rect_L = mHisto.Grid[MI];
            auto E_L = EL_rect_L.template getR<0>().center();
            auto l_L = EL_rect_L.template getR<1>().center();
            auto L_L =l_L*L_E(E_L);

            double L00_L = EL_rect_L.template get<1,0b00>()*L_E(EL_rect_L.left());
            double L01_L = EL_rect_L.template get<1,0b10>()*L_E(EL_rect_L.left());
            double L10_L = EL_rect_L.template get<1,0b01>()*L_E(EL_rect_L.right());
            double L11_L = EL_rect_L.template get<1,0b11>()*L_E(EL_rect_L.right());

            //double T00_L = CalculatePeriod(phi,EL_rect_L.left(),L00_L,10);
            //double T01_L = CalculatePeriod(phi,EL_rect_L.left(),L01_L,10);
            //double T10_L = CalculatePeriod(phi,EL_rect_L.right(),L10_L,10);
            //double T11_L = CalculatePeriod(phi,EL_rect_L.right(),L11_L,10);

            TrajectoryPreInfo T_L = CalculatePeriod(phi,E_L,L_L,10);
            double TdEdL2 = EL_volume(EL_rect_L.left(),EL_rect_L.right(),
                                                              L00_L,L01_L,L10_L,L11_L,
                                                              T_L.T);
            //TrajectoryInfo TI = CalculateTrajectory(phi,);
            double rho_L = mHisto[MI]/TdEdL2;
            //double rho_L = mHisto[MI]/(4*M_PI*M_PI*T_L.T*(EL_rect_L.first().volume()*(L11_L+L01_L-L10_L-L00_L)));
            //COMPARE(TdEdL2,0.5*T_L.T*(EL_rect_L.first().volume()*(L11_L+L11_L-L00_L-L10_L)));
            const size_t rd_sz = rd_func.Grid.size();
            //double _E = EL_rect_L.first().center();
            //#pragma omp parallel for shared(rd_func,rho_L)
            //double msum = 0;
            //PVAR(EL_rect_L);
            for(size_t j=0;j<rd_sz;++j){
                double r = rd_func.Grid[j];
                double mphi = (r < 1 ? phi(r) : 1/r);
                if(mphi + EL_rect_L.right()>0){
                    double rge_dens = 4*M_PI*r*r;
                    double d3v_L= d_3_v_mes(EL_rect_L,L_E,phi,r);
                    //PVAR(_E+mphi);
                    //PVAR(EL_rect_L.inner().center()*L_E(_E)/r);
                    //PVAR(mphi);
                    //PVAR(_E);
                    //PVAR(L_E(_E));
                    //COMPARE(d3v_L_naive,d3v_L);
                    //msum += rge_dens*d3v_L*rho_L;
                    rd_func.Values[j]+=rge_dens*d3v_L*rho_L;
                }
            }
            //COMPARE(msum*rd_func.Grid[1],mHisto[MI]);
            progress++;

        }
        return rd_func;
}


int main(int argc, char ** argv){
    print("getting distribution in r");
    print("params: distrib.grid [filename of grid], \n",
          "distrib.values [filename of values]\n",
          "body [filename of body model]\n",
          "Rcut [default is 10]\n",
          "Nr [default is 1000]\n");
    boost::property_tree::ptree cmd_params;
    auto ret_key = parse_command_line_v1(argc,argv,cmd_params);
    if(!ret_key.empty()){
        std::cout << "error : " << ret_key<<std::endl;
        return 0;
    }
    boost::filesystem::path programm_path = ptree_gets(cmd_params,"config_path");

    std::string bm_filename;
    double Vesc1;
    try {
        bm_filename = cmd_params.get<std::string>("body");
        Vesc1 = cmd_params.get<double>("Vesc",1);

        if(ptree_gets(cmd_params,"debug") != ""){
            std::cout << "body filename: " << bm_filename<<std::endl;
            std::cout << "Vesc " << Vesc1<<std::endl;
        }
    } catch (std::exception &e) {
        std::cout << "can't get body model" <<std::endl;
        return 0;
    }
    bool is_debug = cmd_params.get<bool>("debug",false);
    if(is_debug){
        boost::property_tree::write_json(std::cout,cmd_params);
    }
    auto BM = BodyModel::fromFile(Vesc1,
                             config_path_from(bm_filename,programm_path).string());

    //auto grid_fname = cmd_params.pgets("grid");
    //auto values_fname =cmd_params.pgets("values");
    auto _distrib = loadHisto<grob::MultiGridHisto<
                                    grob::GridVectorHisto<double>,
                                    std::vector<grob::GridVectorHisto<double>>
                                >,
                                double>
            (cmd_params.get_child("distrib"),programm_path);

    double VescMax_nd = BM.VescMax()/BM.VescMin();

    grob::GridUniform<double> r_grid(0,1,BM["Vesc"].size());
    auto VescR = grob::make_function(r_grid,BM["Vesc"]);
    auto phiR = grob::make_function(r_grid,BM["phi"]);
    double Emin = -phiR.Values[0];
    auto LE_func = grob::make_function_f(
                    grob::GridUniform<double>(Emin ,0,1000),
                [&phiR](double e){return maxLndf(phiR,e);}
                );
    double Rcut = cmd_params.get<double>("Rcut",10);
    size_t Nr = cmd_params.get<int>("Nr",1000);
    if(is_debug){
        PVAR(vector_sum(_distrib.Values));
    }
    auto rd_func = r_distrib(_distrib,LE_func,phiR,Nr,Rcut);
    if(is_debug){
        PVAR(rd_func.Grid[1]*vector_sum(rd_func.Values));
    }
    grob::as_csv(rd_func).save(std::ofstream(cmd_params.pgets("out")),6,std::defaultfloat);

    /*
    if(ptree_contain(cmd_params,"csv")){

    }*/



}
