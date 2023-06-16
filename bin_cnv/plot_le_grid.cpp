#include <iostream>
#include <fstream>
#include <grob/multigrid.hpp>
#include <grob/serialization.hpp>
#include <grob/grid_objects.hpp>
#include <grob/csv_io.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#define DBG
#undef  DBG
int main(int argc, char ** _argv){
	#ifdef DBG
	std::vector<const char*> argv = {"progname","D:/Important/work/Distrib/tests/100_dm/capt_light/l_grid_100_dm.txt",
				"D:/Important/work/Distrib/tests/100_dm/capt_light/LE.json.txt","100","100","1"};
	argc = 6;
	#else
	char ** argv = _argv;
	#endif
	if(argc < 6){
		std::cout << "need arguments: [grid] [LE] [W] [H] [dx]"<< std::endl;
		return 0;
	}
	boost::property_tree::ptree P;
	boost::property_tree::read_json(argv[1],P);

	boost::property_tree::ptree P_LE_f;
	boost::property_tree::read_json(argv[2],P_LE_f);

	stools::PtreeSerializator<decltype(P)> S{};
	auto mGrid = stools::DeSerialize<
		grob::MultiGridHisto<
			grob::GridVectorHisto<double>,
			std::vector<grob::GridVectorHisto<double>>
		>
	>(P,S);

	auto LE_func = stools::DeSerialize<grob::GridFunction<
				grob::linear_interpolator,
				grob::GridUniform<double>,
				std::vector<double>
			>
		>(P_LE_f,S);
	size_t W = std::stoi(argv[3]);
	size_t H = std::stoi(argv[4]);
	size_t dx = std::stoi(argv[5]);
	grob::GridUniform<double> Egrid(mGrid.first().front(),mGrid.first().back(),W+1);

	std::cout << "Egrid: " << Egrid << std::endl; 

	double dE = (Egrid.back()-Egrid.front());

	grob::GridUniform<double> L1_grid(0,mGrid.inner().back().back(),H+1);
	
	std::cout << "Egrid: " << L1_grid << std::endl;
	double dl = L1_grid.back()-L1_grid.front();

	grob::GridUniform<double> L_grid(0,L1_grid.back()*LE_func.Values.back(),H+1);
	std::cout << "Egrid: " << L_grid << std::endl;
	double dL = L_grid.back()-L_grid.front();

	auto PreGridFunc = grob::make_function_f(
		grob::mesh_grids(Egrid,L1_grid),
		[&](double E,double l){
			auto R = mGrid[mGrid.pos(E,l)];
			double de = std::min(
				std::abs(E-R.left()),
				std::abs(E-R.right())
			);
			double dl_ = std::min(
				std::abs(l-R.inner().left()),
				std::abs(l-R.inner().right())
			);
			if(dl_/dl*H*2 < dx || de/dE*W*2 < dx){
				return 1;
			}else {
				return 0;
			}
		}
	);
	grob::as_csv(PreGridFunc).save(std::ofstream(argv[1] + std::string("-pre.txt")),6,std::defaultfloat);
	auto GridFunc = grob::make_function_f(
		grob::mesh_grids(Egrid,L_grid),
		[&](double E,double L){
			double Lmax = LE_func(E);
			auto R = mGrid[mGrid.pos(E,L/Lmax)];
			double de = std::min(
				std::abs(E-R.left()),
				std::abs(E-R.right())
			);
			double dL_ = std::min(
				std::abs(L-R.inner().left()*Lmax),
				std::abs(L-R.inner().right()*Lmax)
			);
			if(L > Lmax){
				if(dL_/dL*H*2 < dx){
					return 1;
				}
				else{
					return 0;
				}
			}
			if(dL_/dL*H*2 < dx || de/dE*W*2 < dx){
				return 1;
			}else {
				return 0;
			}
		}
	);
	grob::as_csv(GridFunc).save(std::ofstream(argv[1] + std::string(".txt")),6,std::defaultfloat);

	return 0;
}