#include <iostream>
#include <fstream>
#include <grob/multigrid.hpp>
#include <grob/serialization.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

template <typename T>
void get_slice(size_t i,std::istream &ifs,std::ostream & ofs){
	size_t N,M;
	
	
	ifs.read((char*)&N,sizeof(N));
	ifs.read((char*)&M,sizeof(M));

	
	size_t N1 = 1;

	ofs.write((char*)&N1,sizeof(N1));
	ofs.write((char*)&M,sizeof(M));

	std:: cout << "N = " << N << ", M = " << M << ", i = " << i << std::endl;
	if(i>= N){
		throw std::runtime_error("too big i" );
	}

	size_t offset = M*i*sizeof(T);
	ifs.ignore(offset);
	T tmp;
	for(size_t j=0;j<M;++j){
		ifs.read((char*)&tmp,sizeof(T));
		ofs.write((char*)&tmp,sizeof(T));
	}

}

int main(int argc, char ** argv){
	if(argc < 3){
		std::cout << "need arguments: [input file] [index] or\n";
		std::cout << "[input file] [E] [l]" << std::endl;
		return 0;
	}
	boost::property_tree::ptree P;
	boost::property_tree::read_json(argv[1],P);
	stools::PtreeSerializator<decltype(P)> S{};
	auto mGrid = stools::DeSerialize<
		grob::MultiGridHisto<
			grob::GridVectorHisto<double>,
			std::vector<grob::GridVectorHisto<double>>
		>
	>(P,S);

	if(argc == 3){
		auto Index=  mGrid.FromLinear(std::stoi(argv[2]));
		std::cout << "Index = " << Index << ", Rect = " << mGrid[Index] << std::endl;
	}else{
		auto [Index,Rect] = mGrid.FindElement(
			std::stod(argv[2]),
			std::stod(argv[3])
		);

		std::cout << "Index = " << Index << ", Rect = " << Rect << std::endl;
		std::cout << "LinearIndex = " << mGrid.LinearIndex(Index) << std::endl;
	}
	return 0;
}