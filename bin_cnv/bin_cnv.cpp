#include <iostream>
#include <fstream>

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
	if(argc < 4){
		std::cout << "need arguments: [input file] [slice index] [output file]";
		return 0;
	}
	std::ifstream ifs(argv[1],std::ios::binary);
	if(!ifs.is_open()){
		throw std::runtime_error("no input file");
	}
	size_t i  = std::stoi(argv[2]);
	std::ofstream ofs(argv[3],std::ios::binary);
	if(!ofs.is_open()){
		throw std::runtime_error("no output file");
	}
	get_slice<double>(i,ifs,ofs);
	return 0;
}