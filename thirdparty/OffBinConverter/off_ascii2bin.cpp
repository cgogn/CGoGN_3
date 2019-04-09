#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string>

#include <cgogn/core/utils/endian.h>
#include <cgogn/core/utils/numerics.h>

using namespace cgogn::numerics;

int main(int argc, char **argv)
{
	if (argc < 3)
	{
		std::cerr << argv[0] << "input_ascii.off output_bin.off" << std::endl;
		return 1;
	}

	std::ifstream ifs(argv[1],std::ios::in);
	std::ofstream ofs(argv[2],std::ios::out|std::ofstream::binary);

	std::string str;
	uint32 nv;
	uint32 np;
	uint32 ne;

	ifs >> str;

	if (str != "OFF")
	{
		std::cerr << "no OFF comment at begin"<< std::endl;
		return 1;
	}

	ifs >> nv;
	ifs >> np;
	ifs >> ne;

	uint32 nv_be = cgogn::swap_endianness_native_big(nv);
	uint32 np_be = cgogn::swap_endianness_native_big(np);
	uint32 ne_be = cgogn::swap_endianness_native_big(ne);


	ofs << "OFF BINARY"<< std::endl;
	ofs.write(reinterpret_cast<char*>(&nv_be),4);
	ofs.write(reinterpret_cast<char*>(&np_be),4);
	ofs.write(reinterpret_cast<char*>(&ne_be),4);

	float* vertices = new float[3*nv];

	for (uint32 i=0; i<nv; ++i)
	{
		ifs >> vertices[3*i]  >> vertices[3*i+1] >> vertices[3*i+2];
	}

	uint32* ptr = reinterpret_cast<uint32*>(vertices);
	for (uint32 i=0; i<3*nv;++i)
	{
		*ptr = cgogn::swap_endianness_native_big(*ptr);
		ptr++;
	}


	ofs.write(reinterpret_cast<char*>(vertices),4*3*nv);

	delete[] vertices;

	std::vector<uint32> prim;
	prim.reserve(8*1024*1024);

	for (uint32 i=0; i<np; ++i)
	{
		uint32 nb;
		ifs >> nb;
		prim.push_back(nb);
		for (uint32 j=0; j<nb; ++j)
		{
			uint32 ind;
			ifs >> ind;
			prim.push_back(ind);
		}
	}

	ptr = reinterpret_cast<uint32*>(&(prim[0]));
	for (uint32 i=0; i<prim.size();++i)
	{
		*ptr = cgogn::swap_endianness_native_big(*ptr);
		ptr++;
	}

	ofs.write(reinterpret_cast<char*>(&(prim[0])),std::streamsize(4*prim.size()));

	ofs.close();
	ifs.close();

	return 0;
}
