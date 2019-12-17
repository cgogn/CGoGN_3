/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/rendering/types.h>
#define STB_IMAGE_IMPLEMENTATION
#include <cgogn/rendering/stb_image.h>
#include <fstream>

namespace cgogn
{

namespace rendering
{

// class Resource
//{
//	uint8* ptr_;
//	virtual uint8* code64_buffer();
//	virtual std::size_t code64_length();
// public:
//	void encode64(const std::string &src, const std::string& dest, const std::string& name);
//	void* get_buffer();
//	void release_buffer();
//	inline std::size_t size() {return (code64_length()*3)/4;};
//};

// void* Resource::get_buffer()
//{
//	std::size_t nb = code64_length();
//	uint8* dst = new uint8[(nb*3)/4];
//	uint8* src = code64_buffer();

//	while(nb>2)
//	{
//		*dst   |= (*src++)<<2;
//		*dst++ |= (*src)>>4;
//		*dst   |= (*src++)&15<<4;
//		*dst++ |= (*src)>>2;
//		*dst   |= (*src++)&3<<6;
//		*dst++ |= (*src++);
//		nb -=3;
//	}
//	if (nb>1)
//	{
//		*dst   |= (*src++)<<2;
//		*dst++ |= (*src)>>4;
//		--nb;
//	}
//	if (nb>0)
//	{
//		*dst   |= (*src++)&15<<4;
//		*dst++ |= (*src)>>2;
//	}
//	return dst;
//}

// void Resource::encode64(const std::string &src, const std::string& dest, const std::string& name)
//{
//	std::ifstream fin(src.c_str(), std::ios::in | std::ios_base::binary);
//	std::ofstream fout(dest.c_str(), std::ios::out);

//	const std::streamsize N = 3*1024*1024;
//	char* buff_in = new char[N];

//	std::streamsize size = 0;
//	std::streamsize nb = 0;
//	do
//	{
//		fout <<'\\"';
//		nb = fin.readsome(buff_in,N);
//		for (std::streamsize i=0; i<nb; ++i)
//		{
//			char c1 = ' '+(buff_in[i] >> 2);
//			char c2 = ' '+(((buff_in[i] & 3 )<<4) || (buff_in[i+1]>>4));
//			char c3 = ' '+(((buff_in[i+1] & 15 )<<2) || (buff_in[i+2]>>6));
//			char c4 = ' '+(buff_in[i+2] & 63);
//			fout <<  c1 <<c2 <<c3 <<c4;

//			if (i%64 ==63)
//				fout <<'\\"'<<std::endl;
//		}
//		size += nb;
//	} while (nb == N);

//}

//// 64 : 543210 543210 543210 543210
//// BIN: 765432 107654 321076 543210

GLImage::GLImage(int32 w, int32 h, int32 d) : width_(w), height_(h), bpp_(d), stb_(false)
{
	data_ = new uint8[d * h * w];
}

GLImage::GLImage(const std::string& filename) : stb_(true)
{
	data_ = stbi_load(filename.c_str(), &width_, &height_, &bpp_, 0);
}

GLImage::~GLImage()
{
	if (stb_)
		stbi_image_free(data_);
	else
		delete[] data_;
}

} // namespace rendering

} // namespace cgogn
