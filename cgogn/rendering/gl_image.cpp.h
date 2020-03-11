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

namespace cgogn
{
namespace rendering
{

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
