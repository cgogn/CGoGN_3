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

#include <cgogn/rendering/gl_image.h>

#define STB_IMAGE_IMPLEMENTATION
#include <cgogn/rendering/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <cgogn/rendering/stb_image_write.h>

#include <cgogn/core/utils/string.h>

namespace cgogn
{

namespace rendering
{

GLImage::GLImage(int32 width, int32 height, int32 channels)
	: width_(width), height_(height), channels_(channels), stb_loaded_(false)
{
	data_ = new uint8[channels_ * height_ * width_];
}

GLImage::GLImage(const std::string& filename) : stb_loaded_(true)
{
	data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_, 0);
}

GLImage::~GLImage()
{
	if (stb_loaded_)
		stbi_image_free(data_);
	else
		delete[] data_;
}

void GLImage::save(const std::string& filename, bool flip_y) const
{
	stbi_flip_vertically_on_write(flip_y);
	std::string ext = extension(filename);
	if (ext == "jpg")
		stbi_write_jpg(filename.c_str(), width_, height_, channels_, data_, 90);
	else if (ext == "png")
		stbi_write_png(filename.c_str(), width_, height_, channels_, data_, width_ * channels_);
	else if (ext == "bmp")
		stbi_write_bmp(filename.c_str(), width_, height_, channels_, data_);
}

} // namespace rendering

} // namespace cgogn
