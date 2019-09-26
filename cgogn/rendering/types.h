/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_RENDERING_TYPES_H_
#define CGOGN_RENDERING_TYPES_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/core/utils/numerics.h>

// #include <Eigen/Core>
#include <Eigen/Dense>
// #include <Eigen/Eigen>
// #include <Eigen/Geometry>
// #include <Eigen/SVD>

#include <string>

namespace cgogn
{

namespace rendering
{

using GLVec2 = Eigen::Vector2f;
using GLVec3 = Eigen::Vector3f;
using GLVec4 = Eigen::Vector4f;

using GLVec2d = Eigen::Vector2d;
using GLVec3d = Eigen::Vector3d;
using GLVec4d = Eigen::Vector4d;

using GLMat3 = Eigen::Matrix3f;
using GLMat4 = Eigen::Matrix4f;

using GLMat3d = Eigen::Matrix3d;
using GLMat4d = Eigen::Matrix4d;

using GLColor = Eigen::Vector4f;

using Transfo3d = Eigen::Affine3d;

// inline GLColor col4i(uint8 R, uint8 G, uint8 B, uint8 A)
// {
// 	return GLColor(R/255.0f, G/255.0f, B/255.0f, A/255.0f);
// }

// inline GLColor col3i(uint8 R, uint8 G, uint8 B)
// {
// 	return GLColor(R/255.0f, G/255.0f, B/255.0f, 1.0f);
// }

// inline GLColor col1i(uint8 R)
// {
// 	return GLColor(R/255.0f, R/255.0f, R/255.0f, 1.0f);
// }

class CGOGN_RENDERING_EXPORT GLImage
{
	uint8* data_;
	int32 width_;
	int32 height_;
	int32 bpp_; // 1:grey, 3 RGB, 4 RGBA
	bool stb_;

public:

	GLImage(int32 w, int32 h, int32 d);
	GLImage(const std::string& filename);
	~GLImage();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(GLImage);

	inline int32 width() const { return width_; }
	inline int32 height() const { return height_; }
	inline int32 depth() const { return bpp_; }
	inline const uint8* data() const { return data_; }
	inline void set_pixel(int32 x, int32 y, const GLColor& col)
	{
		uint8* ptr = data_+ (bpp_ * (y * width_ + x));
		for (int32 i = 0; i < bpp_; ++i)
			*ptr++ = uint8(255 * col[i]);
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_TYPES_H_
