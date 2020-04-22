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

#ifndef CGOGN_RENDERING_SHADERS_FLAT_COLOR_PER_FACE_H_
#define CGOGN_RENDERING_SHADERS_FLAT_COLOR_PER_FACE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(FlatColorPerFace, CGOGN_STR(FlatColorPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamFlatColorPerFace : public ShaderParam
{
	void set_uniforms() override;

public:
	VBO* vbo_pos_;
	VBO* vbo_color_;
	GLColor ambiant_color_;
	GLVec3 light_position_;
	bool double_side_;

	template <typename... Args>
	void fill(Args&&... args)
	{
		auto a = std::forward_as_tuple(args...);
		ambiant_color_ = std::get<0>(a);
		light_position_ = std::get<1>(a);
		double_side_ = std::get<2>(a);
	}

	using LocalShader = ShaderFlatColorPerFace;

	ShaderParamFlatColorPerFace(LocalShader* sh)
		: ShaderParam(sh), vbo_pos_(nullptr), ambiant_color_(0.05f, 0.05f, 0.05f, 1), light_position_(10, 100, 1000),
		  double_side_(true)
	{
	}

	inline ~ShaderParamFlatColorPerFace() override
	{
	}

	inline void set_vbos(const std::vector<VBO*>& vbos) override
	{
		vbo_pos_ = vbos[0];
		vbo_color_ = vbos[1];
		if (vbo_pos_ && vbo_color_)
			vao_initialized_ = true;
		else
			vao_initialized_ = false;
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_FLAT_COLOR_PER_FACE_H_
