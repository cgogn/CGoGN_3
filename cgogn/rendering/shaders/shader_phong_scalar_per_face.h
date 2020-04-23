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

#ifndef CGOGN_RENDERING_SHADERS_PHONG_SCALAR_PER_FACE_H_
#define CGOGN_RENDERING_SHADERS_PHONG_SCALAR_PER_FACE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_function_color_maps.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(PhongScalarPerFace, CGOGN_STR(PhongScalarPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamPhongScalarPerFace : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	VBO* vbo_pos_;
	VBO* vbo_norm_;
	VBO* vbo_scalar_;
	ColorMap color_map_;
	int32 expansion_;
	float32 min_value_;
	float32 max_value_;
	GLColor ambiant_color_;
	GLColor specular_color_;
	float32 specular_coef_;
	GLVec3 light_position_;
	bool double_side_;

	template <typename... Args>
	void fill(Args&&... args)
	{
		auto a = std::forward_as_tuple(args...);
		ambiant_color_ = std::get<0>(a);
		specular_color_ = std::get<1>(a);
		specular_coef_ = std::get<2>(a);
		light_position_ = std::get<3>(a);
		double_side_ = std::get<4>(a);
	}

	using ShaderType = ShaderPhongScalarPerFace;

	ShaderParamPhongScalarPerFace(ShaderType* sh)
		: ShaderParam(sh), vbo_pos_(nullptr), vbo_norm_(nullptr), vbo_scalar_(nullptr), color_map_(BWR), expansion_(0),
		  min_value_(.0f), max_value_(1.0f), ambiant_color_(color_ambiant_default), specular_color_(1, 1, 1, 1),
		  specular_coef_(250), light_position_(10, 100, 1000), double_side_(true)
	{
	}

	inline void set_vbos(const std::vector<VBO*>& vbos) override
	{
		vbo_pos_ = vbos[0];
		vbo_norm_ = vbos[1];
		vbo_scalar_ = vbos[2];
		if (vbo_pos_ && vbo_norm_ && vbo_scalar_)
			vao_initialized_ = true;
		else
			vao_initialized_ = false;
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_PHONG_SCALAR_PER_FACE_H_
