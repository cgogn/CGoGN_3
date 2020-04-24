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

#ifndef CGOGN_RENDERING_SHADERS_PHONG_COLOR_PER_FACE_H_
#define CGOGN_RENDERING_SHADERS_PHONG_COLOR_PER_FACE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{
DECLARE_SHADER_CLASS(PhongColorPerFace, true, CGOGN_STR(PhongColorPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamPhongColorPerFace : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	std::array<VBO*, 3> vbos_;
	GLColor ambiant_color_;
	GLColor specular_color_;
	float32 specular_coef_;
	GLVec3 light_position_;
	bool double_side_;

	inline void pick_parameters(const PossibleParameters& pp) override
	{
		ambiant_color_ = pp.ambiant_color_;
		specular_color_ = pp.specular_color_;
		specular_coef_ = pp.specular_coef_;
		light_position_ = pp.light_position_;
		double_side_ = pp.double_side_;
	}
	using ShaderType = ShaderPhongColorPerFace;

	ShaderParamPhongColorPerFace(ShaderType* sh) : ShaderParam(sh)
	{
		for (auto& v : vbos_)
			v = nullptr;
	}

	inline VBO** vbo_tb(uint32 i) override
	{
		return &vbos_[i];
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_PHONG_COLOR_PER_FACE_H_
