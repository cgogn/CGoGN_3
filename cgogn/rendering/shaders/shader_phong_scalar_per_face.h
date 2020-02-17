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

#ifndef CGOGN_RENDERING_SHADERS_PHONG_SCALAR_PERFACE_H_
#define CGOGN_RENDERING_SHADERS_PHONG_SCALAR_PERFACE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{
DECLARE_SHADER_CLASS(PhongScalarPerFace,CGOGN_STR(PhongScalarPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamPhongScalarPerFace : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	VBO* vbo_pos_;
	VBO* vbo_norm_;
	VBO* vbo_scalar_; // perface
	GLVec3 light_position_;
	GLColor ambiant_color_;
	GLColor specular_color_;
	float32 specular_coef_;
	bool double_side_;

	using ShaderType = ShaderPhongScalarPerFace;

	ShaderParamPhongScalarPerFace(ShaderType* sh)
		: ShaderParam(sh),
		  vbo_pos_(nullptr),
		  vbo_norm_(nullptr),
		  vbo_scalar_(nullptr),
		  light_position_(), ambiant_color_(), specular_color_(),
		  specular_coef_(), double_side_()
	{
	}

	inline void set_vbos(const std::vector<VBO*>& vbos) override
	{
		vbo_pos_ = vbos[0];
		vbo_norm_ = vbos[1];
		vbo_scalar_ = vbos[2];
	}
};
} // namespace rendering
} // namespace cgogn
#endif
