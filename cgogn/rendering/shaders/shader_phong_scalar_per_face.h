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
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/shaders/shader_function_color_maps.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(PhongScalarPerFace, true, CGOGN_STR(PhongScalarPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamPhongScalarPerFace : public ShaderParam
{
	void set_uniforms() override;

	std::array<VBO*, 3> vbos_;
	inline void set_texture_buffer_vbo(uint32 i, VBO* vbo) override
	{
		vbos_[i] = vbo;
	}
	void bind_texture_buffers() override;
	void release_texture_buffers() override;

	enum VBOName : uint32
	{
		VERTEX_POSITION = 0,
		VERTEX_NORMAL,
		FACE_SCALAR
	};

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	bool double_side_;
	GLColor specular_color_;
	float32 specular_coef_;
	shader_function::ColorMap::Uniforms color_map_;

	using ShaderType = ShaderPhongScalarPerFace;

	ShaderParamPhongScalarPerFace(ShaderType* sh)
		: ShaderParam(sh), ambiant_color_(0.05f, 0.05f, 0.05f, 1), light_position_(10, 100, 1000), double_side_(true),
		  specular_color_(1, 1, 1, 1), specular_coef_(250)
	{
		for (auto& v : vbos_)
			v = nullptr;
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_PHONG_SCALAR_PER_FACE_H_
