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

#ifndef CGOGN_RENDERING_SHADERS_NO_ILLUM_SCALAR_PER_FACE_H_
#define CGOGN_RENDERING_SHADERS_NO_ILLUM_SCALAR_PER_FACE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/shaders/shader_function_color_maps.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(NoIllumScalarPerFace, true, CGOGN_STR(NoIllumScalarPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamNoIllumScalarPerFace : public ShaderParam
{
	void set_uniforms() override;

	std::array<VBO*, 2> vbos_;
	inline void set_texture_buffer_vbo(uint32 i, VBO* vbo) override
	{
		vbos_[i] = vbo;
	}
	void bind_texture_buffers() override;
	void release_texture_buffers() override;

	enum VBOName : uint32
	{
		VERTEX_POSITION = 0,
		FACE_SCALAR
	};

public:
	bool double_side_;
	shader_function::ColorMap::Uniforms color_map_;

	using ShaderType = ShaderNoIllumScalarPerFace;

	ShaderParamNoIllumScalarPerFace(ShaderType* sh) : ShaderParam(sh)
	{
		for (auto& v : vbos_)
			v = nullptr;
	}

	inline ~ShaderParamNoIllumScalarPerFace() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_NO_ILLUM_SCALAR_PER_FACE_H_
