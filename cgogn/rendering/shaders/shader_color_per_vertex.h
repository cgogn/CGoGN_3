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

#ifndef CGOGN_RENDERING_SHADER_COLOR_PER_VERTEX_H_
#define CGOGN_RENDERING_SHADER_COLOR_PER_VERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(FlatColorPerVertex, CGOGN_STR(FlatColorPerVertex))

class CGOGN_RENDERING_EXPORT ShaderParamFlatColorPerVertex : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(light_position_);
	}

public:
	GLVec3 light_position_;

	using LocalShader = ShaderFlatColorPerVertex;

	ShaderParamFlatColorPerVertex(LocalShader* sh) : ShaderParam(sh), light_position_(10, 100, 1000)
	{
	}

	inline ~ShaderParamFlatColorPerVertex() override
	{
	}
};

DECLARE_SHADER_CLASS(PhongColorPerVertex, CGOGN_STR(PhongColorPerVertex))

class CGOGN_RENDERING_EXPORT ShaderParamPhongColorPerVertex : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(light_position_);
	}

public:
	GLVec3 light_position_;

	using LocalShader = ShaderPhongColorPerVertex;

	ShaderParamPhongColorPerVertex(LocalShader* sh) : ShaderParam(sh), light_position_(10, 100, 1000)
	{
	}

	inline ~ShaderParamPhongColorPerVertex() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_COLOR_PER_VERTEX_H_
