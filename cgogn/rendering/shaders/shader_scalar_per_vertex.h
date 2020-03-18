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

#ifndef CGOGN_RENDERING_SHADER_SCALAR_PER_VERTEX_H_
#define CGOGN_RENDERING_SHADER_SCALAR_PER_VERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_function_color_maps.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(ScalarPerVertex, CGOGN_STR(ScalarPerVertex))

class CGOGN_RENDERING_EXPORT ShaderParamScalarPerVertex : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(color_map_, expansion_, min_value_, max_value_, show_iso_lines_, nb_iso_levels_,
									 light_position_);
	}

public:
	ColorMap color_map_;
	int32 expansion_;
	float32 min_value_;
	float32 max_value_;
	bool show_iso_lines_;
	int32 nb_iso_levels_;
	GLVec3 light_position_;

	using LocalShader = ShaderScalarPerVertex;

	ShaderParamScalarPerVertex(LocalShader* sh)
		: ShaderParam(sh), color_map_(BWR), expansion_(0), min_value_(.0f), max_value_(1.0f), show_iso_lines_(false),
		  nb_iso_levels_(10), light_position_(10, 100, 1000)
	{
	}

	inline ~ShaderParamScalarPerVertex() override
	{
	}
};

DECLARE_SHADER_CLASS(ScalarPerVertexGouraud, CGOGN_STR(ScalarPerVertexGouraud))

class CGOGN_RENDERING_EXPORT ShaderParamScalarPerVertexGouraud : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(color_map_, expansion_, min_value_, max_value_, show_iso_lines_, nb_iso_levels_,
									 light_position_);
	}

public:
	ColorMap color_map_;
	int32 expansion_;
	float32 min_value_;
	float32 max_value_;
	bool show_iso_lines_;
	int32 nb_iso_levels_;
	GLVec3 light_position_;

	using LocalShader = ShaderScalarPerVertexGouraud;

	ShaderParamScalarPerVertexGouraud(LocalShader* sh)
		: ShaderParam(sh), color_map_(BWR), expansion_(0), min_value_(.0f), max_value_(1.0f), show_iso_lines_(false),
		  nb_iso_levels_(10), light_position_(10, 100, 1000)
	{
	}

	inline ~ShaderParamScalarPerVertexGouraud() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_SCALAR_PER_VERTEX_H_
