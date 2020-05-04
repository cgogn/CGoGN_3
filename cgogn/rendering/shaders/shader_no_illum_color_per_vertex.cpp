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

#include <cgogn/rendering/shaders/shader_no_illum_color_per_vertex.h>

namespace cgogn
{

namespace rendering
{

static char const* vertex_shader_source =
	R"(#version 330
	in vec3 vertex_pos;
	in vec3 vertex_color;
	uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;
	out vec3 color;
	void main()
	{
		gl_Position = projection_matrix * model_view_matrix * vec4(vertex_pos,1.0);
		color = vertex_color;
	}
	)";

static char const* fragment_shader_source =
	R"(#version 330
	out vec3 fragColor;
	in vec3 color;
	uniform bool double_side;

	void main()
	{
		if (gl_FrontFacing || double_side)
			fragColor = color;
		else
			discard;
	}
	)";

ShaderNoIllumColorPerVertex* ShaderNoIllumColorPerVertex::instance_ = nullptr;

ShaderNoIllumColorPerVertex::ShaderNoIllumColorPerVertex()
{
	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos", "vertex_color");
	add_uniforms("double_side");
}

void ShaderParamNoIllumColorPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(double_side_);
}

} // namespace rendering

} // namespace cgogn
