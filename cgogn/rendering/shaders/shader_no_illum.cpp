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

#include <iostream>

#include <cgogn/rendering/shaders/shader_no_illum.h>

namespace cgogn
{

namespace rendering
{

ShaderNoIllum* ShaderNoIllum::instance_ = nullptr;

ShaderNoIllum::ShaderNoIllum()
{
	char const* vertex_shader_source_ = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		in vec3 vertex_position;
		
		void main()
		{
			gl_Position = projection_matrix * model_view_matrix * vec4(vertex_position, 1.0);
		}
	)";

	char const* fragment_shader_source_ = R"(
		#version 330
		uniform bool double_side;
		uniform vec4 color;
		
		out vec4 frag_out;

		void main()
		{
			if (double_side || gl_FrontFacing)
				frag_out = color;
			else
				discard;
		}
	)";

	load2_bind(vertex_shader_source_, fragment_shader_source_, "vertex_position");
	add_uniforms("color", "double_side");
}

void ShaderParamNoIllum::set_uniforms()
{
	shader_->set_uniforms_values(color_, double_side_);
}

} // namespace rendering

} // namespace cgogn
