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

#include <cgogn/rendering/shaders/shader_phong.h>
#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderPhong* ShaderPhong::instance_ = nullptr;

ShaderPhong::ShaderPhong()
{
	const char* vertex_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform mat3 normal_matrix;

		in vec3 vertex_position;
		in vec3 vertex_normal;
		
		out vec3 position;
		out vec3 normal;
		
		void main ()
		{
			vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
			position = position4.xyz;
			normal = normal_matrix * vertex_normal;
			gl_Position = projection_matrix * position4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform vec4 front_color;
		uniform vec4 back_color;
		uniform vec4 ambiant_color;
		uniform vec3 light_position;
		uniform bool double_side;
		uniform vec4 specular_color;
		uniform float specular_coef;
		uniform bool ghost_mode;

		in vec3 position;
		in vec3 normal;
		
		out vec4 frag_out;
		
		void main()
		{
			vec3 N = normalize(normal);
			vec3 L = normalize(light_position - position);
			vec4 final_color = ambiant_color;
			vec4 current_color = front_color;
			if (!gl_FrontFacing)
			{
				if (!double_side)
					discard;
				N *= -1.0;
				current_color = back_color;
			}
			float lambert_term = clamp(dot(N, L), 0.0, 1.0);
			if (ghost_mode)
				lambert_term = 0.4 * pow(1.0 - lambert_term, 2);
			final_color += current_color * lambert_term;
			vec3 E = normalize(-position);
			vec3 R = reflect(-L, N);
			float specular = pow(max(dot(R, E), 0.0), specular_coef);
			final_color += specular_color * specular;
			frag_out = vec4(final_color.rgb, 1.0);
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position", "vertex_normal");
	get_uniforms("light_position", "front_color", "back_color", "ambiant_color", "specular_color", "specular_coef",
				 "double_side", "ghost_mode");
}

void ShaderParamPhong::set_uniforms()
{
	shader_->set_uniforms_values(light_position_, front_color_, back_color_, ambiant_color_, specular_color_,
								 specular_coef_, double_side_, ghost_mode_);
}

} // namespace rendering

} // namespace cgogn
