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

#define CGOGN_RENDER_SHADERS_PHONG_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_phong_scalar_per_vertex.h>

namespace cgogn
{

namespace rendering
{

ShaderPhongScalarPerVertex* ShaderPhongScalarPerVertex::instance_ = nullptr;

ShaderPhongScalarPerVertex::ShaderPhongScalarPerVertex()
{
	const char* vertex_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform mat3 normal_matrix;

		in vec3 vertex_position;
		in vec3 vertex_normal;
		in float vertex_scalar;

		out vec3 position;
		out vec3 normal;
		out float value;
		out vec3 color;

		//_insert_colormap_function_here

		void main()
		{
			vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
			position = position4.xyz;
			normal = normal_matrix * vertex_normal;
			value = transform_value(vertex_scalar);
			color = value2color(value);
			gl_Position = projection_matrix * position4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform vec4 ambiant_color;
		uniform vec3 light_position;
		uniform bool double_side;
		uniform vec4 specular_color;
		uniform float specular_coef;
		uniform bool show_iso_lines;
		uniform int nb_iso_lines;

		in vec3 position;
		in vec3 normal;
		in float value;
		in vec3 color;
		
		out vec4 frag_out;
		
		void main()
		{
			vec3 N = normalize(normal);
			if (!gl_FrontFacing)
			{
				if (!double_side)
					discard;
				N *= -1.0;
			}
			vec3 L = normalize(light_position - position);
			vec4 final_color = ambiant_color;
			float lambert_term = clamp(dot(N, L), 0.0, 1.0);
			final_color += vec4(color, 1.0) * lambert_term;
			vec3 E = normalize(-position);
			vec3 R = reflect(-L, N);
			float specular = pow(max(dot(R, E), 0.0), specular_coef);
			final_color += specular_color * specular;
			frag_out = vec4(final_color.rgb, 1.0);
			if (show_iso_lines)
			{
				float s = value * float(nb_iso_lines);
				if (s - floor(s) < 0.05)
					frag_out = vec4(0.0);
			}
		}
	)";

	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_function_here"), shader_function::ColorMap::source);
	load2_bind(v_src, fragment_shader_source, "vertex_position", "vertex_normal", "vertex_scalar");
	get_uniforms("light_position", "ambiant_color", "specular_color", "specular_coef", "double_side", "show_iso_lines",
				 "nb_iso_lines", shader_function::ColorMap::uniform_names[0],
				 shader_function::ColorMap::uniform_names[1], shader_function::ColorMap::uniform_names[2],
				 shader_function::ColorMap::uniform_names[3]);
}

void ShaderParamPhongScalarPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(light_position_, ambiant_color_, specular_color_, specular_coef_, double_side_,
								 show_iso_lines_, nb_iso_lines_, color_map_.color_map_, color_map_.expansion_,
								 color_map_.min_value_, color_map_.max_value_);
}

} // namespace rendering

} // namespace cgogn
