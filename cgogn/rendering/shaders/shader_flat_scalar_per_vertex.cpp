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

#include <cgogn/rendering/shaders/shader_flat_scalar_per_vertex.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatScalarPerVertex* ShaderFlatScalarPerVertex::instance_ = nullptr;

ShaderFlatScalarPerVertex::ShaderFlatScalarPerVertex()
{
	const char* vertex_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		in vec3 vertex_position;
		in float vertex_scalar;

		out vec3 position;
		out float value;
		out vec3 color;

		//_insert_colormap_function_here

		void main()
		{
			vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
			position = position4.xyz;
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
		uniform bool show_iso_lines;
		uniform int nb_iso_lines;
		
		in vec3 position;
		in float value;
		in vec3 color;
		
		out vec4 frag_out;
		
		void main()
		{
			if (!gl_FrontFacing && !double_side)
				discard;
			else
			{
				vec3 N = normalize(cross(dFdx(position), dFdy(position)));
				vec3 L = normalize(light_position - position);
				float lambert = dot(N, L);
				vec3 frag_color = ambiant_color.rgb + lambert * color;
				frag_out = vec4(frag_color, 1.0);
				if (show_iso_lines)
				{
					float s = value * float(nb_iso_lines);
					if (s - floor(s) < 0.05)
						frag_out = vec4(0.0);
				}
			}
		}
	)";

	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_function_here"), shader_function::ColorMap::source);
	load2_bind(v_src, fragment_shader_source, "vertex_position", "vertex_scalar");
	get_uniforms("ambiant_color", "light_position", "double_side", "show_iso_lines", "nb_iso_lines",
				 shader_function::ColorMap::uniform_names[0], shader_function::ColorMap::uniform_names[1],
				 shader_function::ColorMap::uniform_names[2], shader_function::ColorMap::uniform_names[3]);
}

void ShaderParamFlatScalarPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(ambiant_color_, light_position_, double_side_, show_iso_lines_, nb_iso_lines_,
								 color_map_.color_map_, color_map_.expansion_, color_map_.min_value_,
								 color_map_.max_value_);
}

} // namespace rendering

} // namespace cgogn
