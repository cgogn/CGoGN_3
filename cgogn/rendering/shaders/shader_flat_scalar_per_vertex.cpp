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

static const char* vertex_shader_source =
	R"(#version 150
	in vec3 vertex_pos;
	in float vertex_scalar;

	uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;

	out vec3 pos;
	out vec3 color;

	//_insert_colormap_funcion_here

	void main()
	{
		vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);
		color = scalar2color(vertex_scalar);
		pos = pos4.xyz;
		gl_Position = projection_matrix * pos4;
	}
	)";

static const char* fragment_shader_source =
	R"(#version 150
	out vec3 fragColor;
	uniform vec4 ambiant_color;
	uniform vec3 light_position;
	uniform bool double_side;
	in vec3 pos;
	in vec3 color;
	void main()
	{
		vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));
		vec3 L = normalize(light_position-pos);
		float lambert = dot(N,L);
		if (!gl_FrontFacing && !double_side)
			 discard;
		else
			fragColor = ambiant_color.rgb + lambert*color.rgb;
	}
	)";

ShaderFlatScalarPerVertex* ShaderFlatScalarPerVertex::instance_ = nullptr;

ShaderFlatScalarPerVertex::ShaderFlatScalarPerVertex()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"), shader_funcion::ColorMap::source);
	load2_bind(v_src, fragment_shader_source, "vertex_pos", "vertex_scalar");

	add_uniforms("ambiant_color", "light_position", "double_side", shader_funcion::ColorMap::name[0],
				 shader_funcion::ColorMap::name[1], shader_funcion::ColorMap::name[2],
				 shader_funcion::ColorMap::name[3]);
}

void ShaderParamFlatScalarPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(ambiant_color_, light_position_, double_side_, cm_.color_map_, cm_.expansion_,
								 cm_.min_value_, cm_.max_value_);
}

} // namespace rendering

} // namespace cgogn
