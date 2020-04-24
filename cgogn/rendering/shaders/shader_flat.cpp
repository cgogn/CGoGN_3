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

#include <cgogn/rendering/shaders/shader_flat.h>

namespace cgogn
{

namespace rendering
{
const char* vertex_shader_source =
	R"(#version 150
in vec3 vertex_pos;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
out vec3 pos;
void main()
{
vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);
pos = pos4.xyz;
  gl_Position = projection_matrix * pos4;
}
)";

const char* fragment_shader_source =
	R"(#version 150
out vec4 fragColor;
uniform vec4 front_color;
uniform vec4 back_color;
uniform vec4 ambiant_color;
uniform vec3 light_position;
uniform bool double_side;
in vec3 pos;
void main()
{
	vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));
	vec3 L = normalize(light_position-pos);
	float lambert = dot(N,L);
	if (gl_FrontFacing)
		fragColor = vec4(ambiant_color.rgb+lambert*front_color.rgb, front_color.a);
	else
		if (!double_side) discard;
		else fragColor = vec4(ambiant_color.rgb+lambert*back_color.rgb, back_color.a);
}
)";

ShaderFlat* ShaderFlat::instance_ = nullptr;

ShaderFlat::ShaderFlat()
{

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos");
	add_uniforms("front_color", "back_color", "ambiant_color", "light_position", "double_side");
}

void ShaderParamFlat::set_uniforms()
{
	shader_->set_uniforms_values(front_color_, back_color_, ambiant_color_, light_position_, double_side_);
}

} // namespace rendering

} // namespace cgogn
