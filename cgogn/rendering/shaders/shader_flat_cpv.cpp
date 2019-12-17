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

#include <cgogn/rendering/shaders/shader_flat_cpv.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatColor* ShaderFlatColor::instance_ = nullptr;

ShaderFlatColor::ShaderFlatColor()
{
	const char* vertex_shader_source = "#version 150\n"
									   "in vec3 vertex_pos;\n"
									   "in vec3 vertex_col;\n"
									   "uniform mat4 projection_matrix;\n"
									   "uniform mat4 model_view_matrix;\n"
									   "out vec3 pos;\n"
									   "out vec3 col;\n"
									   "void main()\n"
									   "{\n"
									   "	vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);\n"
									   "	pos = pos4.xyz;\n"
									   "	col = vertex_col;\n"
									   "   gl_Position = projection_matrix * pos4;\n"
									   "}\n";

	const char* fragment_shader_source = "#version 150\n"
										 "out vec4 fragColor;\n"
										 "uniform vec4 ambiant_color;\n"
										 "uniform vec3 lightPosition;\n"
										 "uniform bool cull_back_face;\n"
										 "in vec3 pos;\n"
										 "in vec3 col;\n"
										 "void main()\n"
										 "{\n"
										 "	vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));\n"
										 "	vec3 L = normalize(lightPosition-pos);\n"
										 "	float lambert = dot(N,L);\n"
										 "	if ((gl_FrontFacing==false) && cull_back_face) discard;\n"
										 "	else fragColor = ambiant_color+vec4(lambert*col,1.0);\n"
										 "}\n";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos", "vertex_col");

	add_uniforms("ambiant_color", "lightPosition", "cull_back_face");
}

} // namespace rendering

} // namespace cgogn
