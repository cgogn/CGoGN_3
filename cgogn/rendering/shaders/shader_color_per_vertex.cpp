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

#include <cgogn/rendering/shaders/shader_color_per_vertex.h>

#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderFlatColorPerVertex* ShaderFlatColorPerVertex::instance_ = nullptr;
ShaderPhongColorPerVertex* ShaderPhongColorPerVertex::instance_ = nullptr;

static const char* vertex_shader_source = "in vec3 vertex_pos;\n"
										  "in vec3 vertex_color;\n"
										  "uniform mat4 projection_matrix;\n"
										  "uniform mat4 model_view_matrix;\n"
										  "out vec3 pos_v;\n"
										  "out vec3 color_v;\n"
										  "#if WITH_NORMAL==1\n"
										  "uniform mat3 normal_matrix;\n"
										  "in vec3 vertex_normal;\n"
										  "out vec3 normal_v;\n"
										  "#endif\n"
										  "void main()\n"
										  "{\n"
										  "	color_v = vertex_color;\n"
										  "	vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);\n"
										  "	pos_v = pos4.xyz;\n"
										  "#if WITH_NORMAL==1\n"
										  "	normal_v = normal_matrix * vertex_normal;\n"
										  "#endif\n"
										  "   gl_Position = projection_matrix * pos4;\n"
										  "}\n";

static const char* fragment_shader_source = "in vec3 pos_v;\n"
											"in vec3 color_v;\n"
											"#if WITH_NORMAL==1\n"
											"in vec3 normal_v;\n"
											"#endif\n"
											"uniform vec3 light_position;\n"
											"out vec3 fragColor;\n"
											"void main()\n"
											"{\n"
											"#if WITH_NORMAL==1\n"
											"	vec3 N = normal_v;\n"
											"#else\n"
											"	vec3 N = normalize(cross(dFdx(pos_v),dFdy(pos_v)));\n"
											"#endif\n"
											"	vec3 L = normalize(light_position-pos_v);\n"
											"	float lambert = dot(N,L);\n"
											"	fragColor = lambert*color_v;\n"
											"}\n";

ShaderFlatColorPerVertex::ShaderFlatColorPerVertex()
{
	std::string bs("#version 150\n#define WITH_NORMAL 0\n");

	std::string vs = bs + std::string(vertex_shader_source);
	std::string fs = bs + std::string(fragment_shader_source);

	load2_bind(vs, fs, "vertex_pos", "vertex_color");
	add_uniforms("light_position");
}

ShaderPhongColorPerVertex::ShaderPhongColorPerVertex()
{
	std::string bs("#version 150\n#define WITH_NORMAL 1\n");

	std::string vs = bs + std::string(vertex_shader_source);
	std::string fs = bs + std::string(fragment_shader_source);

	load2_bind(vs, fs, "vertex_pos", "vertex_normal", "vertex_color");
	add_uniforms("light_position");
}

} // namespace rendering

} // namespace cgogn
