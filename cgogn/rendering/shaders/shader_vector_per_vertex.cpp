/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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


#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>


namespace cgogn
{

namespace rendering
{

ShaderVectorPerVertex* ShaderVectorPerVertex::instance_ = nullptr;

ShaderVectorPerVertex::ShaderVectorPerVertex()
{
	const char* vertex_shader_source =
		"#version 150\n"
		"in vec3 vertex_pos;\n"
		"in vec3 vertex_normal;\n"
		"out vec3 normal;\n"
		"void main()\n"
		"{\n"
		"	normal = vertex_normal;\n"
		"	gl_Position = vec4(vertex_pos,1.0);\n"
		"}\n";

	const char* fragment_shader_source =
		"#version 150\n"
		 "uniform vec4 color;\n"
		 "out vec4 fragColor;\n"
		 "void main()\n"
		 "{\n"
		 "	fragColor = color;\n"
		 "}\n";

	const char* geometry_shader_source =
		"#version 150\n"
		"layout(points) in;\n"
		"layout(line_strip,max_vertices=2) out;\n"
		"in vec3 normal[];\n"
		"uniform mat4 projection_matrix;\n"
		"uniform mat4 model_view_matrix;\n"
		"uniform float length;\n"
		"void main()\n"
		"{\n"
		"	gl_Position = projection_matrix * model_view_matrix * gl_in[0].gl_Position;\n"
		"	EmitVertex();\n"
		"	vec4 end_point = gl_in[0].gl_Position + vec4(length * normal[0], 0.0);\n"
		"	gl_Position = projection_matrix * model_view_matrix * end_point;\n"
		"	EmitVertex();\n"
		"	EndPrimitive();\n"
		"}\n";

	load3_bind(vertex_shader_source,fragment_shader_source,geometry_shader_source,
		 "vertex_pos","vertex_normal");
	add_uniforms("color","length");
}


} // namespace rendering

} // namespace cgogn
