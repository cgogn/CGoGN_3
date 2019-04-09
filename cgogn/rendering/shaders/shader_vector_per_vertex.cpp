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
#include <cgogn/core/utils/unique_ptr.h>

#include <QColor>
#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderVectorPerVertex* ShaderVectorPerVertex::instance_ = nullptr;

std::unique_ptr<ShaderVectorPerVertex::Param> ShaderVectorPerVertex::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderVectorPerVertex();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

const char* ShaderVectorPerVertex::vertex_shader_source_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"in vec3 vertex_normal;\n"
"out vec3 normal;\n"
"void main()\n"
"{\n"
"	normal = vertex_normal;\n"
"	gl_Position = vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderVectorPerVertex::geometry_shader_source_ =
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

const char* ShaderVectorPerVertex::fragment_shader_source_ =
"#version 150\n"
"uniform vec4 color;\n"
"out vec4 fragColor;\n"
"void main()\n"
"{\n"
"	fragColor = color;\n"
"}\n";

ShaderVectorPerVertex::ShaderVectorPerVertex()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Geometry, geometry_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
	prg_.bindAttributeLocation("vertex_normal", ATTRIB_NORMAL);
	prg_.link();

	get_matrices_uniforms();

	unif_color_ = prg_.uniformLocation("color");
	unif_length_ = prg_.uniformLocation("length");

	//default param
	set_color(QColor(255, 255, 255));
	set_length(1.0);
}

void ShaderVectorPerVertex::set_color(const QColor& rgb)
{
	prg_.setUniformValue(unif_color_, rgb);
}

void ShaderVectorPerVertex::set_length(float32 l)
{
	prg_.setUniformValue(unif_length_, l);
}

} // namespace rendering

} // namespace cgogn
