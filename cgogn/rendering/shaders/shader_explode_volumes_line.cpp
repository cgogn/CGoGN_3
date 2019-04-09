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

#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>
#include <cgogn/core/utils/unique_ptr.h>

#include <iostream>
#include <QColor>
#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

ShaderExplodeVolumesLine* ShaderExplodeVolumesLine::instance_ = nullptr;

std::unique_ptr<ShaderExplodeVolumesLine::Param> ShaderExplodeVolumesLine::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderExplodeVolumesLine();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

const char* ShaderExplodeVolumesLine::vertex_shader_source_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"void main()\n"
"{\n"
"   gl_Position = vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderExplodeVolumesLine::geometry_shader_source_ =
"#version 150\n"
"layout (triangles) in;\n"
"layout (line_strip, max_vertices=2) out;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform float explode_vol;\n"
"uniform vec4 plane_clip;\n"
"uniform vec4 plane_clip2;\n"
"void main()\n"
"{\n"
"	float d = dot(plane_clip,gl_in[0].gl_Position);\n"
"	float d2 = dot(plane_clip2,gl_in[0].gl_Position);\n"
"	if ((d<=0.0)&&(d2<=0.0))\n"
"	{\n"
"		for (int i=1; i<=2; i++)\n"
"		{\n"
"			vec4 Q = explode_vol *  gl_in[i].gl_Position  + (1.0-explode_vol) * gl_in[0].gl_Position;\n"
"			gl_Position = projection_matrix * model_view_matrix *  Q;\n"
"			EmitVertex();\n"
"		}\n"
"		EndPrimitive();\n"
"	}\n"
"}\n";

const char* ShaderExplodeVolumesLine::fragment_shader_source_ =
"#version 150\n"
"uniform vec4 color;\n"
"out vec4 fragColor;\n"
"void main()\n"
"{\n"
"   fragColor = color;\n"
"}\n";

ShaderExplodeVolumesLine::ShaderExplodeVolumesLine()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Geometry, geometry_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
	prg_.link();
	get_matrices_uniforms();
	unif_expl_v_ = prg_.uniformLocation("explode_vol");
	unif_plane_clip_ = prg_.uniformLocation("plane_clip");
	unif_plane_clip2_ = prg_.uniformLocation("plane_clip2");
	unif_color_ = prg_.uniformLocation("color");

	// default param
	bind();
	set_explode_volume(0.8f);
	set_color(QColor(255, 255, 255));
	set_plane_clip(QVector4D(0, 0, 0, 0));
	set_plane_clip2(QVector4D(0, 0, 0, 0));
	release();
}

void ShaderExplodeVolumesLine::set_color(const QColor& rgb)
{
	if (unif_color_ >= 0)
		prg_.setUniformValue(unif_color_, rgb);
}

void ShaderExplodeVolumesLine::set_explode_volume(float32 x)
{
	prg_.setUniformValue(unif_expl_v_, x);
}

void ShaderExplodeVolumesLine::set_plane_clip(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip_, plane);
}

void ShaderExplodeVolumesLine::set_plane_clip2(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip2_, plane);
}

} // namespace rendering

} // namespace cgogn
