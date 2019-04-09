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


#include <cgogn/rendering/shaders/shader_texture.h>
#include <QOpenGLFunctions>
#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderTexture* ShaderTexture::instance_ = nullptr;

const char* ShaderTexture::vertex_shader_source_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"in vec2 vertex_tc;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"out vec2 tc;\n"
"void main()\n"
"{\n"
"	tc = vertex_tc;\n"
"   gl_Position = projection_matrix * model_view_matrix * vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderTexture::fragment_shader_source_ =
"#version 150\n"
"out vec4 frag_color;\n"
"uniform sampler2D texture_unit;\n"
"in vec2 tc;\n"
"void main()\n"
"{\n"
"	frag_color = texture(texture_unit,tc);\n"
"}\n";

ShaderTexture::ShaderTexture()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
	prg_.bindAttributeLocation("vertex_tc", ATTRIB_TC);
	prg_.link();
	get_matrices_uniforms();
	prg_.setUniformValue("texture_unit", 0);
}

ShaderParamTexture::ShaderParamTexture(ShaderTexture* sh) :
	ShaderParam(sh),
	texture_(nullptr)
{}

void ShaderParamTexture::set_uniforms()
{
	if (texture_)
	{
		QOpenGLContext::currentContext()->functions()->glActiveTexture(GL_TEXTURE0);
		texture_->bind();
	}
}

void ShaderParamTexture::set_vbo(VBO* vbo_pos, VBO* vbo_tc)
{
	QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();

	shader_->bind();
	vao_->bind();

	// position vbo
	vbo_pos->bind();
	ogl->glEnableVertexAttribArray(ShaderTexture::ATTRIB_POS);
	ogl->glVertexAttribPointer(ShaderTexture::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	vbo_pos->release();

	// color  vbo
	vbo_tc->bind();
	ogl->glEnableVertexAttribArray(ShaderTexture::ATTRIB_TC);
	ogl->glVertexAttribPointer(ShaderTexture::ATTRIB_TC, vbo_tc->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	vbo_tc->release();

	vao_->release();
	shader_->release();
}

} // namespace rendering

} // namespace cgogn
