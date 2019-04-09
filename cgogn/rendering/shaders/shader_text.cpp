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

#include <cgogn/rendering/shaders/shader_text.h>
#include <cgogn/core/utils/unique_ptr.h>

#include <QOpenGLFunctions>
#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderText* ShaderText::instance_ = nullptr;

const char* ShaderText::vertex_shader_source_ =
"#version 150\n"
"in vec4 vertex_in;\n"
"in float char_in;\n"
"in vec4 colsz_in;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform vec4 italic;\n"
"out vec2 tc;\n"
"out vec3 color;\n"
"const vec2 quad[4] = vec2[](vec2(-0.5,-1.0),vec2(0.5,-1.0),vec2(0.5,1.0),vec2(-0.5,1.));\n"
"const float tc_u[4] = float[](0,0.009,0.009,0);\n"
"const float tc_v[4] = float[](1,1,0,0);\n"
"void main()\n"
"{\n"
"	float size = colsz_in.w;\n"
"   vec4 P4 = model_view_matrix * vec4(vertex_in.xyz,1.0);\n"
"   P4[0] += size*vertex_in.w;\n"
"   P4 += vec4(size*(quad[gl_VertexID]+italic[gl_VertexID]),0,0);\n"
"	tc = vec2( char_in + tc_u[gl_VertexID], tc_v[gl_VertexID]);\n"
"   gl_Position = projection_matrix * P4;\n"
"	color = colsz_in.rgb;\n"
"}\n";

const char* ShaderText::fragment_shader_source_ =
"#version 150\n"
"out vec3 frag_color;\n"
"uniform sampler2D texture_unit;\n"
"in vec3 color;\n"
"in vec2 tc;\n"
"void main()\n"
"{\n"
"	float a = texture(texture_unit,tc).r;\n"
"	if (a==0)\n"
"		discard;\n"
"	else\n"
"		frag_color = a*color;\n"

"}\n";

ShaderText::ShaderText()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_in", ATTRIB_POS);
	prg_.bindAttributeLocation("char_in", ATTRIB_CHAR);
	prg_.bindAttributeLocation("colsz_in", ATTRIB_COLSZ);

	prg_.link();
	get_matrices_uniforms();
	unif_italic_ = prg_.uniformLocation("italic");
	bind();
	prg_.setUniformValue("texture_unit", 0);
	set_italic(0);
	release();
}

void ShaderText::set_italic(float32 i)
{
		prg_.setUniformValue(unif_italic_, QVector4D(-i,-i,i,i));
}

std::unique_ptr<ShaderText::Param> ShaderText::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderText();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<ShaderText::Param>(instance_);
}

ShaderParamText::ShaderParamText(ShaderText* sh) :
	ShaderParam(sh),
	texture_(nullptr)
{}

void ShaderParamText::set_uniforms()
{
	if (texture_)
	{
		QOpenGLContext::currentContext()->functions()->glActiveTexture(GL_TEXTURE0);
		texture_->bind();
	}

	ShaderText* sh = static_cast<ShaderText*>(this->shader_);
	sh->set_italic(italic_);

}

void ShaderParamText::set_vbo(VBO* vbo_pos, VBO* vbo_char, VBO* vbo_colsize)
{
	QOpenGLFunctions_3_3_Core * ogl = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();

	shader_->bind();
	vao_->bind();

	vbo_pos->bind();
	ogl->glEnableVertexAttribArray(ShaderText::ATTRIB_POS);
	ogl->glVertexAttribPointer(ShaderText::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	ogl->glVertexAttribDivisor(ShaderText::ATTRIB_POS,1u);
	vbo_pos->release();

	vbo_char->bind();
	ogl->glEnableVertexAttribArray(ShaderText::ATTRIB_CHAR);
	ogl->glVertexAttribPointer(ShaderText::ATTRIB_CHAR, vbo_char->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	ogl->glVertexAttribDivisor(ShaderText::ATTRIB_CHAR,1u);
	vbo_char->release();

	vbo_colsize->bind();
	ogl->glEnableVertexAttribArray(ShaderText::ATTRIB_COLSZ);
	ogl->glVertexAttribPointer(ShaderText::ATTRIB_COLSZ, vbo_colsize->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	ogl->glVertexAttribDivisor(ShaderText::ATTRIB_COLSZ, 1u);
	vbo_colsize->release();

	vao_->release();
	shader_->release();
}

} // namespace rendering

} // namespace cgogn
