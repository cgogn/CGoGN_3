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

#define CGOGN_RENDER_FLAT_TR_DR_CPP_

#include <iostream>

#include <cgogn/rendering/transparency_shaders/shader_copy_depth.h>

namespace cgogn
{

namespace rendering
{

const char* ShaderCopyDepth::vertex_shader_source_ =
"#version 330\n"
"out vec2 tc;\n"
"void main(void)\n"
"{\n"
"	vec4 vertices[4] = vec4[4](vec4(-1.0, -1.0, 0.0, 1.0), vec4(1.0, -1.0, 0.0, 1.0), vec4(-1.0, 1.0, 0.0, 1.0), vec4(1.0, 1.0, 0.0, 1.0));\n"
"	gl_Position = vertices[gl_VertexID];\n"
"	tc = (gl_Position.xy + vec2(1,1))/2;\n"
"}\n";

const char* ShaderCopyDepth::fragment_shader_source_ =
"#version 330\n"
"uniform sampler2D depth_texture;\n"
"in vec2 tc;\n"
"void main(void)\n"
"{\n"
"	gl_FragDepth = texture(depth_texture,tc).r;\n"
"}\n";

ShaderCopyDepth* ShaderCopyDepth::instance_ = nullptr;

ShaderCopyDepth::ShaderCopyDepth()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.link();

	unif_depth_texture_sampler_ = prg_.uniformLocation("depth_texture");
}


void ShaderCopyDepth::set_depth_sampler(GLuint depth_samp)
{
	prg_.setUniformValue(unif_depth_texture_sampler_, depth_samp);
}

std::unique_ptr< ShaderCopyDepth::Param> ShaderCopyDepth::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderCopyDepth();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

ShaderParamCopyDepth::ShaderParamCopyDepth(ShaderCopyDepth* sh) :
	ShaderParam(sh),
	texture_(nullptr)
{}

void ShaderParamCopyDepth::set_uniforms()
{
	ShaderCopyDepth* sh = static_cast<ShaderCopyDepth*>(this->shader_);
	sh->set_depth_sampler(0);
	QOpenGLContext::currentContext()->functions()->glActiveTexture(GL_TEXTURE0);
	if (texture_ != nullptr)
		texture_->bind();
}

} // namespace rendering

} // namespace cgogn
