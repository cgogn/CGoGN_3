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

#include <cgogn/rendering/transparency_shaders/shader_transparent_quad.h>


namespace cgogn
{

namespace rendering
{

const char* ShaderTranspQuad::vertex_shader_source_ =
"#version 330\n"
"out vec2 tc;\n"
"void main(void)\n"
"{\n"
"	vec4 vertices[4] = vec4[4](vec4(-1.0, -1.0, 0.0, 1.0), vec4(1.0, -1.0, 0.0, 1.0), vec4(-1.0, 1.0, 0.0, 1.0), vec4(1.0, 1.0, 0.0, 1.0));\n"
"	gl_Position = vertices[gl_VertexID];\n"
"	tc = (gl_Position.xy + vec2(1,1))/2;\n"
"}\n";

const char* ShaderTranspQuad::fragment_shader_source_ =
"#version 330\n"
"uniform sampler2D rgba_texture;\n"
"uniform sampler2D depth_texture;\n"
"in vec2 tc;\n"
"out vec4 fragColor;\n"
"void main(void)\n"
"{\n"
"	vec4 color = texture(rgba_texture, tc);"
"	if (color.a <= 0.0) discard;"
"	gl_FragDepth = texture(depth_texture, tc).r;\n"
"	fragColor = color;\n"
"}\n";

ShaderTranspQuad* ShaderTranspQuad::instance_ = nullptr;

ShaderTranspQuad::ShaderTranspQuad()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.link();

	unif_rgba_texture_sampler_ =  prg_.uniformLocation("rgba_texture");
	unif_depth_texture_sampler_ = prg_.uniformLocation("depth_texture");
}

void ShaderTranspQuad::set_rgba_sampler(GLuint rgba_samp)
{
	prg_.setUniformValue(unif_rgba_texture_sampler_, rgba_samp);
}

void ShaderTranspQuad::set_depth_sampler(GLuint depth_samp)
{
	prg_.setUniformValue(unif_depth_texture_sampler_, depth_samp);
}

std::unique_ptr< ShaderTranspQuad::Param> ShaderTranspQuad::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderTranspQuad();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}



ShaderParamTranspQuad::ShaderParamTranspQuad(ShaderTranspQuad* sh) :
	ShaderParam(sh)
{}


void ShaderParamTranspQuad::set_uniforms()
{
	ShaderTranspQuad* sh = static_cast<ShaderTranspQuad*>(this->shader_);
	sh->set_rgba_sampler(rgba_texture_sampler_);
	sh->set_depth_sampler(depth_texture_sampler_);
}

} // namespace rendering

} // namespace cgogn
