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

#include <iostream>

#include <cgogn/rendering/transparency_shaders/shader_transparent_flat.h>


namespace cgogn
{

namespace rendering
{


const char* ShaderFlatTransp::vertex_shader_source_ =
"#version 330\n"
"in vec3 vertex_pos;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"out vec3 pos;\n"
"out vec4 projCoord;\n"
"void main()\n"
"{\n"
"	vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);\n"
"	pos = pos4.xyz;"
"	gl_Position = projection_matrix * pos4;\n"
"   projCoord = gl_Position;\n"
"}\n";

const char* ShaderFlatTransp::fragment_shader_source_ =
"#version 330\n"
"layout(location = 0) out vec4 color_out;\n"
"layout(location = 1) out float depth_out;\n"
"in vec4 projCoord;\n"
"uniform vec4 front_color;\n"
"uniform vec4 back_color;\n"
"uniform vec4 ambiant_color;\n"
"uniform vec3 lightPosition;\n"
"uniform sampler2D rgba_texture ;\n"
"uniform sampler2D depth_texture ;\n"
"uniform int layer;\n"
"uniform bool cull_back_face;\n"
"uniform bool lighted;\n"
"in vec3 pos;\n"
"void main()\n"
"{\n"
"	vec3 tc = 0.5*projCoord.xyz/projCoord.w + vec3(0.5,0.5,0.5);\n"
"	if ((layer>0) && (tc.z <= texture(depth_texture, tc.xy).r)) discard;\n"
"	float lambert = 1.0;\n"
"	if (lighted)\n"
"	{\n"
"		vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));\n"
"		vec3 L = normalize(lightPosition-pos);\n"
"		lambert = dot(N,L);\n"
"	}\n"
"	vec3 color;\n"
"	float alpha;\n"
"	if (gl_FrontFacing)\n"
"		{ color = ambiant_color.rgb+lambert*front_color.rgb; alpha = front_color.a;}\n"
"	else \n"
"		if (cull_back_face) discard;\n"
"		 else { color = ambiant_color.rgb+lambert*back_color.rgb; alpha = back_color.a;}\n"
"	vec4 dst = texture(rgba_texture, tc.xy);\n"
"	color_out = vec4( (dst.a * (alpha * color) + dst.rgb), (1.0-alpha)*dst.a);"
"	depth_out = tc.z;\n"
"}\n";

ShaderFlatTransp* ShaderFlatTransp::instance_ = nullptr;

ShaderFlatTransp::ShaderFlatTransp()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
	prg_.link();
	get_matrices_uniforms();

	unif_front_color_ = prg_.uniformLocation("front_color");
	unif_back_color_ = prg_.uniformLocation("back_color");
	unif_ambiant_color_ = prg_.uniformLocation("ambiant_color");
	unif_light_position_ = prg_.uniformLocation("lightPosition");
	unif_lighted_ = prg_.uniformLocation("lighted");
	unif_layer_ = prg_.uniformLocation("layer");
	unif_bf_culling_ = prg_.uniformLocation("cull_back_face");
	unif_rgba_texture_sampler_ =  prg_.uniformLocation("rgba_texture");
	unif_depth_texture_sampler_ = prg_.uniformLocation("depth_texture");
}

void ShaderFlatTransp::set_light_position(const QVector3D& l)
{
	prg_.setUniformValue(unif_light_position_, l);
}

void ShaderFlatTransp::set_front_color(const QColor& rgb)
{
	if (unif_front_color_ >= 0)
		prg_.setUniformValue(unif_front_color_, rgb);
}

void ShaderFlatTransp::set_back_color(const QColor& rgb)
{
	if (unif_back_color_ >= 0)
		prg_.setUniformValue(unif_back_color_, rgb);
}

void ShaderFlatTransp::set_ambiant_color(const QColor& rgb)
{
	prg_.setUniformValue(unif_ambiant_color_, rgb);
}

void ShaderFlatTransp::set_bf_culling(bool cull)
{
	prg_.setUniformValue(unif_bf_culling_, cull);
}

void ShaderFlatTransp::set_lighted(bool lighted)
{
	prg_.setUniformValue(unif_lighted_, lighted);
}


void ShaderFlatTransp::set_layer(int layer)
{
	prg_.setUniformValue(unif_layer_, layer);
}

void ShaderFlatTransp::set_rgba_sampler(GLuint rgba_samp)
{
	prg_.setUniformValue(unif_rgba_texture_sampler_, rgba_samp);
}

void ShaderFlatTransp::set_depth_sampler(GLuint depth_samp)
{
	prg_.setUniformValue(unif_depth_texture_sampler_, depth_samp);
}


std::unique_ptr< ShaderFlatTransp::Param> ShaderFlatTransp::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderFlatTransp();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

ShaderFlatTransp* ShaderFlatTransp::get_instance()
{
	if (!instance_)
	{
		instance_ = new ShaderFlatTransp();
		ShaderProgram::register_instance(instance_);
	}
	return instance_;
}




ShaderParamFlatTransp::ShaderParamFlatTransp(ShaderFlatTransp* sh) :
	ShaderParam(sh),
	front_color_(250, 0, 0),
	back_color_(0, 250, 0),
	ambiant_color_(5, 5, 5),
	light_pos_(10, 100, 1000),
	bf_culling_(false),
	lighted_(true)
{}

void ShaderParamFlatTransp::set_position_vbo(VBO* vbo_pos)
{
	QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
	shader_->bind();
	vao_->bind();
	// position vbo
	vbo_pos->bind();
	ogl->glEnableVertexAttribArray(ShaderFlatTransp::ATTRIB_POS);
	ogl->glVertexAttribPointer(ShaderFlatTransp::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	vbo_pos->release();
	vao_->release();
	shader_->release();
}

void ShaderParamFlatTransp::set_uniforms()
{
	ShaderFlatTransp* sh = static_cast<ShaderFlatTransp*>(this->shader_);
	sh->set_front_color(front_color_);
	sh->set_back_color(back_color_);
	sh->set_ambiant_color(ambiant_color_);
	sh->set_light_position(light_pos_);
	sh->set_bf_culling(bf_culling_);
	sh->set_lighted(lighted_);
}

void ShaderParamFlatTransp::set_alpha(int alpha)
{
	front_color_.setAlpha(alpha);
	back_color_.setAlpha(alpha);
}


} // namespace rendering

} // namespace cgogn
