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

#include <cgogn/rendering/transparency_shaders/shader_transparent_phong.h>


namespace cgogn
{

namespace rendering
{


const char* ShaderPhongTransp::vertex_shader_source_ =
"#version 330\n"
"in vec3 vertex_pos;\n"
"in vec3 vertex_normal;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform mat3 normal_matrix;\n"
"out vec3 Normal;\n"
"out vec3 pos;\n"
"out vec4 projCoord;\n"
"void main()\n"
"{\n"
"	Normal = normal_matrix * vertex_normal;\n"
"	vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);\n"
"	pos = pos4.xyz;"
"	gl_Position = projection_matrix * pos4;\n"
"   projCoord = gl_Position;\n"
"}\n";

const char* ShaderPhongTransp::fragment_shader_source_ =
"#version 330\n"
"layout(location = 0) out vec4 color_out;\n"
"layout(location = 1) out float depth_out;\n"
"in vec3 pos;\n"
"in vec4 projCoord;\n"
"in vec3 Normal;\n"
"uniform vec4 front_color;\n"
"uniform vec4 back_color;\n"
"uniform vec4 ambiant_color;\n"
"uniform vec4 spec_color;\n"
"uniform float spec_coef;\n"
"uniform vec3 lightPosition;\n"
"uniform sampler2D rgba_texture ;\n"
"uniform sampler2D depth_texture ;\n"
"uniform int layer;\n"
"uniform bool cull_back_face;\n"
"void main()\n"
"{\n"
"	vec3 EyeVector = -pos;"
"	vec3 N = normalize (Normal);\n"
"	vec3 L = normalize (lightPosition - pos);\n"
"	vec3 E = normalize(-pos);\n"
"	vec3 R = reflect(-L, N);\n"
"	float specular = pow( max(dot(R, E), 0.0), spec_coef );\n"

"	vec3 tc = 0.5*projCoord.xyz/projCoord.w + vec3(0.5,0.5,0.5);\n"
"	if ((layer>0) && (tc.z <= texture(depth_texture, tc.xy).r)) discard;\n"
"	if (gl_FrontFacing==false)\n" // do not use ! because of bug on old intel under OS/X
"	{\n"
"		N *= -1.0;\n"
"	}\n"
"	float lambert = clamp(dot(N,L),0.0,1.0);\n"
"	vec3 color;\n"
"	float alpha;\n"
"	if (gl_FrontFacing)\n"
"		{ color = ambiant_color.rgb+lambert*front_color.rgb+spec_color.rgb*specular; alpha = front_color.a;}\n"
"	else \n"
"		if (cull_back_face) discard;\n"
"		 else { color = ambiant_color.rgb+lambert*back_color.rgb+spec_color.rgb*specular; alpha = back_color.a;}\n"
"	vec4 dst = texture(rgba_texture, tc.xy);\n"
"	color_out = vec4( (dst.a * (alpha * color) + dst.rgb), (1.0-alpha)*dst.a);"
"	depth_out = tc.z;\n"
"}\n";

ShaderPhongTransp* ShaderPhongTransp::instance_ = nullptr;

ShaderPhongTransp::ShaderPhongTransp()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
	prg_.link();
	get_matrices_uniforms();

	unif_front_color_ = prg_.uniformLocation("front_color");
	unif_back_color_ = prg_.uniformLocation("back_color");
	unif_ambiant_color_ = prg_.uniformLocation("ambiant_color");
	unif_specular_color_ = prg_.uniformLocation("spec_color");
	unif_specular_coef_ = prg_.uniformLocation("spec_coef");

	unif_light_position_ = prg_.uniformLocation("lightPosition");
	unif_layer_ = prg_.uniformLocation("layer");
	unif_bf_culling_ = prg_.uniformLocation("cull_back_face");
	unif_rgba_texture_sampler_ =  prg_.uniformLocation("rgba_texture");
	unif_depth_texture_sampler_ = prg_.uniformLocation("depth_texture");
}

void ShaderPhongTransp::set_light_position(const QVector3D& l)
{
	prg_.setUniformValue(unif_light_position_, l);
}

void ShaderPhongTransp::set_front_color(const QColor& rgb)
{
	if (unif_front_color_ >= 0)
		prg_.setUniformValue(unif_front_color_, rgb);
}

void ShaderPhongTransp::set_back_color(const QColor& rgb)
{
	if (unif_back_color_ >= 0)
		prg_.setUniformValue(unif_back_color_, rgb);
}

void ShaderPhongTransp::set_ambiant_color(const QColor& rgb)
{
	prg_.setUniformValue(unif_ambiant_color_, rgb);
}

void ShaderPhongTransp::set_specular_color(const QColor& rgb)
{
	prg_.setUniformValue(unif_specular_color_, rgb);
}

void ShaderPhongTransp::set_specular_coef(GLfloat coef)
{
	prg_.setUniformValue(unif_specular_coef_, coef);
}


void ShaderPhongTransp::set_bf_culling(bool cull)
{
	prg_.setUniformValue(unif_bf_culling_, cull);
}


void ShaderPhongTransp::set_layer(int layer)
{
	prg_.setUniformValue(unif_layer_, layer);
}

void ShaderPhongTransp::set_rgba_sampler(GLuint rgba_samp)
{
	prg_.setUniformValue(unif_rgba_texture_sampler_, rgba_samp);
}

void ShaderPhongTransp::set_depth_sampler(GLuint depth_samp)
{
	prg_.setUniformValue(unif_depth_texture_sampler_, depth_samp);
}


std::unique_ptr< ShaderPhongTransp::Param> ShaderPhongTransp::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderPhongTransp();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

ShaderPhongTransp* ShaderPhongTransp::get_instance()
{
	if (!instance_)
	{
		instance_ = new ShaderPhongTransp();
		ShaderProgram::register_instance(instance_);
	}
	return instance_;
}




ShaderParamPhongTransp::ShaderParamPhongTransp(ShaderPhongTransp* sh) :
	ShaderParam(sh),
	front_color_(250, 0, 0),
	back_color_(0, 250, 0),
	ambiant_color_(5, 5, 5),
	specular_color_(255, 255, 255),
	specular_coef_(100),
	light_pos_(10, 100, 1000),
	bf_culling_(false)
{}

void ShaderParamPhongTransp::set_position_vbo(VBO* vbo_pos)
{
	QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
	shader_->bind();
	vao_->bind();
	// position vbo
	vbo_pos->bind();
	ogl->glEnableVertexAttribArray(ShaderPhongTransp::ATTRIB_POS);
	ogl->glVertexAttribPointer(ShaderPhongTransp::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	vbo_pos->release();
	vao_->release();
	shader_->release();
}

void ShaderParamPhongTransp::set_normal_vbo(VBO* vbo_normal)
{
	QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
	shader_->bind();
	vao_->bind();
	vbo_normal->bind();
	ogl->glEnableVertexAttribArray(ShaderPhongTransp::ATTRIB_NORM);
	ogl->glVertexAttribPointer(ShaderPhongTransp::ATTRIB_NORM, vbo_normal->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	vbo_normal->release();
	vao_->release();
	shader_->release();
}



void ShaderParamPhongTransp::set_uniforms()
{
	ShaderPhongTransp* sh = static_cast<ShaderPhongTransp*>(this->shader_);
	sh->set_front_color(front_color_);
	sh->set_back_color(back_color_);
	sh->set_ambiant_color(ambiant_color_);
	sh->set_specular_color(specular_color_);
	sh->set_specular_coef(specular_coef_);
	sh->set_light_position(light_pos_);
	sh->set_bf_culling(bf_culling_);
}

void ShaderParamPhongTransp::set_alpha(int alpha)
{
	front_color_.setAlpha(alpha);
	back_color_.setAlpha(alpha);
}

} // namespace rendering

} // namespace cgogn
