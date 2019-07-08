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


#include <cgogn/rendering/transparency_shaders/shader_transparent_volumes.h>

#include <QOpenGLFunctions>
#include <iostream>
#include<QColor>
#include<QImage>

namespace cgogn
{

namespace rendering
{


const char* ShaderTransparentVolumes::vertex_shader_source_ =
"#version 330\n"
"in vec3 vertex_pos;\n"
"void main()\n"
"{\n"
"   gl_Position = vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderTransparentVolumes::geometry_shader_source_ =
"#version 330\n"
"layout (lines_adjacency) in;\n"
"layout (triangle_strip, max_vertices=3) out;\n"
"out vec3 color_f;\n"
"out float alpha;\n"
"out vec4 projCoord;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform mat3 normal_matrix;\n"
"uniform float explode_vol;\n"
"uniform vec3 light_position;\n"
"uniform vec4 color;\n"
"uniform vec4 plane_clip;\n"
"uniform vec4 plane_clip2;\n"
"uniform bool lighted;\n"
"void main()\n"
"{\n"
"	float d = dot(plane_clip,gl_in[0].gl_Position);\n"
"	float d2 = dot(plane_clip2,gl_in[0].gl_Position);\n"
"	if ((d<=0.0)&&(d2<=0.0))\n"
"	{\n"
"		vec4 face_center =  model_view_matrix * gl_in[1].gl_Position;\n"
"		float lambertTerm = 1.0;\n"
"		if (lighted)\n"
"		{\n"
"			vec3 v1 = gl_in[2].gl_Position.xyz - gl_in[1].gl_Position.xyz;\n"
"			vec3 v2 = gl_in[3].gl_Position.xyz - gl_in[1].gl_Position.xyz;\n"
"			vec3 N  = normalize(normal_matrix*cross(v1,v2));\n"
"			vec3 L =  normalize (light_position - face_center.xyz);\n"
"			lambertTerm = abs(dot(N,L));\n"
"		}\n"
"		for (int i=1; i<=3; i++)\n"
"		{\n"
"			vec4 Q = explode_vol *  gl_in[i].gl_Position  + (1.0-explode_vol) * gl_in[0].gl_Position;\n"
"			gl_Position = projection_matrix * model_view_matrix *  Q;\n"
"			color_f = color.rgb * lambertTerm;\n"
"			alpha = color.a;\n"
"			projCoord = gl_Position;\n"
"			EmitVertex();\n"
"		}\n"
"		EndPrimitive();\n"
"	}\n"
"}\n";

const char* ShaderTransparentVolumes::fragment_shader_source_ =
"#version 330\n"
"layout(location = 0) out vec4 color_out;\n"
"layout(location = 1) out float depth_out;\n"
"in vec4 projCoord;\n"
"in vec3 color_f;\n"
"in float alpha;\n"
"uniform sampler2D rgba_texture ;\n"
"uniform sampler2D depth_texture ;\n"
"uniform int layer;\n"
"uniform bool cull_back_face;\n"


"void main()\n"
"{\n"
"	if (!gl_FrontFacing && cull_back_face) discard;\n"
"	vec3 tc = 0.5*projCoord.xyz/projCoord.w + vec3(0.5,0.5,0.5);\n"
"	if ((layer>0) && (tc.z <= texture(depth_texture, tc.xy).r)) discard;\n"
"	vec4 dst = texture(rgba_texture, tc.xy);\n"
"	color_out = vec4( (dst.a * (alpha * color_f) + dst.rgb), (1.0-alpha)*dst.a);"
"	depth_out = tc.z;\n"
"}\n";



ShaderTransparentVolumes* ShaderTransparentVolumes::instance_ = nullptr;


ShaderTransparentVolumes::ShaderTransparentVolumes()
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
	unif_light_position_ = prg_.uniformLocation("light_position");
	unif_color_ = prg_.uniformLocation("color");

	unif_bf_culling_ = prg_.uniformLocation("cull_back_face");
	unif_lighted_ = prg_.uniformLocation("lighted");
	unif_layer_ = prg_.uniformLocation("layer");
	unif_depth_texture_sampler_ = prg_.uniformLocation("depth_texture");
	unif_rgba_texture_sampler_ = prg_.uniformLocation("rgba_texture");

	// default param
	bind();
	set_light_position(QVector3D(10.0f,100.0f,1000.0f));
	set_explode_volume(0.8f);
	set_color(QColor(255,0,0));
	set_plane_clip(QVector4D(0,0,0,0));
	set_plane_clip2(QVector4D(0,0,0,0));
	release();
}


void ShaderTransparentVolumes::set_color(const QColor& rgb)
{
		prg_.setUniformValue(unif_color_, rgb);
}

void ShaderTransparentVolumes::set_light_position(const QVector3D& l)
{
	prg_.setUniformValue(unif_light_position_, l);
}

void ShaderTransparentVolumes::set_explode_volume(float32 x)
{
	prg_.setUniformValue(unif_expl_v_, x);
}

void ShaderTransparentVolumes::set_plane_clip(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip_, plane);
}

void ShaderTransparentVolumes::set_plane_clip2(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip2_, plane);
}

void ShaderTransparentVolumes::set_bf_culling(bool cull)
{
	prg_.setUniformValue(unif_bf_culling_, cull);
}

void ShaderTransparentVolumes::set_lighted(bool lighted)
{
	prg_.setUniformValue(unif_lighted_, lighted);
}

void ShaderTransparentVolumes::set_layer(int layer)
{
	prg_.setUniformValue(unif_layer_, layer);
}

void ShaderTransparentVolumes::set_rgba_sampler(GLuint rgba_samp)
{
	prg_.setUniformValue(unif_rgba_texture_sampler_, rgba_samp);
}

void ShaderTransparentVolumes::set_depth_sampler(GLuint depth_samp)
{
	prg_.setUniformValue(unif_depth_texture_sampler_, depth_samp);
}


std::unique_ptr< ShaderTransparentVolumes::Param> ShaderTransparentVolumes::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderTransparentVolumes();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<ShaderTransparentVolumes::Param>(instance_);
}

ShaderTransparentVolumes* ShaderTransparentVolumes::get_instance()
{
	if (!instance_)
	{
		instance_ = new ShaderTransparentVolumes();
		ShaderProgram::register_instance(instance_);
	}
	return instance_;
}



ShaderParamTransparentVolumes::ShaderParamTransparentVolumes(ShaderTransparentVolumes* sh) :
	ShaderParam(sh),
	color_(255, 0, 0),
	plane_clip_(0, 0, 0, 0),
	plane_clip2_(0, 0, 0, 0),
	light_position_(10.0f, 100.0f, 1000.0f),
	explode_factor_(0.8f),
	bf_culling_(false),
	lighted_(true)
{}

void ShaderParamTransparentVolumes::set_position_vbo(VBO* vbo_pos)
{
	QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
	shader_->bind();
	vao_->bind();
	vbo_pos->bind();
	ogl->glEnableVertexAttribArray(ShaderTransparentVolumes::ATTRIB_POS);
	ogl->glVertexAttribPointer(ShaderTransparentVolumes::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
	vbo_pos->release();
	vao_->release();
	shader_->release();

}

void ShaderParamTransparentVolumes::set_uniforms()
{
	ShaderTransparentVolumes* sh = static_cast<ShaderTransparentVolumes*>(this->shader_);
	sh->set_color(color_);
	sh->set_explode_volume(explode_factor_);
	sh->set_light_position(light_position_);
	sh->set_plane_clip(plane_clip_);
	sh->set_plane_clip2(plane_clip2_);
	sh->set_bf_culling(bf_culling_);
	sh->set_lighted(lighted_);
}


} // namespace rendering

} // namespace cgogn
