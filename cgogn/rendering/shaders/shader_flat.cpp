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

#define CGOGN_RENDER_SHADERS_FLAT_CPP_


#include <iostream>

#include <cgogn/rendering/shaders/shader_flat.h>


namespace cgogn
{

namespace rendering
{

const char* ShaderFlatGen::vertex_shader_source_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"out vec3 pos;\n"
"void main()\n"
"{\n"
"	vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);\n"
"	pos = pos4.xyz;"
"   gl_Position = projection_matrix * pos4;\n"
"}\n";

const char* ShaderFlatGen::fragment_shader_source_ =
"#version 150\n"
"out vec4 fragColor;\n"
"uniform vec4 front_color;\n"
"uniform vec4 back_color;\n"
"uniform vec4 ambiant_color;\n"
"uniform vec3 lightPosition;\n"
"uniform bool cull_back_face;\n"
"in vec3 pos;\n"
"void main()\n"
"{\n"
"	vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));\n"
"	vec3 L = normalize(lightPosition-pos);\n"
"	float lambert = dot(N,L);\n"
"	if (gl_FrontFacing)\n"
"		fragColor = vec4(ambiant_color.rgb+lambert*front_color.rgb, front_color.a);\n"
//"		fragColor = ambiant_color+lambert*front_color;\n"
"	else\n"
"		if (cull_back_face) discard;\n"
//"		else fragColor = ambiant_color+lambert*back_color;\n"
"		else fragColor = vec4(ambiant_color.rgb+lambert*back_color.rgb, back_color.a);\n"
"}\n";

const char* ShaderFlatGen::vertex_shader_source2_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"in vec3 vertex_col;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"out vec3 pos;\n"
"out vec3 col;\n"
"void main()\n"
"{\n"
"	vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);\n"
"	pos = pos4.xyz;\n"
"	col = vertex_col;\n"
"   gl_Position = projection_matrix * pos4;\n"
"}\n";

const char* ShaderFlatGen::fragment_shader_source2_ =
"#version 150\n"
"out vec4 fragColor;\n"
"uniform vec4 ambiant_color;\n"
"uniform vec3 lightPosition;\n"
"uniform bool cull_back_face;\n"
"in vec3 pos;\n"
"in vec3 col;\n"
"void main()\n"
"{\n"
"	vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));\n"
"	vec3 L = normalize(lightPosition-pos);\n"
"	float lambert = dot(N,L);\n"
"	if (gl_FrontFacing)\n"
"		fragColor = ambiant_color+vec4(lambert*col,1.0);\n"
"	else\n"
"		if (cull_back_face) discard;\n"
"		else fragColor = ambiant_color-vec4(lambert*col,1.0);\n"
"}\n";

ShaderFlatGen::ShaderFlatGen(bool color_per_vertex)
{
	if (color_per_vertex)
	{
		prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source2_);
		prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source2_);
		prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
		prg_.bindAttributeLocation("vertex_col", ATTRIB_COLOR);
		prg_.link();
		get_matrices_uniforms();
	}
	else
	{
		prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
		prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
		prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
		prg_.link();
		get_matrices_uniforms();
	}
	unif_front_color_ = prg_.uniformLocation("front_color");
	unif_back_color_ = prg_.uniformLocation("back_color");
	unif_ambiant_color_ = prg_.uniformLocation("ambiant_color");
	unif_light_position_ = prg_.uniformLocation("lightPosition");
	unif_bf_culling_ = prg_.uniformLocation("cull_back_face");
}

void ShaderFlatGen::set_light_position(const QVector3D& l)
{
	prg_.setUniformValue(unif_light_position_, l);
}

void ShaderFlatGen::set_local_light_position(const QVector3D& l, const QMatrix4x4& view_matrix)
{
	QVector4D loc4 = view_matrix.map(QVector4D(l, 1.0));
	prg_.setUniformValue(unif_light_position_, QVector3D(loc4) / loc4.w());
}

void ShaderFlatGen::set_front_color(const QColor& rgb)
{
	if (unif_front_color_ >= 0)
		prg_.setUniformValue(unif_front_color_, rgb);
}

void ShaderFlatGen::set_back_color(const QColor& rgb)
{
	if (unif_back_color_ >= 0)
		prg_.setUniformValue(unif_back_color_, rgb);
}

void ShaderFlatGen::set_ambiant_color(const QColor& rgb)
{
	prg_.setUniformValue(unif_ambiant_color_, rgb);
}

void ShaderFlatGen::set_bf_culling(bool cull)
{
	prg_.setUniformValue(unif_bf_culling_, cull);
}

template class CGOGN_RENDERING_EXPORT ShaderFlatTpl<false>;
template class CGOGN_RENDERING_EXPORT ShaderFlatTpl<true>;
template class CGOGN_RENDERING_EXPORT ShaderParamFlat<false>;
template class CGOGN_RENDERING_EXPORT ShaderParamFlat<true>;

} // namespace rendering

} // namespace cgogn
