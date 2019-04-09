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

#define CGOGN_RENDER_SHADERS_ROUND_POINT_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_round_point.h>


namespace cgogn
{

namespace rendering
{

const char* ShaderRoundPointGen::vertex_shader_source_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"void main()\n"
"{\n"
"   gl_Position =  vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderRoundPointGen::geometry_shader_source_ =
"#version 150\n"
"layout (points) in;\n"
"layout (triangle_strip, max_vertices=4) out;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform vec2 pointSizes;\n"
"uniform vec4 plane_clip;\n"
"uniform vec4 plane_clip2;\n"
"out vec2 local;\n"
"void main()\n"
"{\n"
"	float d = dot(plane_clip,gl_in[0].gl_Position);\n"
"	float d2 = dot(plane_clip2,gl_in[0].gl_Position);\n"
"	if ((d<=0.0)&&(d2<=0.0))\n"
"	{\n"
"		vec4 A = projection_matrix*model_view_matrix * gl_in[0].gl_Position;\n"
"		A = A/A.w;\n"
"		local = vec2(-1.0,-1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(-pointSizes[0],-pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		local = vec2(1.0,-1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(pointSizes[0],-pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		local = vec2(-1.0,1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(-pointSizes[0],pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		local = vec2(1.0,1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(pointSizes[0],pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		EndPrimitive();\n"
"	}\n"
"}\n";

const char* ShaderRoundPointGen::fragment_shader_source_ =
"#version 150\n"
"uniform vec4 color;\n"
"in vec2 local;\n"
"out vec4 fragColor;\n"
"void main()\n"
"{\n"

"	float r2 = dot(local,local);\n"
"   if (r2 > 1.0) discard;\n"
"   fragColor = vec4(color.rgb,(1.0-r2*r2));\n"
"}\n";

const char* ShaderRoundPointGen::vertex_shader_source2_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"in vec3 vertex_color;\n"
"out vec3 color_v;\n"
"void main()\n"
"{\n"
"   color_v = vertex_color;\n"
"   gl_Position = vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderRoundPointGen::geometry_shader_source2_ =
"#version 150\n"
"layout (points) in;\n"
"layout (triangle_strip, max_vertices=4) out;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform vec2 pointSizes;\n"
"uniform vec4 plane_clip;\n"
"uniform vec4 plane_clip2;\n"
"in vec3 color_v[];\n"
"out vec2 local;\n"
"out vec3 color_f;\n"
"void main()\n"
"{\n"
"	float d = dot(plane_clip,gl_in[0].gl_Position);\n"
"	float d2 = dot(plane_clip2,gl_in[0].gl_Position);\n"
"	if ((d<=0.0)&&(d2<=0.0))\n"
"	{\n"
"		vec4 A = projection_matrix*model_view_matrix * gl_in[0].gl_Position;\n"
"		A = A/A.w;\n"
"		color_f = color_v[0];\n"
"		local = vec2(-1.0,-1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(-pointSizes[0],-pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		local = vec2(1.0,-1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(pointSizes[0],-pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		local = vec2(-1.0,1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(-pointSizes[0],pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		local = vec2(1.0,1.0);\n"
"		gl_Position = vec4(A.xyz-vec3(pointSizes[0],pointSizes[1],0.0), 1.0);\n"
"		EmitVertex();\n"
"		EndPrimitive();\n"
"	}\n"
"}\n";

const char* ShaderRoundPointGen::fragment_shader_source2_ =
"#version 150\n"
"in vec2 local;\n"
"in vec3 color_f;\n"
"out vec4 fragColor;\n"
"void main()\n"
"{\n"
"	float r2 = dot(local,local);\n"
"   if (r2 > 1.0) discard;\n"
"   fragColor = vec4(color_f,(1.0-r2*r2));\n"
"}\n";

ShaderRoundPointGen::ShaderRoundPointGen(bool color_per_vertex)
{
	if (color_per_vertex)
	{
		prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source2_);
		prg_.addShaderFromSourceCode(QOpenGLShader::Geometry, geometry_shader_source2_);
		prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source2_);
		prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
		prg_.bindAttributeLocation("vertex_color", ATTRIB_COLOR);
		prg_.link();
	}
	else
	{
		prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
		prg_.addShaderFromSourceCode(QOpenGLShader::Geometry, geometry_shader_source_);
		prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
		prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
		prg_.link();
	}

	get_matrices_uniforms();
	unif_color_ = prg_.uniformLocation("color");
	unif_size_ = prg_.uniformLocation("pointSizes");
	unif_plane_clip_ = prg_.uniformLocation("plane_clip");
	unif_plane_clip2_ = prg_.uniformLocation("plane_clip2");

	set_size(3.0f);
	set_color(QColor(255, 255, 255));
	set_plane_clip(QVector4D(0,0,0,0));
	set_plane_clip2(QVector4D(0,0,0,0));
}

void ShaderRoundPointGen::set_color(const QColor& rgb)
{
	if (unif_color_ >= 0)
		prg_.setUniformValue(unif_color_, rgb);
}

void ShaderRoundPointGen::set_size(float32 wpix)
{
	QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
	int viewport[4];
	ogl->glGetIntegerv(GL_VIEWPORT, viewport);
	QSizeF wd(wpix / float32(viewport[2]), wpix / float32(viewport[3]));
	prg_.setUniformValue(unif_size_, wd);
}


void ShaderRoundPointGen::set_plane_clip(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip_, plane);
}

void ShaderRoundPointGen::set_plane_clip2(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip2_, plane);
}




template class CGOGN_RENDERING_EXPORT ShaderRoundPointTpl<false>;
template class CGOGN_RENDERING_EXPORT ShaderRoundPointTpl<true>;
template class CGOGN_RENDERING_EXPORT ShaderParamRoundPoint<false>;
template class CGOGN_RENDERING_EXPORT ShaderParamRoundPoint<true>;


} // namespace rendering

} // namespace cgogn
