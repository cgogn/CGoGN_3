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

#include <cgogn/rendering/shaders/shader_round_point.h>

#include <iostream>

namespace cgogn
{

namespace rendering
{

static const char* vertex_shader_source =
"#version 150\n"
"in vec3 vertex_pos;\n"
"void main()\n"
"{\n"
"   gl_Position =  vec4(vertex_pos,1.0);\n"
"}\n";

static const char* geometry_shader_source =
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

static const char* fragment_shader_source =
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

ShaderRoundPoint* ShaderRoundPoint::instance_ = nullptr;

ShaderRoundPoint::ShaderRoundPoint()
{
	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source,
		 "vertex_pos");
	add_uniforms("color", "pointSizes", "plane_clip", "plane_clip2");
}

} // namespace rendering

} // namespace cgogn
