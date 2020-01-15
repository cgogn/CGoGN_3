/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
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

#include <cgogn/rendering/shaders/shader_explode_volumes.h>

namespace cgogn
{

namespace rendering
{

static const char* vertex_shader_source = "#version 150\n"
										  "uniform mat4 model_view_matrix;\n"
										  "in vec3 vertex_pos;\n"
										  "out vec3 pos;\n"
										  "void main()\n"
										  "{\n"
										  "	pos = vertex_pos;\n"
										  "   gl_Position = model_view_matrix * vec4(vertex_pos,1.0);\n"
										  "}\n";

static const char* geometry_shader_source =
	"#version 150\n"
	"layout (lines_adjacency) in;\n"
	"layout (triangle_strip, max_vertices=3) out;\n"
	"in vec3 pos[];\n"
	"out vec3 color_f;\n"
	"uniform mat4 projection_matrix;\n"
	"uniform float explode_vol;\n"
	"uniform vec3 light_position;\n"
	"uniform vec4 color;\n"
	"uniform vec4 plane_clip;\n"
	"uniform vec4 plane_clip2;\n"
	"void main()\n"
	"{\n"
	"	float d = dot(plane_clip,vec4(pos[0],1.0));\n"
	"	float d2 = dot(plane_clip2,vec4(pos[0],1.0));\n"
	"	if ((d<=0.0)&&(d2<=0.0))\n"
	"	{\n"
	"		vec3 v1 = gl_in[2].gl_Position.xyz - gl_in[1].gl_Position.xyz;\n"
	"		vec3 v2 = gl_in[3].gl_Position.xyz - gl_in[1].gl_Position.xyz;\n"
	"		vec3 N  = normalize(cross(v1,v2));\n"
	"		vec4 face_center =  gl_in[1].gl_Position;\n"
	"		vec3 L =  normalize (light_position - face_center.xyz);\n"
	"		float lambertTerm = abs(dot(N,L));\n"
	"		for (int i=1; i<=3; i++)\n"
	"		{\n"
	"			vec4 Q = explode_vol *  gl_in[i].gl_Position  + (1.0-explode_vol) * gl_in[0].gl_Position;\n"
	"			gl_Position = projection_matrix * Q;\n"
	"			color_f = color.rgb * lambertTerm;\n"
	"			EmitVertex();\n"
	"		}\n"
	"		EndPrimitive();\n"
	"	}\n"
	"}\n";

static const char* fragment_shader_source = "#version 150\n"
											"in vec3 color_f;\n"
											"uniform vec4 color;\n"
											"out vec4 fragColor;\n"
											"void main()\n"
											"{\n"
											"	if (!gl_FrontFacing) discard;\n"
											"   fragColor = vec4(color_f,color.a);\n"
											"}\n";

ShaderExplodeVolumes* ShaderExplodeVolumes::instance_ = nullptr;

ShaderExplodeVolumes::ShaderExplodeVolumes()
{
	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_pos");
	add_uniforms("color", "light_position", "explode_vol", "plane_clip", "plane_clip2");
}

} // namespace rendering
} // namespace cgogn
