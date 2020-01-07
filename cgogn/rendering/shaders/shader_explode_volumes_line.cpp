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

#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>

namespace cgogn
{

namespace rendering
{
static const char* vertex_shader_source = "#version 150\n"
										  "in vec3 vertex_pos;\n"
										  "void main()\n"
										  "{\n"
										  "   gl_Position = vec4(vertex_pos,1.0);\n"
										  "}\n";

static const char* geometry_shader_source =
	"#version 150\n"
	"layout (triangles) in;\n"
	"layout (line_strip, max_vertices=2) out;\n"
	"uniform mat4 projection_matrix;\n"
	"uniform mat4 model_view_matrix;\n"
	"uniform float explode_vol;\n"
	"uniform vec4 plane_clip;\n"
	"uniform vec4 plane_clip2;\n"
	"void main()\n"
	"{\n"
	"	float d = dot(plane_clip,gl_in[0].gl_Position);\n"
	"	float d2 = dot(plane_clip2,gl_in[0].gl_Position);\n"
	"	if ((d<=0.0)&&(d2<=0.0))\n"
	"	{\n"
	"		for (int i=1; i<=2; i++)\n"
	"		{\n"
	"			vec4 Q = explode_vol *  gl_in[i].gl_Position  + (1.0-explode_vol) * gl_in[0].gl_Position;\n"
	"			gl_Position = projection_matrix * model_view_matrix *  Q;\n"
	"			EmitVertex();\n"
	"		}\n"
	"		EndPrimitive();\n"
	"	}\n"
	"}\n";

static const char* fragment_shader_source = "#version 150\n"
											"uniform vec4 color;\n"
											"out vec4 fragColor;\n"
											"void main()\n"
											"{\n"
											"   fragColor = color;\n"
											"}\n";

ShaderExplodeVolumesLine* ShaderExplodeVolumesLine::instance_ = nullptr;

ShaderExplodeVolumesLine::ShaderExplodeVolumesLine()
{
	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_pos");
	add_uniforms("color", "explode_vol", "plane_clip", "plane_clip2");
}

} // namespace rendering

} // namespace cgogn
