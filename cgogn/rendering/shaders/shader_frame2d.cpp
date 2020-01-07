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

#include <cgogn/rendering/shaders/shader_frame2d.h>

namespace cgogn
{

namespace rendering
{

ShaderFrame2d* ShaderFrame2d::instance_ = nullptr;

ShaderFrame2d::ShaderFrame2d()
{
	const char* vertex_shader_source =
		"#version 150\n"
		"const float vertex_pos[]=float[](-1,-1,1, -1,-1,0, 1,-1,1, 1,-1,0, 1,1,1, 1,1,0, -1,1,1, -1,1,0, -1,-1,1, "
		"-1,-1,0);\n"
		"uniform float w;\n"
		"uniform float h;\n"
		"uniform float sz;\n"
		"out float alpha;\n"
		"void main()\n"
		"{\n"
		"	vec2 mu = vec2(1.0,1.0) - 2.0 * sz* float((gl_VertexID)%2)/vec2(w,h);\n"
		"	vec2 P = vec2(vertex_pos[3*gl_VertexID],vertex_pos[3*gl_VertexID+1]) * mu;\n"
		"   gl_Position = vec4(P,0.0,1.0);\n"
		"	alpha = vertex_pos[3*gl_VertexID +2];\n"
		"}\n";

	const char* fragment_shader_source = "#version 150\n"
										 "in float alpha;\n"
										 "out vec4 frag;\n"
										 "uniform vec4 color;\n"
										 "void main()\n"
										 "{\n"
										 "	frag = vec4(color.rgb,alpha);\n"
										 "}\n";

	load(vertex_shader_source, fragment_shader_source);
	add_uniforms("color", "sz", "w", "h");
}

} // namespace rendering
} // namespace cgogn
