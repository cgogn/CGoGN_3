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

#include <cgogn/rendering/shaders/shader_histo.h>

namespace cgogn
{

namespace rendering
{

ShaderHisto* ShaderHisto::instance_ = nullptr;

ShaderHisto::ShaderHisto()
{
	const char* vertex_shader_source = "#version 150\n"
									   "uniform sampler2D TU;\n"
									   "uniform int width;\n"
									   "uniform float k;\n"
									   "void main()\n"
									   "{\n"
									   "	ivec2 tc  = ivec2(gl_VertexID%width, gl_VertexID/width);\n"
									   "	float p = texelFetch(TU,tc,0).r * k;\n"
									   "	p = p*2.0 - 1.0;\n"
									   "   gl_Position = vec4(p,0.0,0.0,1.0);\n"
									   "}\n";
	const char* fragment_shader_source = "#version 150\n"
										 "out float frag;\n"
										 "void main()\n"
										 "{\n"
										 "	frag = 1.0;\n"
										 "}\n";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("TU", "width", "k");
}

} // namespace rendering
} // namespace cgogn
