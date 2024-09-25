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

#include <cgogn/rendering/shaders/shader_texture.h>

namespace cgogn
{

namespace rendering
{

ShaderTexture* ShaderTexture::instance_ = nullptr;

ShaderTexture::ShaderTexture()
{
	const char* vertex_shader_source =
		"#version 150\n"
		"in vec3 vertex_pos;\n"
		"in vec2 vertex_tc;\n"
		"uniform mat4 projection_matrix;\n"
		"uniform mat4 model_view_matrix;\n"
		"out vec2 tc;\n"
		"void main()\n"
		"{\n"
		"	tc = vertex_tc;\n"
		"   gl_Position = projection_matrix * model_view_matrix * vec4(vertex_pos,1.0);\n"
		"}\n";
	const char* fragment_shader_source = "#version 150\n"
										 "out vec4 frag_color;\n"
										 "uniform sampler2D texture_unit;\n"
										 "in vec2 tc;\n"
										 "void main()\n"
										 "{\n"
										 "	frag_color = texture(texture_unit,tc);\n"
										 "}\n";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos", "vertex_tc");
	get_uniforms("texture_unit");
}

} // namespace rendering

} // namespace cgogn
