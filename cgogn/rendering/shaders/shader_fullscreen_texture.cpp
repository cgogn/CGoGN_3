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

#include <cgogn/rendering/shaders/shader_fullscreen_texture.h>

namespace cgogn
{

namespace rendering
{

ShaderFullScreenTexture* ShaderFullScreenTexture::instance_ = nullptr;

ShaderFullScreenTexture::ShaderFullScreenTexture()
{
	const char* vertex_shader_source = R"(
		#version 150
		out vec2 tc;
		
		void main()
		{
			vec2 p = 2.0 * vec2(gl_VertexID % 2, gl_VertexID / 2);
			tc = p;
			p = 2.0 * p - 1.0;
			gl_Position = vec4(p, 0.0, 1.0);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform sampler2D texture_unit;
		uniform float alpha;

		in vec2 tc;
		out vec4 frag_out;
		
		void main()
		{
			frag_out = vec4(texture(texture_unit, tc).rgb, alpha);
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("texture_unit", "alpha");
}

void ShaderParamFullScreenTexture::set_uniforms()
{
	shader_->set_uniforms_values(texture_->bind(unit_), alpha_);
}

} // namespace rendering

} // namespace cgogn
