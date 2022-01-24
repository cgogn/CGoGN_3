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

#include <cgogn/rendering/shaders/texture_to_vbo.h>

namespace cgogn
{

namespace rendering
{

ShaderTEX2VBO* ShaderTEX2VBO::instance_ = nullptr;

ShaderTEX2VBO::ShaderTEX2VBO()
{
	const char* vertex_shader_source =
			R"(
			#version 330
			uniform sampler2D TUin;
			uniform bool normalization;
			out vec3 vbo_out;

			void main()
			{
				int w = textureSize(TUn,0).x;
				ivec2 icoord = ivec2(gl_VertexID%w,gl_VertexID/w);
				if (normalization)
					vbo_out = normalize(texelFetch(TUn,icoord,0).rgb);
				else
					vbo_out = texelFetch(TUn,icoord,0).rgb;
			}
			)";


	load_tfb1_bind(vertex_shader_source, {"vbo_out"});

	get_uniforms("TUin","normalization");
}

} // namespace rendering

} // namespace cgogn
