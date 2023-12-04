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

#include <cgogn/rendering/shaders/shader_mesh_2d_edges.h>

namespace cgogn
{

namespace rendering
{

ShaderMesh2DEdges* ShaderMesh2DEdges::instance_ = nullptr;

ShaderMesh2DEdges::ShaderMesh2DEdges()
{
char const* vertex_shader_source = R"(
	#version 330
	uniform vec2 ratio;
	uniform usamplerBuffer edge_indices;
	uniform samplerBuffer vertex_tc;

	void main()
	{
		int ind = gl_VertexID;
		int ind_tc = int(texelFetch(edge_indices, ind).r);
		vec2 P2 = vec2(2.0 * texelFetch(vertex_tc, ind_tc).rg - 1.0) * ratio;
		gl_Position = vec4(P2.xy, 0.5, 1);
	}
	)";

	char const* fragment_shader_source = R"(
		#version 330
		uniform vec4 color;
		
		out vec4 frag_out;

		void main()
		{
			frag_out = color;
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("edge_indices", "vertex_tc", "ratio", "color");
}

void ShaderParamMesh2DEdges::set_uniforms()
{

	shader_->set_uniforms_values(10, 11, ratio_, color_);
}

void ShaderParamMesh2DEdges::bind_texture_buffers()
{
	vbos_[VERTEX_TC]->bind_texture_buffer(11);
}

void ShaderParamMesh2DEdges::release_texture_buffers()
{
	vbos_[VERTEX_TC]->release_texture_buffer(11);
}


} // namespace rendering

} // namespace cgogn
