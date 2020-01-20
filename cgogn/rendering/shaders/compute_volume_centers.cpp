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

#include <cgogn/rendering/shaders/compute_volume_centers.h>

namespace cgogn
{

namespace rendering
{
static const char* vertex_shader_source1 =
		R"(
		#version 330
		uniform mat3 tex_matrix;
		uniform usamplerBuffer vertex_volume;
		uniform samplerBuffer pos_vertex;
		out vec4 P;
		void main()
		{
			int vid = gl_VertexID;
			int ind_v = int(texelFetch( vid++).r);
			uint ind_c = texelFetch(vid).r;
			vec2 coord_c = tex_matrix * vec2(float(ind_c%1024),float(ind_c/1024),1);
			gl_Position = vec4(coord_c,0,1);
			P.xyz = texelFetch(pos_vertex, ind_v).rgb;
			P.w = 1.0;
		}
		)";

static const char* fragment_shader_source1 =
		R"(
		#version 330
		out vec4 fragOout;
		in vec4 P;
		void main()
		{
			fragOut = P;
		};
		)";


static const char* vertex_shader_source2 =
		R"(
		#version 330
		uniform sampler2D TUin;
		out vec3 vbo_out;

		const int w = 1024;

		void main()
		{
			ivec2 icoord = ivec2(gl_VertexID%w,gl_VertexID/w);
			vec4 P4 = texelFetch(TUn,icoord,0)
			vbo_out = P4.xyz / P4.w;
		}
		)";



ShaderComputeCenter1* ShaderComputeCenter1::instance_ = nullptr;

ShaderComputeCenter1::ShaderComputeCenter1()
{
	load2_bind(vertex_shader_source1, fragment_shader_source1);
	add_uniforms("tex_matrix","pos_vertex");
}


ShaderComputeCenter2* ShaderComputeCenter2::instance_ = nullptr;

ShaderComputeCenter2::ShaderComputeCenter2()
{
	load_tfb1_bind(vertex_shader_source2, {"vbo_out"});

}


} // namespace rendering

} // namespace cgogn
