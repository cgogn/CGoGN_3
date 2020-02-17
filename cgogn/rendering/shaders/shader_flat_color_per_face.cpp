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

#include <cgogn/rendering/shaders/shader_flat_color_per_face.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatColorPerFace* ShaderFlatColorPerFace::instance_ = nullptr;

ShaderFlatColorPerFace::ShaderFlatColorPerFace()
{
	const char* vertex_shader_source =
			R"(
			#version 330
			uniform mat4 projection_matrix;
			uniform mat4 model_view_matrix;
			uniform usamplerBuffer tri_ind;
			uniform usamplerBuffer tri_emb;
			uniform samplerBuffer pos_vertex;
			uniform samplerBuffer color_tri;
			out vec3 A;
			flat out vec3 N;
			flat out color;
			void main()
			{
				int tri = int(texelFetch(tri_ind, int(gl_InstanceID)).r);
				int color_ind = int(texelFetch(tri_emb, int(gl_InstanceID)).r);
				color = texelFetch(color_tri, color_i).rgb
				int vid = gl_VertexID;
				int ind_a = int(texelFetch(tri_ind, 3*int(gl_InstanceID)+vid).r);
				A = model_view_matrix * vec4(texelFetch(pos_vertex, ind_a).rgb)).rgb;
				vid  = (vide+1)%3
				int ind_b = int(texelFetch(tri_ind, 3*int(gl_InstanceID)+vid).r);
				vec3 B = model_view_matrix * vec4(texelFetch(pos_vertex, (ind_b)).rgb)).rgb;
				vid  = (vide+1)%3
				int ind_c = int(texelFetch(tri_ind, 3*int(gl_InstanceID)+vid).r);
				vec3 C = model_view_matrix * vec4(texelFetch(pos_vertex, ind_c).rgb)).rgb;
				N = normalize(cross(C-A,B-A));
				vec2 coord_N = (-1.0+d) + 2.0 * d * vec2(float(ind_a%1024u),float(ind_a/1024u));
				gl_Position = projection_matrix*vec4(A,1);
			}
			)";

	const char* fragment_shader_source =
		R"(#version 330
		out vec3 fragColor;
		uniform vec4 ambiant_color;
		uniform vec3 light_position;
		uniform bool double_side;
		in vec3 A;
		flat in vec3 N;
		flat in vec3 color;
		void main()
		{
			vec3 No = normalize(cross(N);
			vec3 L = normalize(light_position-A);
			float lambert = dot(No,L);
			if (double_side || gl_FrontFacing)
				fragColor = ambiant_color.rgb+lambert*color);
			else
				discard;
		}
		)";

	load2_bind(vertex_shader_source, fragment_shader_source, "");

	add_uniforms("tri_ind","tri_emb", "pos_vertex","ambiant_color", "light_position", "double_side");
}


void ShaderParamFlatColorPerFace::set_uniforms()
{
	if (vbo_pos_)
		shader_->set_uniforms_values(10,11, vbo_pos_->bind_tb(12), 1.0f/float(height_tex_));
}


} // namespace rendering

} // namespace cgogn
