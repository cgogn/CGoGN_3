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

#include <cgogn/rendering/shaders/shader_no_illum_color_per_face.h>

namespace cgogn
{

namespace rendering
{

ShaderNoIllumColorPerFace* ShaderNoIllumColorPerFace::instance_ = nullptr;

ShaderNoIllumColorPerFace::ShaderNoIllumColorPerFace()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform usamplerBuffer tri_ind;
		uniform usamplerBuffer face_emb;
		uniform samplerBuffer vertex_position;
		uniform samplerBuffer face_color;

		flat out vec3 N;
		flat out vec3 color;

		out vec3 A;
		
		void main()
		{
			int tri = int(texelFetch(tri_ind, gl_InstanceID).r);
			int i_c = int(texelFetch(face_emb, gl_InstanceID).r);
			color = texelFetch(face_color, i_c).rgb;
			int vid = gl_VertexID;
			int ind_a = int(texelFetch(tri_ind, 3*gl_InstanceID+vid).r);
			A = (model_view_matrix * vec4(texelFetch(vertex_position, ind_a).rgb,1.0)).xyz;
			vid  = (vid+1)%3;
			int ind_b = int(texelFetch(tri_ind, 3*gl_InstanceID+vid).r);
			vec3 B = (model_view_matrix * vec4(texelFetch(vertex_position, ind_b).rgb,1.0)).xyz;
			vid  = (vid+1)%3;
			int ind_c = int(texelFetch(tri_ind, 3*gl_InstanceID+vid).r);
			vec3 C = (model_view_matrix * vec4(texelFetch(vertex_position, ind_c).rgb,1.0)).xyz;
			N = normalize(cross(B-A,C-A));
			gl_Position = projection_matrix*vec4(A,1);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform bool double_side;
		
		flat in vec3 color;
		
		out vec3 frag_out;

		void main()
		{
			if (double_side || gl_FrontFacing)
				frag_out = color;
			else
				discard;
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("tri_ind", "face_emb", "vertex_position", "face_color", "double_side");

	nb_attributes_ = 2;
}

void ShaderParamNoIllumColorPerFace::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, 12, 13, double_side_);
}

void ShaderParamNoIllumScalarPerFace::bind_texture_buffers()
{
	vbos_[VERTEX_POSITION]->bind_texture_buffer(12);
	vbos_[FACE_COLOR]->bind_texture_buffer(13);
}

void ShaderParamNoIllumScalarPerFace::release_texture_buffers()
{
	vbos_[VERTEX_POSITION]->release_texture_buffer(12);
	vbos_[FACE_COLOR]->release_texture_buffer(13);
}

} // namespace rendering

} // namespace cgogn
