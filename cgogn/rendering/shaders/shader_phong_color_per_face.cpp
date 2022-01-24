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

#define CGOGN_RENDER_SHADERS_PHONG_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_phong_color_per_face.h>

namespace cgogn
{

namespace rendering
{

ShaderPhongColorPerFace* ShaderPhongColorPerFace::instance_ = nullptr;

ShaderPhongColorPerFace::ShaderPhongColorPerFace()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform mat3 normal_matrix;

		uniform usamplerBuffer vertex_ind;
		uniform usamplerBuffer face_ind;
		uniform samplerBuffer vertex_position;
		uniform samplerBuffer vertex_normal;
		uniform samplerBuffer face_color;

		out vec3 position;
		out vec3 normal;
		flat out vec3 color;

		void main()
		{
			int ind_v = int(texelFetch(vertex_ind, 3 * gl_InstanceID + gl_VertexID).r);
			vec3 position_in = texelFetch(vertex_position, ind_v).rgb;
			vec3 normal_in = texelFetch(vertex_normal, ind_v).rgb;

			int ind_f = int(texelFetch(face_ind, gl_InstanceID).r);
			color = texelFetch(face_color, ind_f).rgb;;

			normal = normal_matrix * normal_in;
			vec4 position4 = model_view_matrix * vec4(position_in, 1.0);
			position = position4.xyz;
			gl_Position = projection_matrix * position4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform vec4 ambiant_color;
		uniform vec3 light_position;
		uniform bool double_side;
		uniform vec4 specular_color;
		uniform float specular_coef;

		in vec3 position;
		in vec3 normal;
		flat in vec3 color;

		out vec4 frag_out;

		void main()
		{
			vec3 N = normalize(normal);
			vec3 L = normalize(light_position - position);
			vec4 final_color = ambiant_color;
			if (!gl_FrontFacing)
			{
				if (!double_side)
					discard;
				N *= -1.0;
			}
			float lambert_term = clamp(dot(N, L), 0.0, 1.0);
			final_color += vec4(color, 1.0) * lambert_term;
			vec3 E = normalize(-position);
			vec3 R = reflect(-L, N);
			float specular = pow(max(dot(R,E), 0.0), specular_coef);
			final_color += specular_color * specular;
			frag_out = vec4(final_color.rgb, 1.0);
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("vertex_ind", "face_ind", "vertex_position", "vertex_normal", "face_color", "light_position",
				 "ambiant_color", "specular_color", "specular_coef", "double_side");

	nb_attributes_ = 3;
}

void ShaderParamPhongColorPerFace::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, 12, 13, 14, light_position_, ambiant_color_, specular_color_, specular_coef_,
								 double_side_);
}

void ShaderParamPhongColorPerFace::bind_texture_buffers()
{
	vbos_[VERTEX_POSITION]->bind_texture_buffer(12);
	vbos_[VERTEX_NORMAL]->bind_texture_buffer(13);
	vbos_[FACE_COLOR]->bind_texture_buffer(14);
}

void ShaderParamPhongColorPerFace::release_texture_buffers()
{
	vbos_[VERTEX_POSITION]->release_texture_buffer(12);
	vbos_[VERTEX_NORMAL]->release_texture_buffer(13);
	vbos_[FACE_COLOR]->release_texture_buffer(14);
}

} // namespace rendering

} // namespace cgogn
