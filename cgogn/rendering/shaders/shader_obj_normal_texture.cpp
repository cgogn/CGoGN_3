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

#define CGOGN_RENDER_SHADERS_OBJ_NORMAL_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_obj_normal_texture.h>

namespace cgogn
{

namespace rendering
{

ShaderObjNormalTexture* ShaderObjNormalTexture::instance_ = nullptr;

ShaderObjNormalTexture::ShaderObjNormalTexture()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform mat3 normal_matrix;

		uniform usamplerBuffer position_ind;
		uniform usamplerBuffer tex_coord_ind;
		uniform usamplerBuffer normal_ind;
		uniform samplerBuffer vertex_position;
		uniform samplerBuffer vertex_tc;
		uniform samplerBuffer vertex_normal;

		out vec3 position;
		out vec2 tc;
		out vec2 normal;

		void main()
		{
			int ind =  3 * gl_InstanceID + gl_VertexID
			int ind_v = int(texelFetch(position_ind,ind).r);
			vec3 position_in = texelFetch(vertex_position, ind_v).rgb;

			int ind_tc = int(texelFetch(tex_coord_ind, ind).r);
			tc = texelFetch(vertex_tc, ind_tc).rg;

			int ind_no = int(texelFetch(normal_ind, ind).r);
			normal = normal_matrix*texelFetch(vertex_normal, ind_no).rgb;

			vec4 position4 = model_view_matrix * vec4(position_in, 1.0);
			position = position4.xyz;
			gl_Position = projection_matrix * position4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330

		uniform sampler2D texture_img_unit;
		uniform vec3 light_position;
		uniform bool draw_param;

		in vec3 position;
		in vec2 tc;
		in vec3 normal
		out vec3 frag_out;

		void main()
		{
			vec3 N = normalize(normal);
			vec3 L = normalize(light_position - position);
			float lambert = 0.25 + 0.75 * max(0.0,dot(N, L));
			frag_out = (draw_param) ? N*0.5+0.5 : texture(texture_img_unit,vec2(tc.x,1.0-tc.y)).rgb*lambert;
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("position_ind", "tex_coord_ind", "vertex_position", "vertex_tc", "texture_img_unit", "light_position",
				 "draw_param");
}

void ShaderParamObjNormalTexture::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, 12, 13, texture_->bind(0), light_position_, draw_param_);
}

void ShaderParamObjNormalTexture::bind_texture_buffers()
{
	vbos_[VERTEX_POSITION]->bind_texture_buffer(12);
	vbos_[VERTEX_TC]->bind_texture_buffer(13);
}

void ShaderParamObjNormalTexture::release_texture_buffers()
{
	vbos_[VERTEX_POSITION]->release_texture_buffer(12);
	vbos_[VERTEX_TC]->release_texture_buffer(13);
}

} // namespace rendering

} // namespace cgogn
