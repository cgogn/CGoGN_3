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

#define CGOGN_RENDER_SHADERS_OBJ_FLAT_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_obj_flat_texture.h>

namespace cgogn
{

namespace rendering
{

ShaderObjFlatTexture* ShaderObjFlatTexture::instance_ = nullptr;

ShaderObjFlatTexture::ShaderObjFlatTexture()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform mat3 normal_matrix;

		uniform usamplerBuffer position_ind;
		uniform usamplerBuffer tex_coord_ind;
		uniform samplerBuffer vertex_position;
		uniform samplerBuffer vertex_tc;

		out vec3 position;
		out vec2 tc;

		void main()
		{
			int ind_v = int(texelFetch(position_ind, 3 * gl_InstanceID + gl_VertexID).r);
			vec3 position_in = texelFetch(vertex_position, ind_v).rgb;

			int ind_tc = int(texelFetch(tex_coord_ind, 3 * gl_InstanceID + gl_VertexID).r);
			tc = texelFetch(vertex_tc, ind_tc).rg;

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
		in vec2 tc;

		out vec3 frag_out;

		void main()
		{
			vec3 N = normalize(cross(dFdx(position), dFdy(position)));
			vec3 L = normalize(light_position - position);

			vec3 color = vec3(tc,0);
			float lambert = 0.25 + 0.75 * max(0.0,dot(N, L));
			frag_out = color*lambert;
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("position_ind", "tex_coord_ind", "vertex_position", "vertex_tc",  "light_position");
}

void ShaderParamObjFlatTexture::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, 12, 13, light_position_);
}

void ShaderParamObjFlatTexture::bind_texture_buffers()
{
	vbos_[VERTEX_POSITION]->bind_texture_buffer(12);
	vbos_[VERTEX_TC]->bind_texture_buffer(13);
}

void ShaderParamObjFlatTexture::release_texture_buffers()
{
	vbos_[VERTEX_POSITION]->release_texture_buffer(12);
	vbos_[VERTEX_TC]->release_texture_buffer(13);
}

} // namespace rendering

} // namespace cgogn
