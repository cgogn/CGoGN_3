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

#include <cgogn/rendering/shaders/shader_explode_volumes_scalar.h>

namespace cgogn
{

namespace rendering
{

ShaderExplodeVolumesScalar* ShaderExplodeVolumesScalar::instance_ = nullptr;

ShaderExplodeVolumesScalar::ShaderExplodeVolumesScalar()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		uniform usamplerBuffer vertex_ind;
		uniform samplerBuffer vertex_position;
		uniform samplerBuffer volume_center;
		uniform samplerBuffer volume_scalar;

		uniform float explode;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;

		out vec3 position;
		flat out vec3 color;

		//_insert_colormap_function_here

		void main()
		{
			int ind_v = int(texelFetch(vertex_ind, 4 * gl_InstanceID + gl_VertexID).r);
			int ind_c = int(texelFetch(vertex_ind, 4 * gl_InstanceID + 3).r);

			vec3 position_in = texelFetch(vertex_position, ind_v).rgb;
			vec3 center = texelFetch(volume_center, ind_c).rgb;
			float value = transform_value(texelFetch(volume_scalar, ind_c).r);
			color = value2color(value);

			float d = dot(plane_clip, vec4(center, 1.0));
			float d2 = dot(plane_clip2, vec4(center, 1.0));
			if (d <= 0.0 && d2 <= 0.0)
			{
				vec3 explode_position = mix(center, position_in, explode);
				vec4 position4 = model_view_matrix * vec4(explode_position, 1);
				position = position4.xyz;
				gl_Position = projection_matrix * position4;
			}
			else
				gl_Position = vec4(0.0, 0.0, 0.0, 1.0); // check
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform vec3 light_position;
		
		in vec3 position;
		flat in vec3 color;

		out vec4 frag_out;

		void main()
		{
			vec3 N = normalize(cross(dFdx(position), dFdy(position)));
			vec3 L = normalize(light_position - position);
			float lambert = 0.2 + 0.8 * (max(0.0, dot(N, L)));
			frag_out = vec4(lambert * color, 1.0);
		}
	)";

	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_function_here"), shader_function::ColorMap::source);
	load(v_src, fragment_shader_source);
	get_uniforms("vertex_ind", "vertex_position", "volume_center", "volume_scalar", "light_position", "explode",
				 "plane_clip", "plane_clip2", shader_function::ColorMap::uniform_names[0],
				 shader_function::ColorMap::uniform_names[1], shader_function::ColorMap::uniform_names[2],
				 shader_function::ColorMap::uniform_names[3]);

	nb_attributes_ = 3;
}

void ShaderParamExplodeVolumesScalar::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, 12, 13, light_position_, explode_, plane_clip_, plane_clip2_,
								 color_map_.color_map_, color_map_.expansion_, color_map_.min_value_,
								 color_map_.max_value_);
}

void ShaderParamExplodeVolumesScalar::bind_texture_buffers()
{
	vbos_[VERTEX_POSITION]->bind_texture_buffer(11);
	vbos_[VOLUME_CENTER]->bind_texture_buffer(12);
	vbos_[VOLUME_SCALAR]->bind_texture_buffer(13);
}

void ShaderParamExplodeVolumesScalar::release_texture_buffers()
{
	vbos_[VERTEX_POSITION]->release_texture_buffer(11);
	vbos_[VOLUME_CENTER]->release_texture_buffer(12);
	vbos_[VOLUME_SCALAR]->release_texture_buffer(13);
}

} // namespace rendering

} // namespace cgogn
