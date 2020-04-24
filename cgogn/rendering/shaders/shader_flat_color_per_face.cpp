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
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		uniform usamplerBuffer vertex_ind;
		uniform usamplerBuffer tri_ind;
		uniform samplerBuffer pos_vertex;
		uniform samplerBuffer color_tri;
		
		out vec3 pos;
		flat out vec3 color;

		void main()
		{
			int ind_v = int(texelFetch(vertex_ind, 3*gl_InstanceID+gl_VertexID).r);
			vec3 position_in = texelFetch(pos_vertex, ind_v).rgb;

			int ind_t = int(texelFetch(tri_ind, int(gl_InstanceID)).r);
			color = texelFetch(color_tri, ind_t).rgb;

			vec4 pos4 = model_view_matrix * vec4(position_in,1.0);
			pos = pos4.xyz;
			gl_Position = projection_matrix * pos4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform vec4 ambiant_color;
		uniform vec3 light_position;
		uniform bool double_side;
		
		in vec3 pos;
		flat in vec3 color;

		out vec3 fragColor;

		void main()
		{
			vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));
			vec3 L = normalize(light_position-pos);
			float lambert = dot(N,L);
			if (double_side || gl_FrontFacing)
				fragColor = ambiant_color.rgb+lambert*color;
			else
				discard;
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source, "");

	add_uniforms("vertex_ind", "tri_ind", "pos_vertex", "color_tri", "ambiant_color", "light_position", "double_side");

	nb_attributes_ = 2;
}

void ShaderParamFlatColorPerFace::set_uniforms()
{
	vbos_[0]->bind_tb(12);
	vbos_[1]->bind_tb(13);
	shader_->set_uniforms_values(10, 11, 12, 13, ambiant_color_, light_position_, double_side_);
}

} // namespace rendering

} // namespace cgogn
