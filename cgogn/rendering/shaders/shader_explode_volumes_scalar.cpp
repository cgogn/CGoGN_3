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

static const char* vertex_shader_source =
	R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		uniform samplerBuffer pos_vertex;
		uniform samplerBuffer center_volume;
		uniform samplerBuffer scalar_volume;
		uniform usamplerBuffer tri_indices;

		uniform float explode;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;

		flat out vec3 color;
		out vec3 Po;

//_insert_colormap_function_here

		void main()
		{
			int ind_v = int(texelFetch(tri_indices,4*gl_InstanceID+gl_VertexID).r);
			int ind_c = int(texelFetch(tri_indices,4*gl_InstanceID+3).r);

			vec3 position_in = texelFetch(pos_vertex, ind_v).rgb;
			vec3 center = texelFetch(center_volume, ind_c).rgb;
			color = scalar2color(texelFetch(scalar_volume, ind_c).r);

			float d = dot(plane_clip, vec4(center,1));
			float d2 = dot(plane_clip2,vec4(center,1));
			if ((d<=0.0)&&(d2<=0.0))
			{
				vec3 P = mix(center, position_in, explode);
				vec4 Po4 = model_view_matrix * vec4(P,1);
				Po = Po4.xyz;
				gl_Position = projection_matrix * Po4;
			}
			else
			{
				gl_Position = vec4(0,0,0,1); // check
			}
		}
		)";

static const char* fragment_shader_source =
	R"(
		#version 330
		out vec3 frag_out;
		uniform vec3 light_position;
		flat in vec3 color;
		in vec3 Po;
		void main()
		{
			vec3 N = normalize(cross(dFdx(Po),dFdy(Po)));
			vec3 L = normalize(light_position-Po);
			float lambert = 0.2+0.8*(max(0.0,dot(N,L)));
			frag_out = lambert * color;
		}
		)";

ShaderExplodeVolumesScalar* ShaderExplodeVolumesScalar::instance_ = nullptr;

ShaderExplodeVolumesScalar::ShaderExplodeVolumesScalar()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"), shader_function::ColorMap::source);
	load2_bind(v_src, fragment_shader_source);
	add_uniforms("tri_indices", "pos_vertex", "center_volume", "scalar_volume", "light_position", "explode",
				 "plane_clip", "plane_clip2", shader_function::ColorMap::name[0], shader_function::ColorMap::name[1],
				 shader_function::ColorMap::name[2], shader_function::ColorMap::name[3]);
	this->nb_attributes_ = 3;
}

void ShaderParamExplodeVolumesScalar::set_uniforms()
{
	shader_->set_uniforms_values(10, vbos_[POS]->bind_tb(11), vbos_[CENTER]->bind_tb(12), vbos_[SCALAR]->bind_tb(13),
								 light_pos_, explode_, plane_clip_, plane_clip2_, cm_.color_map_, cm_.expansion_,
								 cm_.min_value_, cm_.max_value_);
}

} // namespace rendering

} // namespace cgogn
