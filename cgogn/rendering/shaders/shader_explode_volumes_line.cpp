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

#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>

namespace cgogn
{

namespace rendering
{
static const char* vertex_shader_source =
	R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		uniform samplerBuffer center_volume;
		uniform samplerBuffer pos_vertex;
		uniform usamplerBuffer edge_indices;

		uniform float explode;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;


		void main()
		{
			int ind_v = int(texelFetch(edge_indices,3*gl_InstanceID+gl_VertexID).r);
			int ind_c = int(texelFetch(edge_indices,3*gl_InstanceID+2).r);

			vec3 position_in = texelFetch(pos_vertex, ind_v).rgb;
			vec3 center = texelFetch(center_volume, ind_c).rgb;

			float d = dot(plane_clip, vec4(center,1));
			float d2 = dot(plane_clip2,vec4(center,1));
			if ((d<=0.0)&&(d2<=0.0))
			{
				vec3 P = mix(center, position_in, explode);
				vec4 Po4 = model_view_matrix * vec4(P,1);
				gl_Position = projection_matrix * Po4;
			}
			else
			{
				gl_Position = vec4(0,0,0,1); // check
			}
		}
		)";

static const char* fragment_shader_source =
	R"(#version 330
		uniform vec4 color;
		out vec4 frag_out;
		void main()
		{
			frag_out = color;
		}
		)";

ShaderExplodeVolumesLine* ShaderExplodeVolumesLine::instance_ = nullptr;

ShaderExplodeVolumesLine::ShaderExplodeVolumesLine()
{
	load2_bind(vertex_shader_source, fragment_shader_source);
	add_uniforms("edge_indices", "pos_vertex", "center_volume", "color", "explode", "plane_clip", "plane_clip2");
	this->nb_attributes_ = 2;
}

void ShaderParamExplodeVolumesLine::set_uniforms()
{
	shader_->set_uniforms_values(10, vbos_[POS]->bind_tb(11), vbos_[CENTER]->bind_tb(12), color_, explode_, plane_clip_,
								 plane_clip2_);
}

} // namespace rendering

} // namespace cgogn
