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

#include <cgogn/rendering/shaders/shader_flat_scalar_per_face.h>

namespace cgogn
{

namespace rendering
{

static const char* vertex_shader_source =
	R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform usamplerBuffer tri_ind;
uniform usamplerBuffer tri_emb;
uniform samplerBuffer pos_vertex;
uniform samplerBuffer scalar_tri;
out vec3 A;
flat out vec3 N;
flat out vec3 color;

//_insert_colormap_funcion_here

void main()
{
	int tri = int(texelFetch(tri_ind, gl_InstanceID).r);
	int i_col = int(texelFetch(tri_emb, gl_InstanceID).r);
	color = scalar2color(texelFetch(scalar_tri, i_col).r);
	int vid = gl_VertexID;
	int tid = 3*gl_InstanceID;
	int ind_a = int(texelFetch(tri_ind, tid+vid).r);
	A = (model_view_matrix * vec4(texelFetch(pos_vertex, ind_a).rgb,1.0)).xyz;
	vid  = (vid+1)%3;
	int ind_b = int(texelFetch(tri_ind, tid+vid).r);
	vec3 B = (model_view_matrix * vec4(texelFetch(pos_vertex, ind_b).rgb,1.0)).xyz;
	vid  = (vid+1)%3;
	int ind_c = int(texelFetch(tri_ind, tid+vid).r);
	vec3 C = (model_view_matrix * vec4(texelFetch(pos_vertex, ind_c).rgb,1.0)).xyz;
	N = normalize(cross(B-A,C-A));
	gl_Position = projection_matrix*vec4(A,1);
}
)";

static const char* fragment_shader_source =
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
			vec3 No = normalize(N);
			vec3 L = normalize(light_position-A);
			float lambert = dot(No,L);
			if (double_side || gl_FrontFacing)
				fragColor = ambiant_color.rgb + lambert*color;
			else
				discard;
		}
		)";

ShaderFlatScalarPerFace* ShaderFlatScalarPerFace::instance_ = nullptr;

ShaderFlatScalarPerFace::ShaderFlatScalarPerFace()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"), shader_funcion::ColorMap::source);
	load2_bind(v_src, fragment_shader_source, "");

	add_uniforms("tri_ind", "tri_emb", "pos_vertex", "scalar_tri", "ambiant_color", "light_position", "double_side",
				 shader_funcion::ColorMap::name[0], shader_funcion::ColorMap::name[1],
				 shader_funcion::ColorMap::name[2], shader_funcion::ColorMap::name[3]);
	this->nb_attributes_ = 2;
}

void ShaderParamFlatScalarPerFace::set_uniforms()
{
	vbos_[0]->bind_tb(12);
	vbos_[1]->bind_tb(13);
	shader_->set_uniforms_values(10, 11, 12, 13, ambiant_color_, light_position_, double_side_, cm_.color_map_,
								 cm_.expansion_, cm_.min_value_, cm_.max_value_);
}

} // namespace rendering

} // namespace cgogn
