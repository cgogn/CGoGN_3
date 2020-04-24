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

#include <iostream>

#include <cgogn/rendering/shaders/shader_phong_scalar_per_face.h>

namespace cgogn
{

namespace rendering
{

ShaderPhongScalarPerFace* ShaderPhongScalarPerFace::instance_ = nullptr;

static const char* vertex_shader_source = R"(
	#version 330
	uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;
	uniform mat3 normal_matrix;

	uniform usamplerBuffer vertex_ind;
	uniform usamplerBuffer tri_ind;

	uniform samplerBuffer position_vertex;
	uniform samplerBuffer normal_vertex;
	uniform samplerBuffer scalar_face;

	out vec3 pos;
	out vec3 normal;
	flat out vec3 color;

	//_insert_colormap_function_here

	void main()
	{
		int ind_v = int(texelFetch(vertex_ind, 3*gl_InstanceID+gl_VertexID).r);
		int ind_t = int(texelFetch(tri_ind, gl_InstanceID).r);

		vec3 normal_in = texelFetch(normal_vertex, ind_v).rgb;
		vec3 position_in = texelFetch(position_vertex, ind_v).rgb;
		color = scalar2color(texelFetch(scalar_face, ind_t).r);

		normal = normal_matrix * normal_in;
		vec4 pos4 = model_view_matrix * vec4(position_in,1);
		pos = pos4.xyz;
		gl_Position = projection_matrix * pos4;
	}
)";

static const char* fragment_shader_source = R"(
	#version 330
	uniform vec3 light_pos;
	uniform vec4 ambiant_color;
	uniform vec4 spec_color;
	uniform float spec_coef;
	uniform bool double_side;

	in vec3 pos;
	in vec3 normal;
	flat in vec3 color;
	
	out vec3 frag_out;

	void main()
	{
		vec3 N = normalize(normal);
		vec3 L = normalize(light_pos-pos);
		if (!gl_FrontFacing)
		{
			if (!double_side)
				discard;
			N *= -1.0;
		}
		float lambert = max(0.0, dot(N,L));
		vec3 E = normalize(-pos);
		vec3 R = reflect(-L, N);
		float specular = pow(max(dot(R,E), 0.0), spec_coef);
		frag_out = ambiant_color.rgb + lambert*color + specular*spec_color.rgb;
	}
)";

ShaderPhongScalarPerFace::ShaderPhongScalarPerFace()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"), shader_function::ColorMap::source);
	load2_bind(v_src, fragment_shader_source);

	add_uniforms("vertex_ind", "tri_ind", "position_vertex", "normal_vertex", "scalar_face", "light_pos",
				 "ambiant_color", "spec_color", "spec_coef", "double_side", shader_function::ColorMap::name[0],
				 shader_function::ColorMap::name[1], shader_function::ColorMap::name[2],
				 shader_function::ColorMap::name[3]);

	nb_attributes_ = 3;
}

void ShaderParamPhongScalarPerFace::set_uniforms()
{
	vbos_[0]->bind_tb(12);
	vbos_[1]->bind_tb(13);
	vbos_[2]->bind_tb(14);
	shader_->set_uniforms_values(10, 11, 12, 13, 14, light_position_, ambiant_color_, specular_color_, specular_coef_,
								 double_side_, cm_.color_map_, cm_.expansion_, cm_.min_value_, cm_.max_value_);
}

} // namespace rendering

} // namespace cgogn
