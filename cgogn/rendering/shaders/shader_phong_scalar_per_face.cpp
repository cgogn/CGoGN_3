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

#include <cgogn/rendering/shaders/shader_phong_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_function_color_maps.h>

namespace cgogn
{

namespace rendering
{

ShaderPhongScalarPerFace* ShaderPhongScalarPerFace::instance_ = nullptr;

static const char* vertex_shader_source =
R"(
#version 330
uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;
uniform mat3 normalMatrix;

uniform usamplerBuffer tri_indices;
uniform usamplerBuffer face_emb;

uniform samplerBuffer position_vertex;
uniform samplerBuffer normal_vertex;
uniform samplerBuffer scalar_face;

out vec3 Po;
out vec3 No;
flat out vec3 Col;

//_insert_colormap_funcion_here

void main()
{
	int ind_v = int(texelFetch(tri_indices,3*gl_InstanceID+gl_VertexID).r);
	int emb = int(texelFetch(face_emb, gl_InstanceID).r);

	vec3 normal_in = texelFetch(normal_vertex, ind_v).rgb;
	vec3 position_in = texelFetch(position_vertex, ind_v).rgb;
	Col = scalar2color(texelFetch(scalar_face, emb).r);

	No = normalMatrix * normal_in;
	vec4 Po4 = viewMatrix * vec4(position_in,1);
	Po = Po4.xyz;
	gl_Position = projectionMatrix * Po4;
}
)";

static const char* fragment_shader_source =
R"(
#version 330
in vec3 Po;
in vec3 No;
flat in vec3 Col;
out vec3 frag_out;

uniform vec3 light_pos;
uniform vec4 ambiant_color;
uniform vec4 spec_color;
uniform float spec_coef;
uniform bool double_side;

void main()
{
	vec3 N = normalize(No);
	vec3 L = normalize(light_pos-Po);
	if (gl_FrontFacing==false)
	{
		if (!double_side)
			discard;
		N *= -1.0;
	}

	float lamb = max(0.0,dot(N,L));
	vec3 E = normalize(-Po);
	vec3 R = reflect(-L, N);
	float spec = pow( max(dot(R,E), 0.0), spec_coef);
	frag_out = ambiant_color.rgb + lamb*Col + spec*spec_color.rgb;
}
)";


ShaderPhongScalarPerFace::ShaderPhongScalarPerFace()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"),shader_funcion::color_map::source());
	load2_bind(v_src,fragment_shader_source);

	add_uniforms("tri_indices", "face_emb",
				 "position_vertex", "normal_vertex", "scalar_face",
				 "light_pos", "ambiant_color",
				 "spec_color", "spec_coef", "double_side",
				 shader_funcion::color_map::name[0],
				 shader_funcion::color_map::name[1],
				 shader_funcion::color_map::name[2],
				 shader_funcion::color_map::name[3]);
}

void ShaderParamPhongScalarPerFace::set_uniforms()
{
	if (vbo_pos_)
		shader_->set_uniforms_values(10,11,
						vbo_pos_->bind_tb(12),
						vbo_norm_->bind_tb(13),
						vbo_scalar_->bind_tb(14),
						light_position_,ambiant_color_,
						specular_color_,specular_coef_,double_side_,
						cm_.color_map_, cm_.expansion_, cm_.min_value_, cm_.max_value_);

}

void ShaderParamPhongScalarPerFace::set_vbos(const std::vector<VBO*>& vbos)
{
	vbo_pos_ = vbos[0];
	vbo_norm_ = vbos[1];
	vbo_scalar_ = vbos[2];
	vao_initialized_ = vbos[0]!=nullptr && vbos[1]!=nullptr && vbos[2]!=nullptr;
}

} // namespace rendering
} // namespace cgogn
