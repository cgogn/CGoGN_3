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
static const char* vertex_shader_source = R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform usamplerBuffer tri_ind;
uniform usamplerBuffer face_emb;
uniform samplerBuffer pos_vertex;
uniform samplerBuffer color_face;
out vec3 A;
flat out vec3 N;
flat out vec3 color;
void main()
{
	int tri = int(texelFetch(tri_ind, gl_InstanceID).r);
	int i_c = int(texelFetch(face_emb, gl_InstanceID).r);
	color = texelFetch(color_face, i_c).rgb;
	int vid = gl_VertexID;
	int ind_a = int(texelFetch(tri_ind, 3*int(gl_InstanceID)+vid).r);
	A = (model_view_matrix * vec4(texelFetch(pos_vertex, ind_a).rgb,1.0)).xyz;
	vid  = (vid+1)%3;
	int ind_b = int(texelFetch(tri_ind, 3*int(gl_InstanceID)+vid).r);
	vec3 B = (model_view_matrix * vec4(texelFetch(pos_vertex, ind_b).rgb,1.0)).xyz;
	vid  = (vid+1)%3;
	int ind_c = int(texelFetch(tri_ind, 3*int(gl_InstanceID)+vid).r);
	vec3 C = (model_view_matrix * vec4(texelFetch(pos_vertex, ind_c).rgb,1.0)).xyz;
	N = normalize(cross(B-A,C-A));
	gl_Position = projection_matrix*vec4(A,1);
}
)";

static const char* fragment_shader_source = R"(#version 330
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
		fragColor = ambiant_color.rgb+lambert*color;
	else
		discard;
}
)";

ShaderFlatColorPerFace* ShaderFlatColorPerFace::instance_ = nullptr;

ShaderFlatColorPerFace::ShaderFlatColorPerFace()
{

	load2_bind(vertex_shader_source, fragment_shader_source, "");

	add_uniforms("tri_ind", "face_emb", "pos_vertex", "color_face", "ambiant_color", "light_position", "double_side");
}

void ShaderParamFlatColorPerFace::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, vbo_pos_->bind_tb(12), vbo_color_->bind_tb(13), ambiant_color_,
								 light_position_, double_side_);
}

void ShaderParamFlatColorPerFace::set_vbos(const std::vector<VBO*>& vbos)
{
	vbo_pos_ = vbos[0];
	vbo_color_ = vbos[1];
	vao_initialized_ = vbos[0] != nullptr && vbos[1] != nullptr;
}

} // namespace rendering

} // namespace cgogn
