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

ShaderFlatScalarPerFace* ShaderFlatScalarPerFace::instance_ = nullptr;

ShaderFlatScalarPerFace::ShaderFlatScalarPerFace()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		uniform usamplerBuffer vertex_ind;
		uniform usamplerBuffer tri_ind;
		uniform samplerBuffer pos_vertex;
		uniform samplerBuffer scalar_tri;
		
		out vec3 pos;
		flat out vec3 color;

		//_insert_colormap_function_here

		void main()
		{
			int ind_v = int(texelFetch(vertex_ind, 3*gl_InstanceID+gl_VertexID).r);
			vec3 position_in = texelFetch(pos_vertex, ind_v).rgb;

			int ind_t = int(texelFetch(tri_ind, int(gl_InstanceID)).r);
			float scalar = transform_value(texelFetch(scalar_tri, ind_t).r);
			color = scalar2color(scalar);

			vec4 pos4 = model_view_matrix * vec4(position_in, 1.0);
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
				fragColor = ambiant_color.rgb + lambert*color;
			else
				discard;
		}
	)";

	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"), shader_function::ColorMap::source);
	load2_bind(v_src, fragment_shader_source, "");

	add_uniforms("tri_ind", "tri_emb", "pos_vertex", "scalar_tri", "ambiant_color", "light_position", "double_side",
				 shader_function::ColorMap::name[0], shader_function::ColorMap::name[1],
				 shader_function::ColorMap::name[2], shader_function::ColorMap::name[3]);

	nb_attributes_ = 2;
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
