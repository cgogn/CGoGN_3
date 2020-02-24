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

#include <cgogn/rendering/shaders/shader_phong_color_per_vertex.h>

namespace cgogn
{

namespace rendering
{

ShaderPhongColorPerVertex* ShaderPhongColorPerVertex::instance_ = nullptr;

ShaderPhongColorPerVertex::ShaderPhongColorPerVertex()
{
	static const char* vertex_shader_source = R"(
		#version 150
		in vec3 vertex_pos;
		in vec3 vertex_normal;
		in vec3 vertex_color;

		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform mat3 normal_matrix;

		out vec3 P;
		out vec3 Normal;
		out vec3 Color;

		void main ()
		{
			Normal = normal_matrix * vertex_normal;
			Color = vertex_color;
			vec4 P4 = model_view_matrix * vec4 (vertex_pos, 1.0);
			P = P4.xyz;
			gl_Position = projection_matrix * P4;
		};
		)";

	static const char* fragment_shader_source =
		R"(	#version 150
		in vec3 P;
		in vec3 Normal;
		in vec3 Color;

		uniform vec4 spec_color;
		uniform vec4 ambiant_color;
		uniform float spec_coef;
		uniform vec3 light_position;
		uniform bool double_side;

		out vec3 frag_color;

		void main()
		{
			vec3 N = normalize(Normal);
			vec3 L = normalize(light_position-P);
			vec3 finalColor = ambiant_color.rgb;
			if (gl_FrontFacing==false) // do not use ! because of bug on old intel under OS/X
			{
				if (!double_side)
					discard;
				N *= -1.0;
			}
			float lambertTerm = clamp(dot(N,L),0.0,1.0);
			finalColor += Color*lambertTerm ;
			vec3 E = normalize(-P);
			vec3 R = reflect(-L, N);
			float specular = pow( max(dot(R, E), 0.0), spec_coef );
			finalColor += spec_color.rgb * specular;
			frag_color = finalColor;
		};
		)";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos", "vertex_normal", "vertex_color");
	add_uniforms("light_position", "ambiant_color", "spec_color", "spec_coef", "double_side");
}

void ShaderParamPhongColorPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(light_position_, ambiant_color_, specular_color_, specular_coef_, double_side_);
}

} // namespace rendering

} // namespace cgogn
