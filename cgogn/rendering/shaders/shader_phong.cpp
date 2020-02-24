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

#include <cgogn/rendering/shaders/shader_phong.h>
#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderPhong* ShaderPhong::instance_ = nullptr;

static const char* vertex_shader_source =
	R"(#version 150
	in vec3 vertex_pos;
	in vec3 vertex_normal;
	uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;
	uniform mat3 normal_matrix;
	uniform vec3 light_position;
	out vec3 EyeVector;
	out vec3 Normal;
	out vec3 LightDir;
	void main ()
	{
		Normal = normal_matrix * vertex_normal;
		vec3 Position = vec3 (model_view_matrix * vec4 (vertex_pos, 1.0));
		LightDir = light_position - Position;
		EyeVector = -Position;
		gl_Position = projection_matrix * model_view_matrix * vec4 (vertex_pos, 1.0);
	};
	)";

static const char* fragment_shader_source =
	R"(#version 150
	in vec3 EyeVector;
	in vec3 Normal;
	in vec3 LightDir;
	uniform vec4 front_color;
	uniform vec4 spec_color;
	uniform vec4 ambiant_color;
	uniform vec4 back_color;
	uniform float spec_coef;
	uniform bool double_side;
	out vec4 frag_color;
	void main()
	{
		vec3 N = normalize (Normal);
		vec3 L = normalize (LightDir);
		vec4 finalColor = ambiant_color;
		vec4 currentColor = front_color;
		if (gl_FrontFacing==false) // do not use ! because of bug on old intel under OS/X
		{
			if (!double_side)
				discard;
			N *= -1.0;
			currentColor = back_color;
		}
		float lambertTerm = clamp(dot(N,L),0.0,1.0);
		finalColor += currentColor*lambertTerm ;
		vec3 E = normalize(EyeVector);
		vec3 R = reflect(-L, N);
		float specular = pow( max(dot(R, E), 0.0), spec_coef );
		finalColor += spec_color * specular;
		frag_color=vec4(finalColor.rgb,1.0);
	};
	)";

ShaderPhong::ShaderPhong()
{
	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos", "vertex_normal");
	add_uniforms("light_position", "front_color", "back_color", "ambiant_color", "spec_color", "spec_coef",
				 "double_side");
}

void ShaderParamPhong::set_uniforms()
{
	shader_->set_uniforms_values(light_position_, front_color_, back_color_, ambiant_color_, specular_color_,
								 specular_coef_, double_side_);
}

} // namespace rendering

} // namespace cgogn
