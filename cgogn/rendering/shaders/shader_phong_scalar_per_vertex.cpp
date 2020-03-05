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

#include <cgogn/rendering/shaders/shader_phong_scalar_per_vertex.h>

namespace cgogn
{

namespace rendering
{

static const char* vertex_shader_source = R"(#version 150
in vec3 vertex_pos;
in vec3 vertex_normal;
in float vertex_scalar;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;
uniform vec3 light_position;
out vec3 EyeVector;
out vec3 Normal;
out vec3 LightDir;
out vec3 color;

//_insert_colormap_funcion_here

void main()
{
	Normal = normal_matrix * vertex_normal;
	vec3 Position = vec3(model_view_matrix * vec4(vertex_pos, 1.0));
	LightDir = light_position - Position;
	EyeVector = -Position;
	color = scalar2color(vertex_scalar);
	gl_Position = projection_matrix * model_view_matrix * vec4(vertex_pos, 1.0);
}
)";

static const char* fragment_shader_source = R"(#version 150
in vec3 EyeVector;
in vec3 Normal;
in vec3 LightDir;
in vec3 color;
uniform vec4 ambiant_color;
uniform vec4 spec_color;
uniform float spec_coef;
uniform bool double_side;
out vec3 frag_color;
void main()
{
	vec3 N = normalize(Normal);
	vec3 L = normalize(LightDir);
	vec3 finalColor = ambiant_color.rgb;
	if (gl_FrontFacing == false) // do not use ! because of bug on old intel under OS/X
	{
		if (!double_side)
			discard;
		N *= -1.0;
	}
	float lambertTerm = max(dot(N, L), 0.0);
	finalColor += color * lambertTerm;
	vec3 E = normalize(EyeVector);
	vec3 R = reflect(-L, N);
	float specular = pow(max(dot(R, E), 0.0), spec_coef);
	finalColor += spec_color.rgb * specular;
	frag_color = finalColor;
}
)";

ShaderPhongScalarPerVertex* ShaderPhongScalarPerVertex::instance_ = nullptr;

ShaderPhongScalarPerVertex::ShaderPhongScalarPerVertex()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_colormap_funcion_here"), shader_funcion::ColorMap::source);
	load2_bind(v_src, fragment_shader_source, "vertex_pos", "vertex_normal", "vertex_scalar");
	add_uniforms("light_position", "ambiant_color", "spec_color", "spec_coef", "double_side",
				 shader_funcion::ColorMap::name[0], shader_funcion::ColorMap::name[1],
				 shader_funcion::ColorMap::name[2], shader_funcion::ColorMap::name[3]);
}

void ShaderParamPhongScalarPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(light_position_, ambiant_color_, specular_color_, specular_coef_, double_side_,
								 cm_.color_map_, cm_.expansion_, cm_.min_value_, cm_.max_value_);
}

} // namespace rendering

} // namespace cgogn
