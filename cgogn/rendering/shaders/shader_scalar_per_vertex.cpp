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

#include <cgogn/rendering/shaders/shader_scalar_per_vertex.h>

#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderFlatScalarPerVertex* ShaderFlatScalarPerVertex::instance_ = nullptr;
ShaderPhongScalarPerVertex* ShaderPhongScalarPerVertex::instance_ = nullptr;

static const char* vertex_shader_source = R"(
	#version 150

	//_insert_defines_here

	uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;

	in vec3 vertex_pos;
	in float vertex_scalar;

	out vec3 pos;
	out vec3 color;
	out float scalar;

	#if WITH_NORMAL==1
	uniform mat3 normal_matrix;
	in vec3 vertex_normal;
	out vec3 normal;
	#endif

	//_insert_colormap_function_here

	void main()
	{
		scalar = transform_value(vertex_scalar);
		color = value2color(scalar);
		#if WITH_NORMAL==1
		normal = normal_matrix * vertex_normal;
		#endif
		vec4 pos4 = model_view_matrix * vec4(vertex_pos, 1.0);
		pos = pos4.xyz;
		gl_Position = projection_matrix * pos4;
	}
)";

static const char* fragment_shader_source = R"(
	#version 150

	//_insert_defines_here

	uniform vec3 light_pos;
	uniform vec4 ambiant_color;
	uniform bool double_side;
	uniform bool show_iso_lines;
	uniform int nb_iso_levels;

	#if WITH_NORMAL==1
	uniform vec4 spec_color;
	uniform float spec_coef;
	in vec3 normal;
	#endif

	in vec3 pos;
	in vec3 color;
	in float scalar;

	out vec3 frag_out;

	void main()
	{
		#if WITH_NORMAL==1
		vec3 N = normal;
		#else
		vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));
		#endif
		vec3 L = normalize(light_pos-pos);
		if (!gl_FrontFacing)
		{
			if (!double_side)
				discard;
			N *= -1.0;
		}
		float lambert = max(0.0, dot(N,L));
		#if WITH_NORMAL==1
		vec3 E = normalize(-pos);
		vec3 R = reflect(-L, N);
		float specular = pow(max(dot(R,E), 0.0), spec_coef);
		if (show_iso_lines)
		{
			float s = scalar * float(nb_iso_levels);
			if (s - floor(s) < 0.05)
				frag_out = vec3(0.0);
			else
				frag_out = ambiant_color.rgb + lambert*color + specular*spec_color.rgb;
		}
		else
			frag_out = ambiant_color.rgb + lambert*color + specular*spec_color.rgb;
		#else
		if (show_iso_lines)
		{
			float s = scalar * float(nb_iso_levels);
			if (s - floor(s) < 0.05)
				frag_out = vec3(0.0);
			else
				frag_out = ambiant_color.rgb + lambert*color;
		}
		else
			frag_out = ambiant_color.rgb + lambert*color;
		#endif
	}
)";

ShaderFlatScalarPerVertex::ShaderFlatScalarPerVertex()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_defines_here"), "#define WITH_NORMAL 0\n");
	v_src.insert(v_src.find("//_insert_colormap_function_here"), shader_function::color_maps_shader_source());
	std::string f_src(fragment_shader_source);
	f_src.insert(f_src.find("//_insert_defines_here"), "#define WITH_NORMAL 0\n");

	load2_bind(v_src, f_src, "vertex_pos", "vertex_scalar");
	add_uniforms("color_map", "expansion", "min_value", "max_value", "show_iso_lines", "nb_iso_levels",
				 "light_position");
}

ShaderPhongScalarPerVertex::ShaderPhongScalarPerVertex()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_defines_here"), "#define WITH_NORMAL 1\n");
	v_src.insert(v_src.find("//_insert_colormap_function_here"), shader_function::color_maps_shader_source());
	std::string f_src(fragment_shader_source);
	f_src.insert(f_src.find("//_insert_defines_here"), "#define WITH_NORMAL 1\n");

	load2_bind(v_src, f_src, "vertex_pos", "vertex_normal", "vertex_scalar");
	add_uniforms("color_map", "expansion", "min_value", "max_value", "show_iso_lines", "nb_iso_levels",
				 "light_position");
}

} // namespace rendering

} // namespace cgogn
