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

#include <cgogn/rendering/shaders/shader_color_per_vertex.h>

#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderFlatColorPerVertex* ShaderFlatColorPerVertex::instance_ = nullptr;
ShaderPhongColorPerVertex* ShaderPhongColorPerVertex::instance_ = nullptr;

static const char* vertex_shader_source = R"(
	#version 150

	//_insert_defines_here

	uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;

	in vec3 vertex_pos;
	in vec3 vertex_color;

	out vec3 pos;
	out vec3 color;

	#if WITH_NORMAL==1
	uniform mat3 normal_matrix;
	in vec3 vertex_normal;
	out vec3 normal;
	#endif

	void main()
	{
		vec4 pos4 = model_view_matrix * vec4(vertex_pos, 1);
		pos = pos4.xyz;
		color = vertex_color;
		#if WITH_NORMAL==1
		normal = normal_matrix * vertex_normal;
		#endif
		gl_Position = projection_matrix * pos4;
	}
)";

static const char* fragment_shader_source = R"(
	#version 150

	//_insert_defines_here

	uniform vec3 light_pos;
	uniform vec4 ambiant_color;
	uniform bool double_side;

	#if WITH_NORMAL==1
	uniform vec4 spec_color;
	uniform float spec_coef;
	in vec3 normal;
	#endif

	in vec3 pos;
	in vec3 color;

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
		frag_out = ambiant_color.rgb + lambert*color + specular*spec_color.rgb;
		#else
		frag_out = ambiant_color.rgb + lambert*color;
		#endif
	}
)";

ShaderFlatColorPerVertex::ShaderFlatColorPerVertex()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_defines_here"), "#define WITH_NORMAL 0\n");
	std::string f_src(fragment_shader_source);
	f_src.insert(f_src.find("//_insert_defines_here"), "#define WITH_NORMAL 0\n");

	load2_bind(v_src, f_src, "vertex_pos", "vertex_color");
	add_uniforms("light_pos", "ambiant_color", "double_side");
}

ShaderPhongColorPerVertex::ShaderPhongColorPerVertex()
{
	std::string v_src(vertex_shader_source);
	v_src.insert(v_src.find("//_insert_defines_here"), "#define WITH_NORMAL 1\n");
	std::string f_src(fragment_shader_source);
	f_src.insert(f_src.find("//_insert_defines_here"), "#define WITH_NORMAL 1\n");

	load2_bind(v_src, f_src, "vertex_pos", "vertex_normal", "vertex_color");
	add_uniforms("light_pos", "ambiant_color", "spec_color", "spec_coef", "double_side");
}

} // namespace rendering

} // namespace cgogn
