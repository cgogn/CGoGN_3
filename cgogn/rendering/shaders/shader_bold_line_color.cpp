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

#define CGOGN_RENDER_SHADERS_BOLD_LINE_CPP_

#include <cgogn/rendering/shaders/shader_bold_line_color.h>

namespace cgogn
{

namespace rendering
{

ShaderBoldLineColor* ShaderBoldLineColor::instance_ = nullptr;

ShaderBoldLineColor::ShaderBoldLineColor()
{
	const char* vertex_shader_source = R"(
		#version 150
		in vec3 vertex_position;
		in vec3 vertex_color;

		out vec3 color_v;
		
		void main()
		{
			color_v = vertex_color;
			gl_Position = vec4(vertex_position, 1.0);
		}
	)";

	const char* geometry_shader_source = R"(
		#version 150
		layout (lines) in;
		layout (triangle_strip, max_vertices = 6) out;
		
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec2 line_width;
		
		in vec3 color_v[];
		out vec4 color_f;
		out vec4 posi_clip;
		
		void main()
		{
			vec4 A = model_view_matrix * gl_in[0].gl_Position;
			vec4 B = model_view_matrix * gl_in[1].gl_Position;
			float nearZ = 1.0;
			if (projection_matrix[2][2] !=  1.0)
				nearZ = - projection_matrix[3][2] / (projection_matrix[2][2] - 1.0); 
			if (A.z < nearZ || B.z < nearZ)
			{
				if (A.z >= nearZ)
					A = B + (A - B) * (nearZ - B.z) / (A.z - B.z);
				if (B.z >= nearZ)
					B = A + (B - A) * (nearZ - A.z) / (B.z - A.z);
				A = projection_matrix * A;
				B = projection_matrix * B;
				A = A / A.w;
				B = B / B.w;
				vec2 U2 = normalize(vec2(line_width[1], line_width[0]) * (B.xy - A.xy));
				vec2 LWCorr = line_width * max(abs(U2.x), abs(U2.y));
				vec3 U = vec3(LWCorr * U2, 0.0);
				vec3 V = vec3(LWCorr*vec2(U2[1], -U2[0]), 0.0);
				color_f = vec4(color_v[0], 0.0);
				posi_clip = gl_in[0].gl_Position;
				gl_Position = vec4(A.xyz - V, 1.0);
				EmitVertex();
				color_f = vec4(color_v[1], 0.0);
				posi_clip = gl_in[1].gl_Position;
				gl_Position = vec4(B.xyz - V, 1.0);
				EmitVertex();
				color_f = vec4(color_v[0], 1.0);
				posi_clip = gl_in[0].gl_Position;
				gl_Position = vec4(A.xyz - U, 1.0);
				EmitVertex();
				color_f = vec4(color_v[1], 1.0);
				posi_clip = gl_in[1].gl_Position;
				gl_Position = vec4(B.xyz + U, 1.0);
				EmitVertex();
				color_f = vec4(color_v[0], 0.0);
				posi_clip = gl_in[0].gl_Position;
				gl_Position = vec4(A.xyz + V, 1.0);
				EmitVertex();
				color_f = vec4(color_v[1], 0.0);
				posi_clip = gl_in[1].gl_Position;
				gl_Position = vec4(B.xyz + V, 1.0);
				EmitVertex();
				EndPrimitive();
			}
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;

		in vec4 color_f;
		in vec4 posi_clip;
		out vec4 frag_out;
		
		void main()
		{
			float d = dot(plane_clip, posi_clip);
			float d2 = dot(plane_clip2, posi_clip);
			if (d > 0.0 || d2 > 0.0)
				discard;
			frag_out = color_f;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position", "vertex_color");
	get_uniforms("line_width", "plane_clip", "plane_clip2");
}

void ShaderParamBoldLineColor::set_uniforms()
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLVec2 width(width_ / float32(viewport[2]), width_ / float32(viewport[3]));
	shader_->set_uniforms_values(width, plane_clip_, plane_clip2_);
}

} // namespace rendering

} // namespace cgogn
