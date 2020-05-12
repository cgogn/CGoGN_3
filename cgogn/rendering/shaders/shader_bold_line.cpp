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

#include <iostream>

#include <cgogn/rendering/shaders/shader_bold_line.h>

namespace cgogn
{

namespace rendering
{

ShaderBoldLine* ShaderBoldLine::instance_ = nullptr;

ShaderBoldLine::ShaderBoldLine()
{
	const char* vertex_shader_source = R"(
		#version 330
		in vec3 vertex_position;

		void main()
		{
			gl_Position =  vec4(vertex_position, 1.0);
		}
	)";

	const char* geometry_shader_source = R"(
		#version 330
		layout (lines) in;
		layout (triangle_strip, max_vertices = 6) out;
		
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec2 line_width;

		out float Nz;
		out vec4 posi_clip;
		flat out vec3 col;

		void main()
		{
			vec4 A = model_view_matrix * gl_in[0].gl_Position;
			vec4 B = model_view_matrix * gl_in[1].gl_Position;
			float nearZ = 1.0;
			if (projection_matrix[2][2] !=  1.0)
				nearZ = - projection_matrix[3][2] / (projection_matrix[2][2] - 1.0);
			float farZ = - projection_matrix[3][2] / (projection_matrix[2][2] + 1.0);

			if (A.z < nearZ || B.z < nearZ)
			{
				if (A.z >= nearZ)
					A = B + (A - B) * (nearZ - B.z) / (A.z - B.z);
				if (B.z >= nearZ)
					B = A + (B - A) * (nearZ - A.z) / (B.z - A.z);

				vec3 AB = B.xyz / B.w - A.xyz / A.w;
				vec3 Nl = normalize(cross(AB, vec3(0.0, 0.0, 1.0)));
				vec3 Nm = vec3(0.0, 0.0, 1.0);

				A = projection_matrix * A;
				B = projection_matrix * B;
				A = A / A.w;
				B = B / B.w;

				vec2 U2 = normalize(vec2(line_width[1], line_width[0]) * (B.xy - A.xy));
				vec2 LWCorr = line_width * max(abs(U2.x), abs(U2.y));

				vec4 U = vec4(0.5 * LWCorr * U2, 0.0, 0.0);
				vec4 V = vec4(LWCorr * vec2(U2[1], -U2[0]), 0.0, 0.0);
				posi_clip = gl_in[0].gl_Position;
				Nz = 0;
				gl_Position = A - V;
				EmitVertex();
				posi_clip = gl_in[1].gl_Position;
				Nz = 0;
				gl_Position = B - V;
				EmitVertex();
				posi_clip = gl_in[0].gl_Position;
				Nz = 1;
				gl_Position = A - U;
				EmitVertex();
				posi_clip = gl_in[1].gl_Position;
				Nz = 1;
				gl_Position = B + U;
				EmitVertex();
				posi_clip = gl_in[0].gl_Position;
				Nz = 0;
				gl_Position = A + V;
				EmitVertex();
				posi_clip = gl_in[1].gl_Position;
				Nz = 0;
				gl_Position = B + V;
				EmitVertex();
				EndPrimitive();
			}
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;
		uniform vec4 line_color;
		uniform float lighted;

		in vec4 posi_clip;
		in float Nz;
		
		out vec3 frag_out;

		void main()
		{
			float d = dot(plane_clip, posi_clip);
			float d2 = dot(plane_clip2, posi_clip);
			if (d > 0.0 || d2 > 0.0)
				discard;

			float lambert = max(1.0 - lighted, Nz); // Nz = dot(N,0,0,1)
			frag_out = line_color.rgb * lambert;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position");
	get_uniforms("line_color", "line_width", "lighted", "plane_clip", "plane_clip2");
}

void ShaderParamBoldLine::set_uniforms()
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLVec2 width(width_ / float32(viewport[2]), width_ / float32(viewport[3]));
	shader_->set_uniforms_values(color_, width, lighted_, plane_clip_, plane_clip2_);
}

} // namespace rendering

} // namespace cgogn
