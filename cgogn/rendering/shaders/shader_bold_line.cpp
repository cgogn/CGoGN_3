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
		in vec3 clipping_position;
		out vec3 clip_pos_v;
		void main()
		{
			gl_Position =  vec4(vertex_position, 1.0);
			clip_pos_v = clipping_position;
		}
	)";

	const char* geometry_shader_source = R"(
		#version 330
		layout (lines) in;
		layout (triangle_strip, max_vertices = 6) out;
		
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec2 line_width;

		in vec3 clip_pos_v[];
		out float Nz;
		out vec4 posi_clip;

		void main()
		{
			vec4 A = model_view_matrix * gl_in[0].gl_Position;
			vec4 B = model_view_matrix * gl_in[1].gl_Position;
			float nearZ = 1.0;
			if (projection_matrix[2][2] !=  1.0)
				nearZ = - projection_matrix[3][2] / (projection_matrix[2][2] - 1.0);
			// float farZ = - projection_matrix[3][2] / (projection_matrix[2][2] + 1.0);

			if (A.z < nearZ || B.z < nearZ)
			{
				if (A.z >= nearZ)
					A = B + (A - B) * (nearZ - B.z) / (A.z - B.z);
				if (B.z >= nearZ)
					B = A + (B - A) * (nearZ - A.z) / (B.z - A.z);

				// vec3 AB = B.xyz / B.w - A.xyz / A.w;
				// vec3 Nl = normalize(cross(AB, vec3(0.0, 0.0, 1.0)));
				// vec3 Nm = vec3(0.0, 0.0, 1.0);

				A = projection_matrix * A;
				B = projection_matrix * B;
				A = A / A.w;
				B = B / B.w;

				vec2 U2 = normalize(vec2(line_width[1], line_width[0]) * (B.xy - A.xy));
				vec2 LWCorr = line_width * max(abs(U2.x), abs(U2.y));

				vec4 U = vec4(0.5 * LWCorr * U2, 0.0, 0.0);
				vec4 V = vec4(LWCorr * vec2(U2[1], -U2[0]), 0.0, 0.0);
				posi_clip = vec4(clip_pos_v[0],1);
				Nz = 0;
				gl_Position = A - V;
				EmitVertex();
				posi_clip = vec4(clip_pos_v[1],1);
				Nz = 0;
				gl_Position = B - V;
				EmitVertex();
				posi_clip = vec4(clip_pos_v[0],1);
				Nz = 1;
				gl_Position = A - U;
				EmitVertex();
				posi_clip = vec4(clip_pos_v[1],1);
				Nz = 1;
				gl_Position = B + U;
				EmitVertex();
				posi_clip = vec4(clip_pos_v[0],1);
				Nz = 0;
				gl_Position = A + V;
				EmitVertex();
				posi_clip = vec4(clip_pos_v[1],1);
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

		in float Nz;
		in vec4 posi_clip;
		
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

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position",
			   "clipping_position");
	get_uniforms("line_color", "line_width", "lighted", "plane_clip", "plane_clip2");
}

void ShaderParamBoldLine::set_uniforms()
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLVec2 width(width_ / float32(viewport[2]), width_ / float32(viewport[3]));
	shader_->set_uniforms_values(color_, width, lighted_, plane_clip_, plane_clip2_);
}

ShaderBoldLineColor* ShaderBoldLineColor::instance_ = nullptr;

ShaderBoldLineColor::ShaderBoldLineColor()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform usamplerBuffer vertex_ind;
		uniform usamplerBuffer edge_ind;
		uniform samplerBuffer vertex_position;
		uniform samplerBuffer edge_color;
		
		out vec3 color_e;

		void main()
		{
			int ind_v = int(texelFetch(vertex_ind, 2 * gl_InstanceID + gl_VertexID).r);
			vec3 position_in = texelFetch(vertex_position, ind_v).rgb;

			int ind_e = int(texelFetch(edge_ind, int(gl_InstanceID)).r);
			color_e = texelFetch(edge_color, ind_e).rgb;

			gl_Position = vec4(position_in, 1.0);
		}
	)";

	const char* geometry_shader_source = R"(
		#version 330
		layout (lines) in;
		layout (triangle_strip, max_vertices = 6) out;

		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec2 line_width;

		in vec3 color_e[];
		
		out vec3 color;
		out float Nz;
		out vec4 posi_clip;

		void main()
		{
			vec4 A = model_view_matrix * gl_in[0].gl_Position;
			vec4 B = model_view_matrix * gl_in[1].gl_Position;
			float nearZ = 1.0;
			if (projection_matrix[2][2] !=  1.0)
				nearZ = - projection_matrix[3][2] / (projection_matrix[2][2] - 1.0);
			// float farZ = - projection_matrix[3][2] / (projection_matrix[2][2] + 1.0);

			if (A.z < nearZ || B.z < nearZ)
			{
				if (A.z >= nearZ)
					A = B + (A - B) * (nearZ - B.z) / (A.z - B.z);
				if (B.z >= nearZ)
					B = A + (B - A) * (nearZ - A.z) / (B.z - A.z);

				// vec3 AB = B.xyz / B.w - A.xyz / A.w;
				// vec3 Nl = normalize(cross(AB, vec3(0.0, 0.0, 1.0)));
				// vec3 Nm = vec3(0.0, 0.0, 1.0);

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
				color = color_e[0];
				gl_Position = A - V;
				EmitVertex();
				posi_clip = gl_in[1].gl_Position;
				Nz = 0;
				color = color_e[0];
				gl_Position = B - V;
				EmitVertex();
				posi_clip = gl_in[0].gl_Position;
				Nz = 1;
				color = color_e[0];
				gl_Position = A - U;
				EmitVertex();
				posi_clip = gl_in[1].gl_Position;
				Nz = 1;
				color = color_e[0];
				gl_Position = B + U;
				EmitVertex();
				posi_clip = gl_in[0].gl_Position;
				Nz = 0;
				color = color_e[0];
				gl_Position = A + V;
				EmitVertex();
				posi_clip = gl_in[1].gl_Position;
				Nz = 0;
				color = color_e[0];
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
		uniform float lighted;

		in vec3 color;
		in float Nz;
		in vec4 posi_clip;

		out vec3 frag_out;
		
		void main()
		{
			float d = dot(plane_clip, posi_clip);
			float d2 = dot(plane_clip2, posi_clip);
			if (d > 0.0 || d2 > 0.0)
				discard;

			float lambert = max(1.0 - lighted, Nz); // Nz = dot(N,0,0,1)
			frag_out = color * lambert;
		}
	)";

	load3(vertex_shader_source, fragment_shader_source, geometry_shader_source);
	get_uniforms("vertex_ind", "edge_ind", "vertex_position", "edge_color", "line_width", "lighted", "plane_clip",
				 "plane_clip2");

	nb_attributes_ = 2;
}

void ShaderParamBoldLineColor::set_uniforms()
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLVec2 width(width_ / float32(viewport[2]), width_ / float32(viewport[3]));
	shader_->set_uniforms_values(10, 11, 12, 13, width, lighted_, plane_clip_, plane_clip2_);
}

void ShaderParamBoldLineColor::bind_texture_buffers()
{
	vbos_[VERTEX_POSITION]->bind_texture_buffer(12);
	vbos_[EDGE_COLOR]->bind_texture_buffer(13);
}

void ShaderParamBoldLineColor::release_texture_buffers()
{
	vbos_[VERTEX_POSITION]->release_texture_buffer(12);
	vbos_[EDGE_COLOR]->release_texture_buffer(13);
}

} // namespace rendering

} // namespace cgogn
