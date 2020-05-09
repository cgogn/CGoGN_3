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

#include <cgogn/rendering/shaders/frame_manip_drawer.h>

namespace cgogn
{

namespace rendering
{

ShaderRings* ShaderRings::instance_ = nullptr;

ShaderRings::ShaderRings()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform int selected;

		out vec2 lc;
		flat out float scale;
		flat out vec3 color;

		void main()
		{
			lc = vec2(gl_VertexID % 2, gl_VertexID / 2) * 3.0 - 1.5;
			vec4 P = vec4(0, 0, 0, 1);
			P[(gl_InstanceID + 1) % 3] = lc.x;
			P[(gl_InstanceID + 2) % 3] = lc.y;

			vec3 No = vec3(0, 0, 0);
			No[gl_InstanceID] = 1.0;
			mat3 normal_mat = mat3(model_view_matrix[0].xyz, model_view_matrix[1].xyz, model_view_matrix[2].xyz);
			float lambert = 0.2 + 0.5 * abs(normalize(normal_mat * No).z);
			color = mix(lambert, 1.0, float(selected == gl_InstanceID)) * No;
			gl_Position = projection_matrix * model_view_matrix * P;
			scale = 1.0 + 0.5 * float(selected == gl_InstanceID);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform int selected;
		
		in vec2 lc;
		flat in float scale;
		flat in vec3 color;
		
		out vec3 frag_out;
		
		void main()
		{
			float d = dot(lc, lc);
			if (abs(d - 1.0) > 0.18 * scale)
				discard;
			frag_out = color;
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source);
	add_uniforms("selected");
}

void ShaderParamRings::set_uniforms()
{
	shader_->set_uniform_value(0, selected_);
}

ShaderAxis* ShaderAxis::instance_ = nullptr;

ShaderAxis::ShaderAxis()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform int selected;
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec2 lineWidths;

		out vec3 colorVert;
		out vec2 line_widths;

		void main()
		{
			vec3 P = vec3(0,0,0);
			P[gl_InstanceID] = gl_VertexID==0 ? 0.0 : 1.0;
			colorVert = selected == gl_InstanceID ? vec3(1,1,0) : P;
			gl_Position = projection_matrix * model_view_matrix * vec4(P*0.75,1);
			line_widths = selected == gl_InstanceID ? 2.0*lineWidths : lineWidths;
		}
	)";

	const char* geometry_shader_source = R"(
		#version 330
		layout (lines) in;
		layout (triangle_strip, max_vertices=15) out;

		uniform int selected;

		in vec2 line_widths[];
		in vec3 colorVert[];
		
		out vec3 color;

		void main()
		{
			vec4 A = gl_in[0].gl_Position;
			vec4 B = gl_in[1].gl_Position;
			A = A / A.w;
			B = B / B.w;
			vec4 C = 0.75*A+ 0.25*B;

			vec2 U2 = normalize(vec2(line_widths[0][1],line_widths[0][0])*(B.xy - A.xy));
			vec2 LWCorr = line_widths[0] * max(abs(U2.x),abs(U2.y));

			vec4 U = 4.0*vec4(LWCorr*U2,0,0);
			vec4 V = vec4(LWCorr*vec2(U2[1], -U2[0]), 0,0);

			color = vec3(0.8,0.8,0.8);
			gl_Position = (A-V);
			EmitVertex();
			gl_Position = (C-V);
			EmitVertex();
			gl_Position = (A);
			EmitVertex();
			gl_Position = (C);
			EmitVertex();
			gl_Position = (A+V);
			EmitVertex();
			gl_Position = (C+V);
			EmitVertex();
			EndPrimitive();

			color = 0.5*colorVert[1];
			gl_Position = (C-V);
			EmitVertex();
			color = 0.5*colorVert[1];
			gl_Position = (B-V);
			EmitVertex();
			color = colorVert[1];
			gl_Position = (C);
			EmitVertex();
			color = colorVert[1];
			gl_Position = (B);
			EmitVertex();
			color = 0.5*colorVert[1];
			gl_Position = (C+V);
			EmitVertex();
			color = 0.5*colorVert[1];
			gl_Position = (B+V);
			EmitVertex();
			color = 0.5*colorVert[1];
			EndPrimitive();

			V *= 4.0;
			color = colorVert[1];
			gl_Position = (B-V);
			EmitVertex();
			gl_Position = (B+U);
			EmitVertex();
			gl_Position = (B+V);
			EmitVertex();
			EndPrimitive();
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		in vec3 color;
		
		out vec3 frag_out;
		
		void main()
		{
			frag_out = color;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source);
	add_uniforms("lineWidths", "selected");
}

void ShaderParamAxis::set_uniforms()
{
	int viewport[4];
	float32 width = 8.0f;
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLVec2 wd(width / float32(viewport[2]), width / float32(viewport[3]));
	shader_->set_uniforms_values(wd, selected_);
}

ShaderXYGrid* ShaderXYGrid::instance_ = nullptr;

ShaderXYGrid::ShaderXYGrid()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform int nb; // nb -1  done outside
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform float sc;

		void main()
		{
			int dir = gl_InstanceID % 2;
			int line_id = gl_InstanceID / 2;
			vec2 P;
			P[dir] = float(gl_VertexID);
			P[(dir+1)%2] = float(line_id) / float(nb);
			P = P * 2.0 - 1.0;
			gl_Position = projection_matrix * model_view_matrix * vec4(P * sc, 0.001, 1);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform vec4 color;

		out vec4 frag_out;

		void main()
		{
			frag_out = color;
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source);
	add_uniforms("sc", "nb", "color");
}

void ShaderParamXYGrid::set_uniforms()
{
	shader_->set_uniforms_values(sc_, nb_ - 1, color_);
}

FrameManipDrawer* FrameManipDrawer::instance_ = nullptr;

FrameManipDrawer::FrameManipDrawer()
{
	param_rings_ = ShaderRings::generate_param();
	param_axis_ = ShaderAxis::generate_param();
	param_grid_ = ShaderXYGrid::generate_param();
}

FrameManipDrawer::~FrameManipDrawer()
{
}

void FrameManipDrawer::draw_transla(const GLMat4& projection, const GLMat4& view, const GLMat4& frame)
{
	param_axis_->bind(projection, view * frame);
	glDrawArraysInstanced(GL_LINES, 0, 2, 3);
	param_axis_->release();
}

void FrameManipDrawer::draw_rota(const GLMat4& projection, const GLMat4& view, const GLMat4& frame)
{
	param_rings_->bind(projection, view * frame);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 3);
	param_rings_->release();
}

void FrameManipDrawer::draw_grid(const GLMat4& projection, const GLMat4& view, const GLMat4& frame, float32 scale)
{
	param_grid_->sc_ = scale;
	param_grid_->nb_ = GLint(5 * scale);
	param_grid_->bind(projection, view * frame);
	glDrawArraysInstanced(GL_LINES, 0, 2, 2 * (param_grid_->nb_));
	param_grid_->release();
}

} // namespace rendering

} // namespace cgogn
