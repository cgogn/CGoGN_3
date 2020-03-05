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

#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>

namespace cgogn
{

namespace rendering
{

ShaderVectorPerVertex* ShaderVectorPerVertex::instance_ = nullptr;

static const char* vertex_shader_source = R"(
#version 150
in vec3 vertex_pos;
in vec3 vertex_vector;
uniform float length;
out vec3 end_vector;
void main()
{
	end_vector = vertex_pos + length*vertex_vector;
	gl_Position =  vec4(vertex_pos,1.0);
}
)";

static const char* geometry_shader_source = R"(
#version 330
layout (points) in;
layout (triangle_strip, max_vertices=6) out;
in vec3 end_vector[];
out vec3 N;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform vec2 lineWidths;
void main()
{
	vec4 A = model_view_matrix * gl_in[0].gl_Position;
	vec4 B = model_view_matrix * vec4(end_vector[0],1);
	float nearZ = 1.0;
	if (projection_matrix[2][2] !=  1.0)
		nearZ = - projection_matrix[3][2] / (projection_matrix[2][2] - 1.0);
	if ((A.z < nearZ) || (B.z < nearZ))
	{
		if (A.z >= nearZ)
			A = B + (A-B)*(nearZ-B.z)/(A.z-B.z);
		if (B.z >= nearZ)
			B = A + (B-A)*(nearZ-A.z)/(B.z-A.z);

		vec3 AB = B.xyz/B.w - A.xyz/A.w;
		vec3 Nl = normalize(cross(AB,vec3(0,0,1)));
		vec3 Nm = normalize(cross(Nl,AB));

		A = projection_matrix*A;
		B = projection_matrix*B;
		A = A/A.w;
		B = B/B.w;
		vec2 U2 = normalize(vec2(lineWidths[1],lineWidths[0])*(B.xy - A.xy));
		vec2 LWCorr =lineWidths * max(abs(U2.x),abs(U2.y));
		vec3 U = vec3(0.5*LWCorr*U2,0.0);
		vec3 V = vec3(LWCorr*vec2(U2[1], -U2[0]), 0.0);

		N = Nl;
		gl_Position = vec4(A.xyz-V, 1.0);
		EmitVertex();
		N = Nl;
		gl_Position = vec4(B.xyz-V, 1.0);
		EmitVertex();
		N = Nm;
		gl_Position = vec4(A.xyz-U, 1.0);
		EmitVertex();
		N = Nm;
		gl_Position = vec4(B.xyz+U, 1.0);
		EmitVertex();
		N = -Nl;
		gl_Position = vec4(A.xyz+V, 1.0);
		EmitVertex();
		N = -Nl;
		gl_Position = vec4(B.xyz+V, 1.0);
		EmitVertex();
		EndPrimitive();
	}
}
)";

static const char* fragment_shader_source = R"(
#version 330
in vec3 N;
out vec3 fragColor;
uniform vec4 lineColor;
uniform float lighted;
void main()
{
	const vec3 light_dir = normalize(vec3(10,100,1000));
	float lambert = (1.0-lighted) + lighted*max(0.0,dot(N,light_dir));
	fragColor = lineColor.rgb * lambert;
}
)";

ShaderVectorPerVertex::ShaderVectorPerVertex()
{
	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_pos", "vertex_vector");
	add_uniforms("lineColor", "lineWidths", "length", "lighted");
}

void ShaderParamVectorPerVertex::set_uniforms()
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLVec2 wd(width_ / float32(viewport[2]), width_ / float32(viewport[3]));
	shader_->set_uniforms_values(color_, wd, length_, lighted_);
}

} // namespace rendering

} // namespace cgogn
