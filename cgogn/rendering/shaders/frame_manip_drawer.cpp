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

#include <cgogn/rendering/shaders/frame_manip_drawer.h>

namespace cgogn
{

namespace rendering
{

ShaderRings* ShaderRings::instance_ = nullptr;

static const char* vertex_shader_source =
	R"(
#version 330
uniform int selected;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;
out vec2 lc;
out vec3 color;
void main()
{
	lc = vec2(gl_VertexID%2,gl_VertexID/2)*2.0-1.0;
	vec4 P = vec4(0,0,0,1);
	P[(gl_InstanceID+1)%3] = lc.x;
	P[(gl_InstanceID+2)%3] = lc.y;

	vec3 No = vec3(0,0,0);
	No[gl_InstanceID] = 1.0;
	float lambert = 0.6;//abs((normal_matrix*No).z);
	if (selected == gl_InstanceID)
		color = (0.4+lambert) * vec3(1,1,0);
	else
		color = (0.4+lambert) * No;

	gl_Position = projection_matrix * model_view_matrix * P;
}
)";

static const char* fragment_shader_source =
	R"(
#version 330
in vec2 lc;
in vec3 color;
out vec3 frag_out;
void main()
{
	float d = dot(lc,lc);
	if (d>1 || d<0.64)
		discard;
	frag_out = color;
}
)";

ShaderRings::ShaderRings()
{
	load2_bind(vertex_shader_source, fragment_shader_source);
	add_uniforms("selected");
}

void ShaderParamRings::set_uniforms()
{
	shader_->set_uniform_value(0, selected_);
}

FrameManipDrawer* FrameManipDrawer::instance_ = nullptr;

FrameManipDrawer::FrameManipDrawer()
{
	param_rings_ = ShaderRings::generate_param();
}

FrameManipDrawer::~FrameManipDrawer()
{
}

void FrameManipDrawer::draw(const GLMat4& projection, const GLMat4& view, const GLMat4& frame)
{
	param_rings_->bind(projection, view * frame);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 3);
	param_rings_->release();
}

} // namespace rendering

} // namespace cgogn
