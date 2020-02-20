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

#include <cgogn/rendering/shaders/shader_flat_color_per_vertex.h>

namespace cgogn
{

namespace rendering
{

const char *vertex_shader_source =
    R"(#version 150
in vec3 vertex_pos;
in vec3 vertex_color;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
out vec3 pos;
out vec3 color;
void main()
{
    vec4 pos4 = model_view_matrix * vec4(vertex_pos,1.0);
    pos = pos4.xyz;
    gl_Position = projection_matrix * pos4;
};
)";

const char *fragment_shader_source =
    R"(#version 150
out vec3 fragColor;
uniform vec4 ambiant_color;
uniform vec3 light_position;
uniform bool double_side;
in vec3 pos;
in vec3 color;
void main()
{
    vec3 N = normalize(cross(dFdx(pos),dFdy(pos)));
    vec3 L = normalize(light_position-pos);
    float lambert = dot(N,L);
    if (!gl_FrontFacing && !double_side)
         discard;
    else
        fragColor = ambiant_color.rgb + lambert*color.rgb;
};
)";

ShaderFlatColorPerVertex *ShaderFlatColorPerVertex::instance_ = nullptr;

ShaderFlatColorPerVertex::ShaderFlatColorPerVertex()
{
    load2_bind(vertex_shader_source, fragment_shader_source, "vertex_pos");

    add_uniforms("ambiant_color", "light_position", "double_side");
}

} // namespace rendering

} // namespace cgogn
