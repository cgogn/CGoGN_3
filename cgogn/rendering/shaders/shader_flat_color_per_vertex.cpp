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

ShaderFlatColorPerVertex* ShaderFlatColorPerVertex::instance_ = nullptr;

ShaderFlatColorPerVertex::ShaderFlatColorPerVertex()
{
	const char* vertex_shader_source = R"(
        #version 150
        uniform mat4 projection_matrix;
        uniform mat4 model_view_matrix;

        in vec3 vertex_position;
        in vec3 vertex_color;
        
        out vec3 position;
        out vec3 color;
        
        void main()
        {
            vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
            position = position4.xyz;
            color = vertex_color;
            gl_Position = projection_matrix * position4;
        }
    )";

	const char* fragment_shader_source = R"(
        #version 150
        uniform vec4 ambiant_color;
        uniform vec3 light_position;
        uniform bool double_side;
        
        in vec3 position;
        in vec3 color;

        out vec3 frag_out;
        
        void main()
        {
            vec3 N = normalize(cross(dFdx(position), dFdy(position)));
            vec3 L = normalize(light_position - position);
            float lambert = dot(N, L);
            if (gl_FrontFacing==false && !double_side)
                discard;
            else
                frag_out = ambiant_color.rgb + lambert * color.rgb;
        }
    )";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position", "vertex_color");
	get_uniforms("ambiant_color", "light_position", "double_side");
}

void ShaderParamFlatColorPerVertex::set_uniforms()
{
	shader_->set_uniforms_values(ambiant_color_, light_position_, double_side_);
}

} // namespace rendering

} // namespace cgogn
