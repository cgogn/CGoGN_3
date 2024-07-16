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

#define CGOGN_RENDER_SHADERS_FLAT_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_flat_texture.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatTexture* ShaderFlatTexture::instance_ = nullptr;

ShaderFlatTexture::ShaderFlatTexture()
{
    	const char* vertex_shader_source = R"(
        #version 330
        uniform mat4 projection_matrix;
        uniform mat4 model_view_matrix;

        in vec3 vertex_position;
        in vec2 vertex_tc;
        
        out vec3 position;
        out vec2 tc;
        void main()
        {
            vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
            position = position4.xyz;
            tc = vertex_tc;
            gl_Position = projection_matrix * position4;
        }
    )";

	const char* fragment_shader_source = R"(
		#version 330

		uniform sampler2D texture_img_unit;
		uniform vec3 light_position;
		uniform bool draw_param;

		in vec3 position;
		in vec2 tc;

		out vec3 frag_out;

		void main()
		{
			vec3 N = normalize(cross(dFdx(position), dFdy(position)));
			vec3 L = normalize(light_position - position);
			float lambert = 0.25 + 0.75 * max(0.0,dot(N, L));
			frag_out = (draw_param) ? vec3(tc,0) : texture(texture_img_unit,vec2(tc.x,1.0-tc.y)).rgb*lambert;
		}
	)";

    load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position", "vertex_tc");
	get_uniforms("texture_img_unit","light_position", "draw_param");
}

void ShaderParamFlatTexture::set_uniforms()
{
	shader_->set_uniforms_values(texture_->bind(0), light_position_, draw_param_);
}

} // namespace rendering

} // namespace cgogn
