*******************************************************************************
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
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/rendering/shaders/shader_flat__per_face.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatPerFace* ShaderFlatPerFace::instance_ = nullptr;

ShaderFlatPerFacce::ShaderFlat()
{
    
	const char* vertex_shader_source =
    R"(
    #version 330
    uniform mat4 projection_matrix;
	uniform mat4 model_view_matrix;
    uniform usamplerBuffer tri_indices;
    uniform samplerBuffer position_vertex;
    uniform samplerBuffer color_face;

    out vec3 No;
    out vec3 Col;

    void main()
    {
        int t = 3*gl_InstanceID;
        int iiA = t + gl_VertexID;
        int iiB = t + (gl_VertexID+1)%3;
        int iiC = t + (gl_VertexID+2)%3;
        int iA = int(texelFetch(tri_indices,iiA).r);
        int iB = int(texelFetch(tri_indices,iiB).r);
        int iC = int(texelFetch(tri_indices,iiC).r);
        vec3 A = texelFetch(position_vertex, iA).rgb;
        vec3 B = texelFetch(position_vertex, iB).rgb;
        vec3 C = texelFetch(position_vertex, iC).rgb;
        No = cross (B-A, C-A);
        Col = texelFetch(color_face, t).rgb;
        gl_Position = projection_matrix * model_view_matrix * vec4(A,1);
    }
    )";

	const char* fragment_shader_source =
    R"(
    #version 430
    uniform vec3 light_direction;
    in vec3 No;
    in vec3 Col;
    out vec3 frag_out;
    void main()
    {
        vec3 N = normalize(No);
        float lambert = dot(N,light_direction);
        frag_out = lambert*col*
    )";

	load2_bind(vertex_shader_source, fragment_shader_source);

	add_uniforms("light_direction");
}

} // namespace rendering

} // namespace cgogn
