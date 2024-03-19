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

#include <cgogn/rendering/shape_drawer.h>

namespace cgogn
{

namespace rendering
{
const char* skel_shape_fragment_shader_source = R"(
#version 330

uniform vec4 color;
uniform float roughness;
uniform float shininess;
uniform vec3 light_position;

in vec3 Po;
in vec3 No;
out vec4 frag_out;

void main()
{
	vec3 N = normalize(No);
	vec3 L = normalize(light_position-Po);
	float lamb = 0.2+0.8*max(0.0,dot(N,L));
	vec3 E = normalize(-Po);
	vec3 R = reflect(-L,N);
	float spec = pow(max(0.0,dot(E,R)),roughness);

	vec3 colamb = color.rgb*lamb;
	vec3 specol = mix(colamb,vec3(1),shininess);
	vec3 renderColor = mix(colamb,specol,spec);
	frag_out = vec4(renderColor, 1);
}
)";


// ShaderSkelCone* ShaderSkelCone::instance_ = nullptr;

// ShaderSkelCone::ShaderSkelCone()
// {
// const char* vertex_shader_source = R"(
// #version 330
// uniform mat4 projection_matrix;
// uniform mat4 model_view_matrix;
// uniform mat3 normal_matrix;
// uniform int nb;
// out vec3 Po;
// out vec3 No;
// const float PI = acos(-1.0);
// void main()
// {
// 	int iid = gl_InstanceID-1;
// 	int id_h = gl_VertexID%2;
// 	float a = 2.0*PI*float(gl_VertexID/2)/float(nb);
// 	float h = (iid==0) ? 2.0*float(id_h) - 1.0 : float(iid);
// 	vec2 Pc = (iid==0)||(id_h==0)?vec2(cos(a),sin(a)):vec2(0);
// 	vec4 Pc4 = vec4(Pc,h,1);
// 	vec4 P4 = model_view_matrix * Pc4;
// 	Po = P4.xyz;
// 	if (iid ==0)
// 		No = normal_matrix * vec3(Pc,0);
// 	else 
// 	 	No = normal_matrix[2]*float(iid);
// 	gl_Position = projection_matrix * P4;
// }

// )";

// 	load2_bind(vertex_shader_source, shape_fragment_shader_source);
// 	get_uniforms("nb","color","roughness","shininess","light_position");
// }

ShaderParamSkelShape::ShaderParamSkelShape(ShaderProgram* prg):
nbs_(16),color_(0.8f, 0, 0, 1), roughness_(150), shininess_(0.5), light_position_(10, 100, 1000)
{}

void ShaderParamSkelShape::set_uniforms()
{
	shader_->set_uniforms_values(nb_, color_, roughness_, shininess_, light_position_);
}


ShaderSkelSphere* ShaderSkelSphere::instance_ = nullptr;

ShaderSkelSphere::ShaderSkelSphere()
{
	const char* vertex_shader_source = R"(
#version 330
uniform int nbs;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;

in vec4 centers;

out vec3 Po;
out vec3 No;

const float PI = acos(-1.0);

void main()
{
	float db = PI/float(nbs);
	float b = -PI/2.0 + db*float(gl_InstanceID%nbs + gl_VertexID%2);
	float a = 2.0*PI*float(gl_VertexID/2)/float(nbs);
	vec2 Cir = vec2(cos(a),sin(a));
	vec3 Ps = vec3(cos(b)*Cir,sin(b));
	vec4 P4 = model_view_matrix * vec4(Ps,1);
	Po = P4.xyz;
	No = normal_matrix * Ps;
	gl_Position = projection_matrix * P4;
}
)";

	load2_bind(vertex_shader_source, shape_fragment_shader_source,"centers");
	get_uniforms("nbs", "color", "roughness", "shininess", "light_position");
}

void ShaderParamSphere::draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2*nbs_+2, nbs_*nbi);	
	this->release();


}

