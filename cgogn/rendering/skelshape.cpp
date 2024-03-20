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

#include <cgogn/rendering/skelshape.h>

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

ShaderParamSkelShape::ShaderParamSkelShape(ShaderProgram* prg)
	: nbs_(16), color_(0.8f, 0, 0, 1), roughness_(150), shininess_(0.5), light_position_(10, 100, 1000),
	  ShaderParam(prg)
{
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

in vec4 center;

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
	vec4 P4 = model_view_matrix * vec4(center.xyz+center.w*Ps,1);
	Po = P4.xyz;
	No = normal_matrix * Ps;
	gl_Position = projection_matrix * P4;
}
)";

	load2_bind(vertex_shader_source, skel_shape_fragment_shader_source, "center");
	get_uniforms("nbs", "color", "roughness", "shininess", "light_position");
}

void ShaderParamSkelSphere::set_uniforms()
{
	shader_->set_uniforms_values(nbs_, color_, roughness_, shininess_, light_position_);
}

void ShaderParamSkelSphere::draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2 * nbs_ + 2, nbs_ * nbi);
	this->release();
}

SkelSphereDrawer* SkelSphereDrawer::instance_ = nullptr;

SkelSphereDrawer::SkelSphereDrawer()
{
	param_sphere_ = ShaderSkelSphere::generate_param();
	vbo_.set_divisor(param_sphere_->nbs_);
}

SkelSphereDrawer* SkelSphereDrawer::instance()
{
	if (instance_ == nullptr)
		instance_ = new SkelSphereDrawer();
	return instance_;
}

void SkelSphereDrawer::set_subdiv(uint32 ssub)
{
	param_sphere_->set_subdiv(ssub);
	vbo_.set_divisor(param_sphere_->nbs_);
}

void SkelSphereDrawer::set_spheres(const std::vector<GLVec4>& spheres)
{
	update_vbo(spheres, &vbo_);
	param_sphere_->set_vbos({&vbo_});
}

void SkelSphereDrawer::draw(const GLMat4& projection, const GLMat4& view)
{
	param_sphere_->draw_inst(projection, view, vbo_.size());
}


ShaderSkelCone* ShaderSkelCone::instance_ = nullptr;

ShaderSkelCone::ShaderSkelCone()
{
const char* vertex_shader_source = R"(
#version 330
uniform int nbs;
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;

in vec4 P1;
in vec4 P2;
in vec3 N1;
in vec3 N2;

out vec3 Po;
out vec3 No;

const float PI = acos(-1.0);

void main()
{
	float a = 2.0*PI*float(gl_VertexID/2)/float(nbs);
	vec3 Q = cos(a) * N1 + sin(a) * N2;

	vec4 Px = mix(P1,P2,float(gl_VertexID%2));
	vec4 P4 = model_view_matrix * vec4(Px.xyz + Px.w * Q, 1);
	
	Po = P4.xyz;
	float b = (P1.w-P2.w)/length(P2.xyz-P1.xyz);
	No = normal_matrix * normalize(Q + b*cross(N1,N2));
	gl_Position = projection_matrix * P4;
}
)";

	load2_bind(vertex_shader_source, skel_shape_fragment_shader_source, "P1","P2","N1","N2");
	get_uniforms("nbs", "color", "roughness", "shininess", "light_position");
}

void ShaderParamSkelCone::set_uniforms()
{
	shader_->set_uniforms_values(nbs_, color_, roughness_, shininess_, light_position_);
}

void ShaderParamSkelCone::draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2 * nbs_ + 2, nbs_ * nbi);
	this->release();
}

SkelConeDrawer* SkelConeDrawer::instance_ = nullptr;

SkelConeDrawer::SkelConeDrawer()
{
	param_cone_ = ShaderSkelCone::generate_param();
	vbo1_.set_divisor(1);
	vbo2_.set_divisor(1);
	vbo3_.set_divisor(1);
	vbo4_.set_divisor(1);
}

SkelConeDrawer* SkelConeDrawer::instance()
{   
	if (instance_ == nullptr)
		instance_ = new SkelConeDrawer();
	return instance_;
}

void SkelConeDrawer::set_subdiv(uint32 ssub)
{
	param_cone_->set_subdiv(ssub);
}

void SkelConeDrawer::set_cones(const std::vector<GLVec4>& P1, const std::vector<GLVec4>& P2, const std::vector<GLVec3>& N1, const std::vector<GLVec3>& N2)
{
	update_vbo(P1, &vbo1_);
	update_vbo(P2, &vbo2_);	
	update_vbo(N1, &vbo3_);
	update_vbo(N2, &vbo4_);
	param_cone_->set_vbos({&vbo1_,&vbo2_,&vbo3_,&vbo4_});
}

void SkelConeDrawer::draw(const GLMat4& projection, const GLMat4& view)
{
	param_cone_->draw_inst(projection, view, vbo1_.size());
}

void SkelConeDrawer::compute_skel_cone(const GLVec4& A, const GLVec4& B, GLVec4& P1, GLVec4& P2, GLVec3& N1, GLVec3& N2)
{
	GLVec4 AB = B-A;

	if (abs(AB.z())<0.7f)
		N1 = AB.topRows<3>().cross(GLVec3(0,0,1)).normalized();
	else
		N1 = AB.topRows<3>().cross(GLVec3(1,0,0)).normalized();
	N2 = AB.topRows<3>().cross(N1).normalized();



	const GLVec3& AB3 =  AB.topRows<3>();
	float d =  AB3.norm();
	float k = AB[3]/(d*d) ;
	P1.topRows<3>() = A.topRows<3>() - AB3 * (k*A[3]);
	P2.topRows<3>() = B.topRows<3>() - AB3 * (k*B[3]);
	float kk = std::sqrt(1.0f-(AB[3]*AB[3]/(d*d)));
	P1[3] = A[3] * kk;
	P2[3] = B[3] * kk;
	// P2[3] = std::sqrt(B[3]*B[3]-(B[3]*B[3]*AB[3]*AB[3]/(d*d)));


}

} // namespace rendering
} // namespace cgogn