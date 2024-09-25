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
	float spec = pow(max(0.0,dot(E,R)),20.0);

	vec3 colamb = color.rgb*lamb;
	vec3 specol = mix(colamb,vec3(1),0.1);
	vec3 renderColor = mix(colamb,specol,spec);
	frag_out = vec4(renderColor, 1);
}
)";


ShaderParamSkelShape::ShaderParamSkelShape(ShaderProgram* prg)
	: nbs_(16), color_(0.8f, 0, 0, 1), light_position_(10, 100, 1000),
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
	get_uniforms("nbs", "color", "light_position");
}

void ShaderParamSkelSphere::set_uniforms()
{
	shader_->set_uniforms_values(nbs_, color_, light_position_);
}

void ShaderParamSkelSphere::draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2 * nbs_ + 2, nbs_ * nbi);
	this->release();
}


SkelSphereDrawer::SkelSphereDrawer()
{
	param_ = ShaderSkelSphere::generate_param();
	vbo_.set_divisor(param_->nbs_);
}

void SkelSphereDrawer::set_vertices_spheres(const std::vector<GLVec4>& spheres)
{
	update_vbo(spheres, &vbo_);
	param_->set_vbos({&vbo_});
}

void SkelSphereDrawer::draw(const GLMat4& projection, const GLMat4& view)
{
	param_->draw_inst(projection, view, vbo_.size());
}

void SkelSphereDrawer::set_subdiv(uint32 ssub)
{
	param_->set_subdiv(ssub);
	vbo_.set_divisor(param_->nbs_);
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
	get_uniforms("nbs", "color", "light_position");
}

void ShaderParamSkelCone::set_uniforms()
{
	shader_->set_uniforms_values(nbs_, color_, light_position_);
}

void ShaderParamSkelCone::draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2 * nbs_ + 2, nbs_ * nbi);
	this->release();
}

SkelConeDrawer::SkelConeDrawer()
{
	param_ = ShaderSkelCone::generate_param();
	vbo1_.set_divisor(1);
	vbo2_.set_divisor(1);
	vbo3_.set_divisor(1);
	vbo4_.set_divisor(1);
}



void SkelConeDrawer::set_real_cones(const std::vector<GLVec4>& P1, const std::vector<GLVec4>& P2, const std::vector<GLVec3>& N1, const std::vector<GLVec3>& N2)
{
	update_vbo(P1, &vbo1_);
	update_vbo(P2, &vbo2_);	
	update_vbo(N1, &vbo3_);
	update_vbo(N2, &vbo4_);
	param_->set_vbos({&vbo1_,&vbo2_,&vbo3_,&vbo4_});
}

//void SkelConeDrawer::set_cones(const std::vector<GLVec4>& P1, const std::vector<GLVec4>& P2)
//{
//	std::vector<GLVec4> A;
//	std::vector<GLVec4> B;
//	std::vector<GLVec3> N1;
//	std::vector<GLVec3> N2;
//	std::size_t n = P1.size();
//	A.resize(n);
//	B.resize(n);
//	N1.resize(n);
//	N2.resize(n);
//	for (std::size_t i = 0; i < n; ++i)
//	{
//		compute_skel_cone(P1[i], P2[i], A[i], B[i], N1[i], N2[i]);
//	}
//	set_real_cones(A, B, N1, N2);
//
//}

void SkelConeDrawer::set_edges_cones(const std::vector<GLVec4>& edges)
{
	std::vector<GLVec4> A;
	std::vector<GLVec4> B;
	std::vector<GLVec3> N1;
	std::vector<GLVec3> N2;
	std::size_t n = edges.size()/2u;
	A.resize(n);
	B.resize(n);
	N1.resize(n);
	N2.resize(n);
	for (std::size_t i = 0; i < n; ++i)
		compute_skel_cone(edges[2 * i], edges[2 * i+1], A[i], B[i], N1[i], N2[i]);
	set_real_cones(A, B, N1, N2);
}

void SkelConeDrawer::draw(const GLMat4& projection, const GLMat4& view)
{
	param_->draw_inst(projection, view, vbo1_.size());
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
	float kk = std::sqrt(1.0f-(AB[3]*k));
	P1[3] = A[3] * kk ;
	P2[3] = B[3] * kk ;
}

ShaderSkelTri* ShaderSkelTri::instance_ = nullptr;

ShaderSkelTri::ShaderSkelTri()
{
	const char* vertex_shader_source = R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;

in vec3 P;
in vec3 N;

out vec3 Po;
out vec3 No;

void main()
{
	vec4 P4 = model_view_matrix * vec4(P, 1);
	Po = P4.xyz;
	No = normal_matrix * N;
	gl_Position = projection_matrix * P4;
}
)";

	load2_bind(vertex_shader_source, skel_shape_fragment_shader_source, "P", "N");
	get_uniforms("color", "light_position");
}

void ShaderParamSkelTri::set_uniforms()
{
	shader_->set_uniforms_values(color_, light_position_);
}

void ShaderParamSkelTri::draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi)
{
	this->bind(projection, view);
	glDrawArrays(GL_TRIANGLES, 0, nbi);
	this->release();
}


SkelTriDrawer::SkelTriDrawer()
{
	param_ = ShaderSkelTri::generate_param();
	vbo_p_.set_divisor(0);
	vbo_n_.set_divisor(0);
}



void SkelTriDrawer::compute_CDisk(const GLVec4& sphA, const GLVec4& sphB, GLVec3& centerA, GLVec3& centerB)
{
	GLVec4 AB = sphB - sphA;
	const GLVec3& AB3 = AB.topRows<3>();
	float d = AB3.norm();
	float k = AB[3] / (d * d);

	centerA = sphA.topRows<3>() - AB3 * (k * sphA[3]);
	centerB = sphB.topRows<3>() - AB3 * (k * sphB[3]);
}

void SkelTriDrawer::interPPS(const GLVec4& planeA, const GLVec4& planeB, const GLVec4& Sph,
							 GLVec3& I1, GLVec3& I2)
{
	const GLVec3& Na = planeA.topRows<3>();
	const GLVec3& Nb = planeB.topRows<3>();

	GLVec3 U = Na.cross(Nb).normalized();

	float dp = Na.dot(Nb);
	float c1 = (-planeA[3] + planeB[3] * dp);
	float c2 = (-planeB[3] + planeA[3] * dp);
	GLVec3 O = (c1 * Na + c2 * Nb) / (1.0f - dp * dp);

	GLVec3 CO = O - Sph.topRows<3>();
	dp = U.dot(CO);
	float delta = dp * dp - (CO.dot(CO) - Sph[3] * Sph[3]);
	float k1 = -dp + std::sqrt(delta);
	float k2 = -dp - std::sqrt(delta);

	I1 = O + k1 * U;
	I2 = O + k2 * U;
}

void SkelTriDrawer::compute_tri2(const GLVec4& Sph1, const GLVec4& Sph2, const GLVec4& Sph3, GLVec3* I)
{
	std::array<GLVec3, 6> centers;

	compute_CDisk(Sph1, Sph2, centers[1], centers[2]);
	GLVec3 N1 = (Sph2.topRows<3>() - Sph1.topRows<3>()).normalized();
	compute_CDisk(Sph2, Sph3, centers[3], centers[4]);
	GLVec3 N2 = (Sph3.topRows<3>() - Sph2.topRows<3>()).normalized();
	compute_CDisk(Sph3, Sph1, centers[5], centers[0]);
	GLVec3 N3 = (Sph1.topRows<3>() - Sph3.topRows<3>()).normalized();

	interPPS({N1[0], N1[1], N1[2], -N1.dot(centers[1])}, {-N3[0], -N3[1], -N3[2], N3.dot(centers[0])}, Sph1, I[0],
			 I[3]);
	interPPS({N2[0], N2[1], N2[2], -N2.dot(centers[3])}, {-N1[0], -N1[1], -N1[2], N1.dot(centers[2])}, Sph2, I[1],
			 I[5]);
	interPPS({N3[0], N3[1], N3[2], -N3.dot(centers[5])}, {-N2[0], -N2[1], -N2[2], N2.dot(centers[4])}, Sph3, I[2],
			 I[4]);
}



void SkelTriDrawer::set_tris_faces(const std::vector<GLVec4>& pts)
{
	assert(pts.size() % 3 == 0);
	std::size_t n = pts.size()/3;

	std::vector<GLVec3> tris;
	tris.resize(6 * n);
	for (std::size_t i = 0; i < n; ++i)
		compute_tri2(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2], &(tris[6 * i]));

	update_vbo(tris, &vbo_p_);

	n *= 2;
	std::vector<GLVec3> no_tris;
	no_tris.resize(3 * n);

	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t j = 3 * i;
		GLVec3 N = (tris[j + 1] - tris[j]).cross(tris[j + 2] - tris[j]).normalized();
		no_tris[j] = N;
		no_tris[j + 1] = N;
		no_tris[j + 2] = N;
	}
	update_vbo(no_tris, &vbo_n_);

	param_->set_vbos({&vbo_p_, &vbo_n_});
}


void SkelTriDrawer::draw(const GLMat4& projection, const GLMat4& view)
{
	param_->draw_inst(projection, view, vbo_p_.size());
}

} // namespace rendering
} // namespace cgogn