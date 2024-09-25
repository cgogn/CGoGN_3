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
const char* shape_fragment_shader_source = R"(
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



ShaderParamShape::ShaderParamShape(ShaderProgram* prg) : ShaderParam(prg),
	color_(0.8f, 0, 0, 1), roughness_(150), shininess_(0.5), light_position_(10, 100, 1000)
{
}

ShaderParamShapeSub::ShaderParamShapeSub(ShaderProgram* prg)
	: ShaderParamShape(prg), nb_(32)
{
}

ShaderCylinder* ShaderCylinder::instance_ = nullptr;

ShaderCylinder::ShaderCylinder()
{
const char* vertex_shader_source = R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;
uniform int nb;
out vec3 Po;
out vec3 No;
const float PI = acos(-1.0);
void main()
{
	int iid = gl_InstanceID-1;
	int id_h = gl_VertexID%2;
	float a = 2.0*PI*float(gl_VertexID/2)/float(nb);
	float h = (iid==0) ? 2.0*float(id_h) - 1.0 : float(iid);
	vec2 Pc = (iid==0)||(id_h==0)?vec2(cos(a),sin(a)):vec2(0);
	vec4 Pc4 = vec4(Pc,h,1);
	vec4 P4 = model_view_matrix * Pc4;
	Po = P4.xyz;
	if (iid ==0)
		No = normal_matrix * vec3(Pc,0);
	else 
	 	No = normal_matrix[2]*float(iid);
	gl_Position = projection_matrix * P4;
}

)";

//const char* fragment_shader_source = R"(
//#version 330
//
//uniform vec4 color;
//uniform float roughness;
//uniform float shininess;
//uniform vec3 light_position;
//
//in vec3 Po;
//in vec3 No;
//out vec4 frag_out;
//
//void main()
//{
//	vec3 N = normalize(No);
//	vec3 L = normalize(light_position-Po);
//	float lamb = 0.2+0.8*max(0.0,dot(N,L));
//	vec3 E = normalize(-Po);
//	vec3 R = reflect(-L,N);
//	float spec = pow(max(0.0,dot(E,R)),roughness);
//
//	vec3 colamb = color.rgb*lamb;
//	vec3 specol = mix(colamb,vec3(1),shininess);
//	vec3 renderColor = mix(colamb,specol,spec);
//	frag_out = vec4(renderColor, 1);
//}
//)";
//
	load2_bind(vertex_shader_source, shape_fragment_shader_source);
	get_uniforms("nb","color","roughness","shininess","light_position");
}

ShaderParamCylinder::ShaderParamCylinder(ShaderType* sh)
	: ShaderParamShapeSub(sh){}


void ShaderParamCylinder::set_uniforms()
{
	shader_->set_uniforms_values(nb_, color_, roughness_, shininess_, light_position_);
}


ShaderParamCylinder::~ShaderParamCylinder()
{
}

ShaderSphere* ShaderSphere::instance_ = nullptr;

ShaderSphere::ShaderSphere()
{
	const char* vertex_shader_source = R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;
uniform int nb;
out vec3 Po;
out vec3 No;
const float PI = acos(-1.0);
void main()
{
	float db = PI/float(nb);
	float b = -PI/2.0 + db*float(gl_InstanceID+gl_VertexID%2);
	float a = 2.0*PI*float(gl_VertexID/2)/float(nb);
	vec2 Cir = vec2(cos(a),sin(a));
	vec3 Ps = vec3(cos(b)*Cir,sin(b));
	vec4 P4 = model_view_matrix * vec4(Ps,1);
	Po = P4.xyz;
	No = normal_matrix * Ps;
	gl_Position = projection_matrix * P4;
}


)";

	load2_bind(vertex_shader_source, shape_fragment_shader_source);
	get_uniforms("nb", "color", "roughness", "shininess", "light_position");
}

ShaderParamSphere::ShaderParamSphere(ShaderType* sh)
	: ShaderParamShapeSub(sh)
{
}

ShaderParamSphere::~ShaderParamSphere()
{}

void ShaderParamSphere::set_uniforms()
{
	shader_->set_uniforms_values(nb_, color_, roughness_, shininess_, light_position_);
}



ShaderCone* ShaderCone::instance_ = nullptr;

ShaderCone::ShaderCone()
{
	const char* vertex_shader_source = R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;
uniform int nb;
out vec3 Po;
out vec3 No;
const float PI = acos(-1.0);
void main()
{
	float a = 2.0*PI*float(gl_VertexID/2)/float(nb);
	vec2 Pc2 = vec2(cos(a),sin(a));
	vec4 P4 = model_view_matrix * ((gl_VertexID%2 == 0) ? vec4(0,0,2.0*float(gl_InstanceID)-1.0,1) : vec4(Pc2,-1,1));
	No = (gl_InstanceID==0) ? -normal_matrix[2] : normal_matrix * vec3(2.0*Pc2,1);
	Po = P4.xyz;
	gl_Position = projection_matrix * P4;
}
)";

	load2_bind(vertex_shader_source, shape_fragment_shader_source);
	get_uniforms("nb", "color", "roughness", "shininess", "light_position");
}



ShaderParamCone::ShaderParamCone(ShaderType* sh)
	: ShaderParamShapeSub(sh)
{
}

ShaderParamCone::~ShaderParamCone()
{}

void ShaderParamCone::set_uniforms()
{
	shader_->set_uniforms_values(nb_, color_, roughness_, shininess_, light_position_);
}




ShaderCube* ShaderCube::instance_ = nullptr;

ShaderCube::ShaderCube()
{
	const char* vertex_shader_source = R"(
#version 330
uniform mat4 projection_matrix;
uniform mat4 model_view_matrix;
uniform mat3 normal_matrix;
out vec3 Po;
out vec3 No;
const uint f_code[6]=uint[](0406602u,0735531u, 0514410u,0627723u, 0321120u,0756546u);

void main()
{
	int f_id = gl_VertexID/6;
	uint vf = (uint(gl_VertexID)%6u)*3u;
	uint xyz = (f_code[f_id] >> vf )&7u; //uint(ind[gl_VertexID]);
	vec3 Q = vec3(xyz&1u, (xyz>>1u)&1u, (xyz>>2u)&1u) * 2.0 - 1.0;
	vec4 P4 = model_view_matrix * vec4(Q,1);
	Po = P4.xyz;

	vec3 N = vec3(0);
	N[f_id/2]=float(f_id%2)*2.0-1.0;
	No = normal_matrix * N;
	gl_Position = projection_matrix * P4;
}

)";

	load2_bind(vertex_shader_source, shape_fragment_shader_source);
	get_uniforms("color", "roughness", "shininess", "light_position");
}

ShaderParamCube::ShaderParamCube(ShaderType* sh)
	: ShaderParamShape(sh)
{
}

ShaderParamCube::~ShaderParamCube()
{
}

void ShaderParamCube::set_uniforms()
{
	shader_->set_uniforms_values(color_, roughness_, shininess_, light_position_);
}

ShapeDrawer* ShapeDrawer::instance_ = nullptr;


ShapeDrawer::ShapeDrawer() : nb_subd_(32)
{
	param_[SHAPE::CYLINDER] = ShaderCylinder::generate_param();
	param_[SHAPE::SPHERE] = ShaderSphere::generate_param();
	param_[SHAPE::CONE] = ShaderCone::generate_param();
	param_[SHAPE::CUBE] = ShaderCube::generate_param();

	for (auto &p: param_)
		p->set_subdiv(nb_subd_);
}


ShapeDrawer::~ShapeDrawer()
{
}

ShapeDrawer* ShapeDrawer::instance()
{
	if (instance_ == nullptr)
		instance_ = new ShapeDrawer();
	return instance_;
}



void ShapeDrawer::update_material(const GLColor& col, float32 roughness, float32 shininess)
{
	for (auto& p : param_)
	{
		p->color_ = col;
		p->roughness_ = roughness;
		p->shininess_ = shininess;
	};
}

void ShapeDrawer::update_light_position(const GLVec3& lp)
{
	for (auto& p : param_)
		p->light_position_ = lp;
}

void ShapeDrawer::update_subdivision(int32 nb)
{
	this->nb_subd_ = nb;

	for (auto& p : param_)
		p->set_subdiv(nb);
}


void ShapeDrawer::draw(SHAPE s, const GLMat4& projection, const GLMat4& view)
{

	param_[s]->draw(projection, view);
}

void ShaderParamCylinder::draw(const GLMat4& projection, const GLMat4& view)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2*nb_+2, 3);
	this->release();
}

void ShaderParamSphere::draw(const GLMat4& projection, const GLMat4& view)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2*nb_+2, nb_);
	this->release();
}

void ShaderParamCone::draw(const GLMat4& projection, const GLMat4& view)
{
	this->bind(projection, view);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2 * nb_+2, 2);
	this->release();
}

void ShaderParamCube::draw(const GLMat4& projection, const GLMat4& view)
{
	this->bind(projection, view);
	glDrawArrays(GL_TRIANGLES, 0, 36);
	this->release();
}

void ShapeDrawer::drawSphere(const GLMat4& projection, const GLMat4& view, const GLVec3& p, float32 radius)
{
	draw(SPHERE, projection, view * (Eigen::Translation3f(p) * Eigen::Scaling(radius)).matrix());
}

GLMat4 ShapeDrawer::points_2_transfo(const GLVec3& p1, const GLVec3& p2, float32 radius)
{
	GLVec3 dir = p2 - p1;
	float l = dir.norm();
	dir /= l;
	GLVec3 axis = GLVec3(-dir.y(), dir.x(), 0.0);
	GLVec3(0, 0, 1).cross(dir);
	float la = axis.norm();
	Eigen::Affine3f tr = Eigen::Affine3f(Eigen::Translation3f((p1 + p2) / 2));
	if (la != 0.0f)
	{
		axis /= la;
		float alpha = std::asin(la);
		if (dir.z() < 0)
			alpha = M_PI - alpha;
		tr *= Eigen::AngleAxisf(alpha, axis);
	}
	tr *= Eigen::Scaling(radius, radius, l / 2);
	return tr.matrix();
}

void ShapeDrawer::drawCylinder(const GLMat4& projection, const GLMat4& view, const GLVec3& p1, const GLVec3& p2, float32 radius)
{
	draw(CYLINDER, projection, view * points_2_transfo(p1,p2,radius));
}


void ShapeDrawer::drawCone(const GLMat4& projection, const GLMat4& view, const GLVec3& p1, const GLVec3& p2, float32 radius)
{
	draw(CONE, projection, view * points_2_transfo(p1,p2,radius));
}

} // rendering
} // cgogn3
