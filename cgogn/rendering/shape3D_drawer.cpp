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

#include <cgogn/rendering/shape3D_drawer.h>

namespace cgogn
{

namespace rendering
{

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

const char* fragment_shader_source = R"(
#version 330

uniform vec3 color;
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

	vec3 colamb = color*lamb;
	vec3 specol = mix(colamb,vec3(1),shininess);
	vec3 renderColor = mix(colamb,specol,spec);
	frag_out = vec4(renderColor, 1);
}
)";

	load2_bind(vertex_shader_source, fragment_shader_source);
	get_uniforms("nb","color","roughness","shininess","light_position");
}

ShaderParamCylinder::ShaderParamCylinder(ShaderType* sh)
	: ShaderParam(sh), nb_(32), color_(0.8, 0, 0, 1), roughness_(150), shininess_(0.5), light_position_(10,50,100)
{}


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

	const char* fragment_shader_source = R"(
#version 330

uniform vec3 color;
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

	vec3 colamb = color*lamb;
	vec3 specol = mix(colamb,vec3(1),shininess);
	vec3 renderColor = mix(colamb,specol,spec);
	frag_out = vec4(renderColor, 1);
}
)";

	load2_bind(vertex_shader_source, fragment_shader_source);
	get_uniforms("nb", "color", "roughness", "shininess", "light_position");
}

ShaderParamSphere::ShaderParamSphere(ShaderType* sh)
	: ShaderParam(sh), nb_(32), color_(0.8, 0, 0, 1), roughness_(150), shininess_(0.5), light_position_(10, 50, 100)
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
	vec4 P4 = model_view_matrix * ((gl_VertexID%2 == 0) ? vec4(0,0,float(gl_InstanceID),1) : vec4(Pc2,0,1));
	No = (gl_InstanceID==0) ? -normal_matrix[2] : normal_matrix * vec3(Pc2,1);
	Po = P4.xyz;
	gl_Position = projection_matrix * P4;
}
)";

	const char* fragment_shader_source = R"(
#version 330

uniform vec3 color;
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

	vec3 colamb = color*lamb;
	vec3 specol = mix(colamb,vec3(1),shininess);
	vec3 renderColor = mix(colamb,specol,spec);
	frag_out = vec4(renderColor, 1);
}
)";

	load2_bind(vertex_shader_source, fragment_shader_source);
	get_uniforms("nb", "color", "roughness", "shininess", "light_position");
}



ShaderParamCone::ShaderParamCone(ShaderType* sh)
	: ShaderParam(sh), nb_(32), color_(0.8, 0, 0, 1), roughness_(150), shininess_(0.5), light_position_(10, 50, 100)
{
}

ShaderParamCone::~ShaderParamCone()
{}

void ShaderParamCone::set_uniforms()
{
	shader_->set_uniforms_values(nb_, color_, roughness_, shininess_, light_position_);
}


Shape3DDrawer::Shape3DDrawer() : nb_subd_(32)
{
	param_cylinder_ = ShaderCylinder::generate_param();
	param_sphere_ = ShaderSphere::generate_param();
	param_cone_ = ShaderCone::generate_param();

	param_cylinder_->nb_ = nb_subd_;
	param_sphere_->nb_ = nb_subd_;
	param_cone_->nb_ = nb_subd_;
}

Shape3DDrawer* Shape3DDrawer::instance_ = nullptr;

Shape3DDrawer::~Shape3DDrawer()
{
}

Shape3DDrawer* Shape3DDrawer::instance()
{
	if (instance_ == nullptr)
		instance_ = new Shape3DDrawer();
	return instance_;
}

void Shape3DDrawer::update_material_cylinder(const GLColor& col, float32 roughness, float32 shininess)
{
	param_cylinder_->color_ = col;
	param_cylinder_->roughness_ = roughness;
	param_cylinder_->shininess_ = shininess;
}

void Shape3DDrawer::update_material_sphere(const GLColor& col, float32 roughness, float32 shininess)
{
	param_sphere_->color_ = col;
	param_sphere_->roughness_ = roughness;
	param_sphere_->shininess_ = shininess;		}

void Shape3DDrawer::update_material_cone(const GLColor& col, float32 roughness, float32 shininess)
{
	param_cone_->color_ = col;
	param_cone_->roughness_ = roughness;
	param_cone_->shininess_ = shininess;
}


void Shape3DDrawer::draw_cylinder(const GLMat4& projection, const GLMat4& view, const GLMat4& transfo)
{
	param_cylinder_->bind(projection, view * transfo);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2*nb_subd_+2, 3);
	param_cylinder_.release();
}

void Shape3DDrawer::draw_sphere(const GLMat4& projection, const GLMat4& view, const GLMat4& transfo)
{
	param_sphere_->bind(projection, view * transfo);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2*nb_subd_+2, nb_subd_);
	param_sphere_.release();
}

void Shape3DDrawer::draw_cone(const GLMat4& projection, const GLMat4& view, const GLMat4& transfo)
{
	param_cone_->bind(projection, view * transfo);
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 2*nb_subd_+2, 2);
	param_cone_.release();
}


void Shape3DDrawer::update_light_position(const GLVec3& lp)
{
	param_cylinder_->light_position_ = lp;
	param_sphere_->light_position_ = lp;
	param_cone_->light_position_ = lp;
}


} // rendering
} // cgogn3
