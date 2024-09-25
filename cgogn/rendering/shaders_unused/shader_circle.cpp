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

#include <cgogn/rendering/shaders/shader_circle.h>

namespace cgogn
{

namespace rendering
{

const char* ShaderCircle::vertex_shader_source_ = "#version 330\n"
												  "precision highp float;\n"
												  "uniform mat4 projection_matrix;\n"
												  "uniform mat4 model_view_matrix;\n"
												  "in vec3 colors;\n"
												  "out vec3 colout;\n"
												  "void main()\n"
												  "{\n"
												  "   float a = 2.0*6.28*float(gl_VertexID)/1000.0;\n"
												  "	float r = 0.7 - a/20.0;\n"
												  "	vec4 P = vec4(r*cos(a), r*sin(a),0.0,1.0);\n"
												  "	gl_Position = projection_matrix * model_view_matrix * P;\n"
												  "	colout = colors;\n"
												  "}\n";

const char* ShaderCircle::fragment_shader_source_ = "#version 330\n"
													"precision highp float;\n"
													"uniform vec4 color;\n"
													"in vec3 colout;\n"
													"out vec4 fragColor;\n"
													"void main()\n"
													"{\n"
													"	fragColor = 0.00001*color+ vec4(colout,1.0);\n"
													"}\n";

ShaderCircle* ShaderCircle::instance_ = nullptr;

std::unique_ptr<ShaderCircle::Param> ShaderCircle::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderCircle();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

void ShaderCircle::set_locations()
{
	glBindAttribLocation(this->id(), 4, "colors");
}

ShaderCircle::ShaderCircle()
{
	this->load(vertex_shader_source_, fragment_shader_source_);
	get_matrices_uniforms();
	unif_color_ = uniform_location("color");
	// default param
	set_color(GLColor(1.0, 1.0, 1., 1.0));
}

void ShaderCircle::set_color(const GLColor& rgba)
{
	set_uniform_value(unif_color_, rgba);
}

} // namespace rendering

} // namespace cgogn
