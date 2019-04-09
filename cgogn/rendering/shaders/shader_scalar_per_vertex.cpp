/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/rendering/shaders/shader_scalar_per_vertex.h>
#include <cgogn/core/utils/unique_ptr.h>

#include <iostream>

namespace cgogn
{

namespace rendering
{

ShaderScalarPerVertex* ShaderScalarPerVertex::instance_ = nullptr;

const char* ShaderScalarPerVertex::vertex_shader_source_ =
"#version 150\n"
"in vec3 vertex_pos;\n"
"in float vertex_scalar;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform float min_value;\n"
"uniform float max_value;\n"
"uniform int color_map;\n"
"uniform int expansion;\n"
"out vec3 color_v;\n"
"out float scalar_v;\n"
"#define M_PI 3.1415926535897932384626433832795\n"
"float scale_and_clamp_to_0_1(float x, float min, float max)\n"
"{\n"
"	float v = (x - min) / (max - min);\n"
"	return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);\n"
"}\n"
"float scale_expand_within_0_1(float x, int n)\n"
"{\n"
"	for (int i = 1; i <= n; i++)\n"
"		x = (1.0 - cos(M_PI * x)) / 2.0;\n"
"	for (int i = -1; i >= n; i--)\n"
"		x = acos(1.0 - 2.0 * x) / M_PI;\n"
"	return x;\n"
"}\n"
"float scale_expand_towards_1(float x, int n)\n"
"{\n"
"	for (int i = 1; i <= n; i++)\n"
"		x = sin(x * M_PI / 2.0);\n"
"	for (int i = -1; i >= n; i--)\n"
"		x = asin(x) * 2.0 / M_PI;\n"
"	return x;\n"
"}\n"
"vec3 color_map_blue_white_red(float x)\n"
"{\n"
"	vec3 c = vec3(0);\n"
"	if (x < 0.0)\n"
"		c.b = 1.0;\n"
"	else if (x < 0.5)\n"
"	{\n"
"		c.r = 2.0 * x;\n"
"		c.g = 2.0 * x;\n"
"		c.b = 1.0;\n"
"	}\n"
"	else if (x < 1.0)\n"
"	{\n"
"		c.r = 1.0;\n"
"		c.g = 2.0 - 2.0 * x;\n"
"		c.b = 2.0 - 2.0 * x;\n"
"	}\n"
"	else\n"
"		c.r = 1.0;\n"
"	return c;\n"
"}\n"
"vec3 color_map_cyan_white_red(float x)\n"
"{\n"
"	if (x < 0.0)\n"
"		return vec3(0.0, 0.0, 1.0) ;\n"
"	if (x < 0.5)\n"
"		return vec3(2.0 * x, 1.0 , 1.0);\n"
"	if (x < 1.0)\n"
"		return vec3(1.0, 2.0 - 2.0 * x, 2.0 - 2.0 * x);\n"
"	return vec3(1.0, 0.0, 0.0) ;\n"
"}\n"
"vec3 color_map_BCGYR(float x)\n"
"{\n"
"	if (x < 0.0)\n"
"		return vec3(0.0, 0.0, 1.0) ;\n"
"	if (x < 0.25)\n"
"		return vec3(0.0, 4.0 * x, 1.0);\n"
"	if (x < 0.5)\n"
"		return vec3(0.0, 1.0 , 2.0 - 4.0 * x);\n"
"	if (x < 0.75)\n"
"		return vec3(4.0 * x - 2.0, 1.0, 0.0);\n"
"	if (x < 1.0)\n"
"		return vec3(1.0, 4.0 - 4.0 * x, 0.0);\n"
"	return vec3(1.0, 0.0, 0.0) ;\n"
"}\n"
"vec3 color_map_blue_green_red(float x)\n"
"{\n"
"	if (x < 0.0)\n"
"		return vec3(0.0, 0.0, 1.0) ;\n"
"	if (x < 0.5)\n"
"		return vec3(0.0, 2.0 * x, 1.0 - 2.0 * x);\n"
"	if (x < 1.0)\n"
"		return vec3(2.0 * x - 1.0, 2.0 - 2.0 * x, 0.0);\n"
"	return vec3(1.0, 0.0, 0.0) ;\n"
"}\n"
"void main()\n"
"{\n"
"	float value =\n"
"		scale_expand_within_0_1(\n"
"			scale_and_clamp_to_0_1(\n"
"				vertex_scalar,\n"
"				min_value,\n"
"				max_value\n"
"			),\n"
"			expansion\n"
"		);\n"
"	switch(color_map)\n"
"	{\n"
"		case 0 : color_v = color_map_blue_white_red(value); break;\n"
"		case 1 : color_v = color_map_cyan_white_red(value); break;\n"
"		case 2 : color_v = color_map_BCGYR(value); break;\n"
"		case 3 : color_v = color_map_blue_green_red(value); break;\n"
"	}\n"
"	scalar_v = value;\n"
"   gl_Position = projection_matrix * model_view_matrix * vec4(vertex_pos, 1.0);\n"
"}\n";

const char* ShaderScalarPerVertex::fragment_shader_source_ =
"#version 150\n"
"in vec3 color_v;\n"
"in float scalar_v;\n"
"uniform bool show_iso_lines;\n"
"uniform int nb_iso_levels;\n"
"out vec3 fragColor;\n"
"void main()\n"
"{\n"
"	if (show_iso_lines)\n"
"	{\n"
"		float s = scalar_v * float(nb_iso_levels);\n"
"		if (s - floor(s) < 0.05)\n"
"			fragColor = vec3(0.0);\n"
"		else\n"
"			fragColor = color_v;\n"
"	}\n"
"	else\n"
"		fragColor = color_v;\n"
"}\n";

ShaderScalarPerVertex::ShaderScalarPerVertex()
{
	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_shader_source_);
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fragment_shader_source_);
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);
	prg_.bindAttributeLocation("vertex_scalar", ATTRIB_SCALAR);
	prg_.link();
	get_matrices_uniforms();
	unif_color_map_ = prg_.uniformLocation("color_map");
	unif_expansion_ = prg_.uniformLocation("expansion");
	unif_min_value_ = prg_.uniformLocation("min_value");
	unif_max_value_ = prg_.uniformLocation("max_value");
	unif_show_iso_lines_ = prg_.uniformLocation("show_iso_lines");
	unif_nb_iso_levels_ = prg_.uniformLocation("nb_iso_levels");
}

void ShaderScalarPerVertex::set_color_map(ColorMap cm)
{
	if (unif_color_map_ >= 0)
		prg_.setUniformValue(unif_color_map_, cm);
}

void ShaderScalarPerVertex::set_expansion(int32 expansion)
{
	if (unif_expansion_ >= 0)
		prg_.setUniformValue(unif_expansion_, expansion);
}

void ShaderScalarPerVertex::set_min_value(float32 value)
{
	if (unif_min_value_ >= 0)
		prg_.setUniformValue(unif_min_value_, value);
}

void ShaderScalarPerVertex::set_max_value(float32 value)
{
	if (unif_max_value_ >= 0)
		prg_.setUniformValue(unif_max_value_, value);
}

void ShaderScalarPerVertex::set_show_iso_lines(bool b)
{
	if (unif_show_iso_lines_ >= 0)
		prg_.setUniformValue(unif_show_iso_lines_, b);
}

void ShaderScalarPerVertex::set_nb_iso_levels(int32 nb)
{
	if (unif_nb_iso_levels_ >= 0)
		prg_.setUniformValue(unif_nb_iso_levels_, nb);
}

std::unique_ptr<ShaderScalarPerVertex::Param> ShaderScalarPerVertex::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderScalarPerVertex();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

} // namespace rendering

} // namespace cgogn
