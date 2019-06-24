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


#include <cgogn/rendering_pureGL/shaders/shader_text.h>

namespace cgogn
{

namespace rendering_pgl
{


ShaderText* ShaderText::instance_ = nullptr;

ShaderText::ShaderText()
{
	const char* vertex_shader_source =
			"#version 150\n"
			"in vec4 vertex_in;\n"
			"in float char_in;\n"
			"in vec4 colsz_in;\n"
			"uniform mat4 projection_matrix;\n"
			"uniform mat4 model_view_matrix;\n"
			"uniform vec4 italic;\n"
			"out vec2 tc;\n"
			"out vec3 color;\n"
			"const vec2 quad[4] = vec2[](vec2(-0.5,-1.0),vec2(0.5,-1.0),vec2(0.5,1.0),vec2(-0.5,1.));\n"
			"const float tc_u[4] = float[](0,0.009,0.009,0);\n"
			"const float tc_v[4] = float[](1,1,0,0);\n"
			"void main()\n"
			"{\n"
			"	float size = colsz_in.w;\n"
			"   vec4 P4 = model_view_matrix * vec4(vertex_in.xyz,1.0);\n"
			"   P4[0] += size*vertex_in.w;\n"
			"   P4 += vec4(size*(quad[gl_VertexID]+italic[gl_VertexID]),0,0);\n"
			"	tc = vec2( char_in + tc_u[gl_VertexID], tc_v[gl_VertexID]);\n"
			"   gl_Position = projection_matrix * P4;\n"
			"	color = colsz_in.rgb;\n"
			"}\n";

	const char* fragment_shader_source =
			"#version 150\n"
			"out vec3 frag_color;\n"
			"uniform sampler2D texture_unit;\n"
			"in vec3 color;\n"
			"in vec2 tc;\n"
			"void main()\n"
			"{\n"
			"	float a = texture(texture_unit,tc).r;\n"
			"	if (a==0)\n"
			"		discard;\n"
			"	else\n"
			"		frag_color = a*color;\n"
			"}\n";

	load2_bind(vertex_shader_source,fragment_shader_source,
			  "vertex_in", "char_in", "colsz_in");

	add_uniforms("texture_unit","italic");
}

}
} // namespaces
