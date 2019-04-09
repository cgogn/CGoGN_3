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

#define CGOGN_RENDER_SHADERS_POINT_SPRITE_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <QOpenGLFunctions>
#include <QColor>

namespace cgogn
{

namespace rendering
{

const char* ShaderPointSpriteGen::vertex_shader_source_ =
"in vec3 vertex_pos;\n"
"#if WITH_COLOR == 1\n"
"in vec3 vertex_col;\n"
"out vec3 color_v;\n"
"#endif\n"
"#if WITH_SIZE == 1\n"
"in float vertex_size;\n"
"out float size_v;\n"
"#endif\n"
"void main()\n"
"{\n"
"	#if WITH_COLOR == 1\n"
"	color_v = vertex_col;\n"
"	#endif\n"
"	#if WITH_SIZE == 1\n"
"	size_v = vertex_size;\n"
"	#endif\n"
"   gl_Position =  vec4(vertex_pos,1.0);\n"
"}\n";

const char* ShaderPointSpriteGen::geometry_shader_source_ =
"layout (points) in;\n"
"layout (triangle_strip, max_vertices=4) out;\n"
"uniform mat4 projection_matrix;\n"
"uniform mat4 model_view_matrix;\n"
"uniform vec4 plane_clip;\n"
"uniform vec4 plane_clip2;\n"
"#if WITH_COLOR == 1\n"
"in vec3 color_v[];\n"
"out vec3 color_f;\n"
"#endif\n"

"#if WITH_SIZE == 1\n"
"in float size_v[];\n"
"out float size_f;\n"
"#else\n"
"uniform float point_size;\n"
"#endif\n"

"out vec2 spriteCoord;\n"
"out vec3 sphereCenter;\n"

"#if ((WITH_COLOR == 1) && (WITH_SIZE == 1)) \n"
"void corner(vec4 center, float x, float y)\n"
"{\n"
"	spriteCoord = vec2(x,y);\n"
"	vec4 pos = center + vec4(size_v[0]*x, size_v[0]*y, 0.0, 0.0);\n"
"	size_f = size_v[0];\n"
"	color_f = color_v[0];\n"
"	gl_Position = projection_matrix *  pos;\n"
"	EmitVertex();\n"
"}\n"
"#endif\n"
"#if ((WITH_COLOR == 1) && (WITH_SIZE == 0)) \n"
"void corner(vec4 center, float x, float y)\n"
"{\n"
"	spriteCoord = vec2(x,y);\n"
"	vec4 pos = center + vec4(point_size*x, point_size*y, 0.0, 0.0);\n"
"	color_f = color_v[0];\n"
"	gl_Position = projection_matrix *  pos;\n"
"	EmitVertex();\n"
"}\n"
"#endif\n"
"#if ((WITH_COLOR == 0) && (WITH_SIZE == 1)) \n"
"void corner(vec4 center, float x, float y)\n"
"{\n"
"	spriteCoord = vec2(x,y);\n"
"	vec4 pos = center + vec4(size_v[0]*x, size_v[0]*y, 0.0, 0.0);\n"
"	size_f = size_v[0];\n"
"	gl_Position = projection_matrix *  pos;\n"
"	EmitVertex();\n"
"}\n"
"#endif\n"
"#if ((WITH_COLOR == 0) && (WITH_SIZE == 0)) \n"
"void corner(vec4 center, float x, float y)\n"
"{\n"
"	spriteCoord = vec2(x,y);\n"
"	vec4 pos = center + vec4(point_size*x, point_size*y, 0.0, 0.0);\n"
"	gl_Position = projection_matrix *  pos;\n"
"	EmitVertex();\n"
"}\n"
"#endif\n"
"void main()\n"
"{\n"
"	float d = dot(plane_clip,gl_in[0].gl_Position);\n"
"	float d2 = dot(plane_clip2,gl_in[0].gl_Position);\n"
"	if ((d<=0.0)&&(d2<=0.0))\n"
"	{\n"
"		vec4 posCenter = model_view_matrix * gl_in[0].gl_Position;\n"
"		sphereCenter = posCenter.xyz;\n"
"		corner(posCenter, -1.4, 1.4);\n"
"		corner(posCenter, -1.4,-1.4);\n"
"		corner(posCenter,  1.4, 1.4);\n"
"		corner(posCenter,  1.4,-1.4);\n"
"		EndPrimitive();\n"
"	}\n"
"}\n";

const char* ShaderPointSpriteGen::fragment_shader_source_ =
"uniform mat4 projection_matrix;\n"
"uniform vec4 ambiant;\n"
"uniform vec3 lightPos;\n"
"#if WITH_SIZE == 1\n"
"in float size_f;\n"
"#else\n"
"uniform float point_size;\n"
"#endif\n"
"#if WITH_COLOR == 1\n"
"in vec3 color_f;\n"
"#else\n"
"uniform vec4 color;\n"
"#endif\n"
"in vec2 spriteCoord;\n"
"in vec3 sphereCenter;\n"
"out vec4 fragColor;\n"

"void main()\n"
"{\n"
"	#if WITH_SIZE == 1\n"
"	float point_size=size_f;\n"
"	#endif\n"
"	vec3 billboard_frag_pos = sphereCenter + vec3(spriteCoord, 0.0) * point_size;\n"
"	vec3 ray_direction = normalize(billboard_frag_pos);\n"
"	float TD = -dot(ray_direction,sphereCenter);\n"
"	float c = dot(sphereCenter, sphereCenter) - point_size * point_size;\n"
"	float arg = TD * TD - c;\n"
"	if (arg < 0.0)\n"
"		discard;\n"
"	float t = -c / (TD - sqrt(arg));\n"
"	vec3 frag_position_eye = ray_direction * t ;\n"
"	vec4 pos = projection_matrix * vec4(frag_position_eye, 1.0);\n"
"	gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;\n"
"	vec3 N = normalize(frag_position_eye - sphereCenter);\n"
"	vec3 L = normalize (lightPos - frag_position_eye);\n"
"	float lambertTerm = dot(N,L);\n"
"	#if WITH_COLOR == 1\n"
"	vec4 result = vec4(color_f*lambertTerm, 1.0);\n"
"	#else\n"
"	vec4 result = vec4(color.rgb*lambertTerm, color.a);\n"
"	#endif\n"
"	result += vec4(ambiant.rgb, 0.0);\n"
"	fragColor = result.rgba;\n"
"}\n";

ShaderPointSpriteGen::ShaderPointSpriteGen(bool color_per_vertex, bool size_per_vertex)
{
	std::string vs("#version 150\n");
	std::string fs("#version 150\n");
	std::string gs("#version 150\n");

	if (color_per_vertex)
	{
		vs += std::string("#define WITH_COLOR 1\n");
		gs += std::string("#define WITH_COLOR 1\n");
		fs += std::string("#define WITH_COLOR 1\n");
	}
	else
	{
		vs += std::string("#define WITH_COLOR 0\n");
		gs += std::string("#define WITH_COLOR 0\n");
		fs += std::string("#define WITH_COLOR 0\n");
	}

	if (size_per_vertex)
	{
		vs += std::string("#define WITH_SIZE 1\n");
		gs += std::string("#define WITH_SIZE 1\n");
		fs += std::string("#define WITH_SIZE 1\n");
	}
	else
	{
		vs += std::string("#define WITH_SIZE 0\n");
		gs += std::string("#define WITH_SIZE 0\n");
		fs += std::string("#define WITH_SIZE 0\n");
	}

	vs += std::string(vertex_shader_source_);
	gs += std::string(geometry_shader_source_);
	fs += std::string(fragment_shader_source_);

	prg_.addShaderFromSourceCode(QOpenGLShader::Vertex, vs.c_str());
	prg_.addShaderFromSourceCode(QOpenGLShader::Geometry, gs.c_str());
	prg_.addShaderFromSourceCode(QOpenGLShader::Fragment, fs.c_str());
	prg_.bindAttributeLocation("vertex_pos", ATTRIB_POS);

	if (color_per_vertex)
		prg_.bindAttributeLocation("vertex_color", ATTRIB_COLOR);

	if (size_per_vertex)
		prg_.bindAttributeLocation("vertex_size", ATTRIB_SIZE);

	prg_.link();
	get_matrices_uniforms();

	unif_color_ = prg_.uniformLocation("color");
	unif_ambiant_ = prg_.uniformLocation("ambiant");
	unif_light_pos_ = prg_.uniformLocation("lightPos");
	unif_size_ = prg_.uniformLocation("point_size");
	unif_plane_clip_ = prg_.uniformLocation("plane_clip");
	unif_plane_clip2_ = prg_.uniformLocation("plane_clip2");


	if (!color_per_vertex)
		set_color(QColor(250, 0, 0));

	set_ambiant(QColor(5, 5, 5));

	if (!size_per_vertex)
		set_size(1.0f);

	set_light_position(QVector3D(10, 10, 1000));
}

void ShaderPointSpriteGen::set_color(const QColor& rgb)
{
	if (unif_color_ >= 0)
		prg_.setUniformValue(unif_color_, rgb);
}

/**
* @brief set ambiant color
* @param rgb
*/
void ShaderPointSpriteGen::set_ambiant(const QColor& rgb)
{
	if (unif_ambiant_ >= 0)
		prg_.setUniformValue(unif_ambiant_, rgb);
}

/**
* @brief set light position relative to screen
* @param l
*/
void ShaderPointSpriteGen::set_light_position(const QVector3D& l)
{
	prg_.setUniformValue(unif_light_pos_, l);
}

/**
* @brief set light position relative to world
* @param l
* @param view_matrix
*/
void ShaderPointSpriteGen::set_local_light_position(const QVector3D& l, const QMatrix4x4& view_matrix)
{
	QVector4D loc4 = view_matrix.map(QVector4D(l, 1.0));
	prg_.setUniformValue(unif_light_pos_, QVector3D(loc4) / loc4.w());
}

/**
* @brief set the size of sphere (call before each draw)
* @param w size ofs phere
*/
//	template <typename std::enable_if<!SPV>::type* = nullptr>
void ShaderPointSpriteGen::set_size(float32 w)
{
	if (unif_size_ >= 0)
		prg_.setUniformValue(unif_size_, w);
}

void ShaderPointSpriteGen::set_plane_clip(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip_, plane);
}

void ShaderPointSpriteGen::set_plane_clip2(const QVector4D& plane)
{
	prg_.setUniformValue(unif_plane_clip2_, plane);
}


template class CGOGN_RENDERING_EXPORT ShaderPointSpriteTpl<false, false>;
template class CGOGN_RENDERING_EXPORT ShaderPointSpriteTpl<true, false>;
template class CGOGN_RENDERING_EXPORT ShaderPointSpriteTpl<false, true>;
template class CGOGN_RENDERING_EXPORT ShaderPointSpriteTpl<true, true>;
template class CGOGN_RENDERING_EXPORT ShaderParamPointSprite<false, false>;
template class CGOGN_RENDERING_EXPORT ShaderParamPointSprite<true, false>;
template class CGOGN_RENDERING_EXPORT ShaderParamPointSprite<false, true>;
template class CGOGN_RENDERING_EXPORT ShaderParamPointSprite<true, true>;

} // namespace rendering

} // namespace cgogn
