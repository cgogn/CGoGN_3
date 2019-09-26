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

#include <cgogn/rendering/shaders/shader_point_sprite.h>

namespace cgogn
{

namespace rendering
{

ShaderPointSprite* ShaderPointSprite::instance_ = nullptr;
ShaderPointSpriteColor* ShaderPointSpriteColor::instance_ = nullptr;
ShaderPointSpriteSize* ShaderPointSpriteSize::instance_ = nullptr;
ShaderPointSpriteColorSize* ShaderPointSpriteColorSize::instance_ = nullptr;

static const char* vertex_shader_source =
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
"   gl_Position = vec4(vertex_pos,1.0);\n"
"}\n";

static const char* geometry_shader_source =
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
"		corner(posCenter, -1.4,  1.4);\n"
"		corner(posCenter, -1.4, -1.4);\n"
"		corner(posCenter,  1.4,  1.4);\n"
"		corner(posCenter,  1.4, -1.4);\n"
"		EndPrimitive();\n"
"	}\n"
"}\n";

static const char* fragment_shader_source =
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


ShaderPointSprite::ShaderPointSprite()
{
	std::string bs("#version 150\n#define WITH_COLOR 0\n#define WITH_SIZE 0\n");

	std::string vs = bs + std::string(vertex_shader_source);
	std::string gs = bs + std::string(geometry_shader_source);
	std::string fs = bs + std::string(fragment_shader_source);

	load3_bind(vs, fs, gs, "vertex_pos");
	add_uniforms("color", "ambiant", "lightPos", "point_size", "plane_clip", "plane_clip2");
}

ShaderPointSpriteColor::ShaderPointSpriteColor()
{
	std::string bs("#version 150\n#define WITH_COLOR 1\n#define WITH_SIZE 0\n");

	std::string vs = bs + std::string(vertex_shader_source);
	std::string gs = bs + std::string(geometry_shader_source);
	std::string fs = bs + std::string(fragment_shader_source);

	load3_bind(vs, fs, gs, "vertex_pos", "vertex_col");
	add_uniforms("ambiant", "lightPos", "point_size", "plane_clip", "plane_clip2");
}

ShaderPointSpriteSize::ShaderPointSpriteSize()
{
	std::string bs("#version 150\n#define WITH_COLOR 0\n#define WITH_SIZE 1\n");

	std::string vs = bs + std::string(vertex_shader_source);
	std::string gs = bs + std::string(geometry_shader_source);
	std::string fs = bs + std::string(fragment_shader_source);

	load3_bind(vs, fs, gs, "vertex_pos", "vertex_size");
	add_uniforms("color", "ambiant", "lightPos", "plane_clip", "plane_clip2");
}

ShaderPointSpriteColorSize::ShaderPointSpriteColorSize()
{
	std::string bs("#version 150\n#define WITH_COLOR 1\n#define WITH_SIZE 1\n");

	std::string vs = bs + std::string(vertex_shader_source);
	std::string gs = bs + std::string(geometry_shader_source);
	std::string fs = bs + std::string(fragment_shader_source);

	load3_bind(vs, fs, gs, "vertex_pos", "vertex_col", "vertex_size");
	add_uniforms("ambiant", "lightPos", "plane_clip", "plane_clip2");
}

} // namespace rendering

} // namespace cgogn
