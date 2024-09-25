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

#include <cgogn/rendering/shaders/shader_point_sprite.h>

namespace cgogn
{

namespace rendering
{

ShaderPointSprite* ShaderPointSprite::instance_ = nullptr;

ShaderPointSprite::ShaderPointSprite()
{
	const char* vertex_shader_source = R"(
		#version 150
		in vec3 vertex_position;
		in vec3 clipping_position;
		out vec3 clip_pos_v;

		void main()
		{
			gl_Position = vec4(vertex_position, 1.0);
			clip_pos_v = clipping_position;
		}
	)";

	const char* geometry_shader_source = R"(
		#version 150
		layout (points) in;
		layout (triangle_strip, max_vertices=4) out;
		
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;
		uniform float point_size;
		
		in vec3 clip_pos_v[];

		out vec2 spriteCoord;
		out vec3 sphereCenter;

		void corner(vec4 center, float x, float y)
		{
			spriteCoord = vec2(x, y);
			vec4 pos = center + vec4(point_size * x, point_size * y, 0.0, 0.0);
			gl_Position = projection_matrix *  pos;
			EmitVertex();
		}
		
		void main()
		{
			float d = dot(plane_clip, vec4(clip_pos_v[0],1));
			float d2 = dot(plane_clip2, vec4(clip_pos_v[0],1));
			if (d <= 0.0 && d2 <= 0.0)
			{
				vec4 posCenter = model_view_matrix * gl_in[0].gl_Position;
				sphereCenter = posCenter.xyz;
				corner(posCenter, -1.4,  1.4);
				corner(posCenter, -1.4, -1.4);
				corner(posCenter,  1.4,  1.4);
				corner(posCenter,  1.4, -1.4);
				EndPrimitive();
			}
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform vec4 ambiant;
		uniform vec3 light_position;
		uniform float point_size;
		uniform vec4 color;

		in vec2 spriteCoord;
		in vec3 sphereCenter;

		out vec4 frag_out;

		void main()
		{
			vec3 billboard_frag_pos = sphereCenter + vec3(spriteCoord, 0.0) * point_size;
			vec3 ray_direction = normalize(billboard_frag_pos);
			float TD = -dot(ray_direction, sphereCenter);
			float c = dot(sphereCenter, sphereCenter) - point_size * point_size;
			float arg = TD * TD - c;
			if (arg < 0.0)
				discard;
			float t = -c / (TD - sqrt(arg));
			vec3 frag_position_eye = ray_direction * t ;
			vec4 pos = projection_matrix * vec4(frag_position_eye, 1.0);
			gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;
			vec3 N = normalize(frag_position_eye - sphereCenter);
			vec3 L = normalize (light_position - frag_position_eye);
			float lambertTerm = dot(N, L);
			vec4 result = vec4(color.rgb * lambertTerm, color.a);
			result += vec4(ambiant.rgb, 0.0);
			frag_out = result.rgba;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position",
			   "clipping_position");
	get_uniforms("color", "ambiant", "light_position", "point_size", "plane_clip", "plane_clip2");
}

void ShaderParamPointSprite::set_uniforms()
{
	shader_->set_uniforms_values(color_, ambiant_color_, light_position_, point_size_, plane_clip_, plane_clip2_);
}

ShaderPointSpriteColor* ShaderPointSpriteColor::instance_ = nullptr;

ShaderPointSpriteColor::ShaderPointSpriteColor()
{
	const char* vertex_shader_source = R"(
		#version 150
		in vec3 vertex_position;
		in vec4 vertex_color;
		in vec3 clipping_position;
		out vec3 clip_pos_v;
		out vec4 color_v;
		
		void main()
		{
			color_v = vertex_color;
			gl_Position = vec4(vertex_position, 1.0);
			clip_pos_v = clipping_position;
		}
	)";

	const char* geometry_shader_source = R"(
		#version 150
		layout (points) in;
		layout (triangle_strip, max_vertices=4) out;

		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;
		uniform float point_size;
		
		in vec4 color_v[];
		in vec3 clip_pos_v[];

		out vec4 color_f;
		out vec2 spriteCoord;
		out vec3 sphereCenter;

		void corner(vec4 center, float x, float y)
		{
			spriteCoord = vec2(x, y);
			vec4 pos = center + vec4(point_size * x, point_size * y, 0.0, 0.0);
			color_f = color_v[0];
			gl_Position = projection_matrix *  pos;
			EmitVertex();
		}
		
		void main()
		{
			float d = dot(plane_clip, vec4(clip_pos_v[0],1));
			float d2 = dot(plane_clip2, vec4(clip_pos_v[0],1));
			if (d <= 0.0 && d2 <= 0.0)
			{
				vec4 posCenter = model_view_matrix * gl_in[0].gl_Position;
				sphereCenter = posCenter.xyz;
				corner(posCenter, -1.4,  1.4);
				corner(posCenter, -1.4, -1.4);
				corner(posCenter,  1.4,  1.4);
				corner(posCenter,  1.4, -1.4);
				EndPrimitive();
			}
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform vec4 ambiant;
		uniform vec3 light_position;
		uniform float point_size;

		in vec4 color_f;
		in vec2 spriteCoord;
		in vec3 sphereCenter;

		out vec4 frag_out;

		void main()
		{
			vec3 billboard_frag_pos = sphereCenter + vec3(spriteCoord, 0.0) * point_size;
			vec3 ray_direction = normalize(billboard_frag_pos);
			float TD = -dot(ray_direction, sphereCenter);
			float c = dot(sphereCenter, sphereCenter) - point_size * point_size;
			float arg = TD * TD - c;
			if (arg < 0.0)
				discard;
			float t = -c / (TD - sqrt(arg));
			vec3 frag_position_eye = ray_direction * t ;
			vec4 pos = projection_matrix * vec4(frag_position_eye, 1.0);
			gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;
			vec3 N = normalize(frag_position_eye - sphereCenter);
			vec3 L = normalize (light_position - frag_position_eye);
			float lambertTerm = dot(N, L);
			vec4 result = vec4(color_f.rgb * lambertTerm, color_f.a);
			result += vec4(ambiant.rgb, 0.0);
			frag_out = result.rgba;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position", "vertex_color",
			   "clipping_position");
	get_uniforms("ambiant", "light_position", "point_size", "plane_clip", "plane_clip2");
}

void ShaderParamPointSpriteColor::set_uniforms()
{
	shader_->set_uniforms_values(ambiant_color_, light_position_, point_size_, plane_clip_, plane_clip2_);
}

ShaderPointSpriteSize* ShaderPointSpriteSize::instance_ = nullptr;

ShaderPointSpriteSize::ShaderPointSpriteSize()
{
	const char* vertex_shader_source = R"(
		#version 150
		in vec3 vertex_position;
		in float vertex_size;
		in vec3 clipping_position;

		out float size_v;
		out vec3 clip_pos_v;

		void main()
		{
			size_v = vertex_size;
			gl_Position = vec4(vertex_position, 1.0);
			clip_pos_v = clipping_position;
		}
	)";

	const char* geometry_shader_source = R"(
		#version 150
		layout (points) in;
		layout (triangle_strip, max_vertices=4) out;

		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;

		in float size_v[];
		in vec3 clip_pos_v[];

		out float size_f;
		out vec2 spriteCoord;
		out vec3 sphereCenter;

		void corner(vec4 center, float x, float y)
		{
			spriteCoord = vec2(x, y);
			vec4 pos = center + vec4(size_v[0] * x, size_v[0] * y, 0.0, 0.0);
			size_f = size_v[0];
			gl_Position = projection_matrix *  pos;
			EmitVertex();
		}
		
		void main()
		{
			float d = dot(plane_clip, vec4(clip_pos_v[0],1));
			float d2 = dot(plane_clip2, vec4(clip_pos_v[0],1));
			if (d <= 0.0 && d2 <= 0.0)
			{
				vec4 posCenter = model_view_matrix * gl_in[0].gl_Position;
				sphereCenter = posCenter.xyz;
				corner(posCenter, -1.4,  1.4);
				corner(posCenter, -1.4, -1.4);
				corner(posCenter,  1.4,  1.4);
				corner(posCenter,  1.4, -1.4);
				EndPrimitive();
			}
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform vec4 ambiant;
		uniform vec3 light_position;
		uniform vec4 color;

		in float size_f;
		in vec2 spriteCoord;
		in vec3 sphereCenter;
		
		out vec4 frag_out;

		void main()
		{
			float point_size=size_f;
			vec3 billboard_frag_pos = sphereCenter + vec3(spriteCoord, 0.0) * point_size;
			vec3 ray_direction = normalize(billboard_frag_pos);
			float TD = -dot(ray_direction, sphereCenter);
			float c = dot(sphereCenter, sphereCenter) - point_size * point_size;
			float arg = TD * TD - c;
			if (arg < 0.0)
				discard;
			float t = -c / (TD - sqrt(arg));
			vec3 frag_position_eye = ray_direction * t ;
			vec4 pos = projection_matrix * vec4(frag_position_eye, 1.0);
			gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;
			vec3 N = normalize(frag_position_eye - sphereCenter);
			vec3 L = normalize (light_position - frag_position_eye);
			float lambertTerm = dot(N, L);
			vec4 result = vec4(color.rgb * lambertTerm, color.a);
			result += vec4(ambiant.rgb, 0.0);
			frag_out = result.rgba;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position", "vertex_size",
			   "clipping_position");
	get_uniforms("color", "ambiant", "light_position", "plane_clip", "plane_clip2");
}

void ShaderParamPointSpriteSize::set_uniforms()
{
	shader_->set_uniforms_values(color_, ambiant_color_, light_position_, plane_clip_, plane_clip2_);
}

ShaderPointSpriteColorSize* ShaderPointSpriteColorSize::instance_ = nullptr;

ShaderPointSpriteColorSize::ShaderPointSpriteColorSize()
{
	const char* vertex_shader_source = R"(
		#version 150
		in vec3 vertex_position;
		in vec4 vertex_color;
		in vec3 clipping_position;
		out vec4 color_v;
		in float vertex_size;
		out float size_v;
		out vec3 clip_pos_v;

		void main()
		{
			color_v = vertex_color;
			size_v = vertex_size;
			gl_Position = vec4(vertex_position, 1.0);
			clip_pos_v = clipping_position;
		}
	)";

	const char* geometry_shader_source = R"(
		#version 150
		layout (points) in;
		layout (triangle_strip, max_vertices=4) out;
		
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
		uniform vec4 plane_clip;
		uniform vec4 plane_clip2;
		
		in vec4 color_v[];
		in float size_v[];
		in vec3 clip_pos_v[];

		out vec4 color_f;
		out float size_f;
		out vec2 spriteCoord;
		out vec3 sphereCenter;

		void corner(vec4 center, float x, float y)
		{
			spriteCoord = vec2(x, y);
			vec4 pos = center + vec4(size_v[0] * x, size_v[0] * y, 0.0, 0.0);
			size_f = size_v[0];
			color_f = color_v[0];
			gl_Position = projection_matrix *  pos;
			EmitVertex();
		}
		
		void main()
		{
			float d = dot(plane_clip, vec4(clip_pos_v[0],1));
			float d2 = dot(plane_clip2, vec4(clip_pos_v[0],1));
			if (d <= 0.0 && d2 <= 0.0)
			{
				vec4 posCenter = model_view_matrix * gl_in[0].gl_Position;
				sphereCenter = posCenter.xyz;
				corner(posCenter, -1.4,  1.4);
				corner(posCenter, -1.4, -1.4);
				corner(posCenter,  1.4,  1.4);
				corner(posCenter,  1.4, -1.4);
				EndPrimitive();
			}
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform vec4 ambiant;
		uniform vec3 light_position;

		in float size_f;
		in vec4 color_f;
		in vec2 spriteCoord;
		in vec3 sphereCenter;

		out vec4 frag_out;

		void main()
		{
			float point_size=size_f;
			vec3 billboard_frag_pos = sphereCenter + vec3(spriteCoord, 0.0) * point_size;
			vec3 ray_direction = normalize(billboard_frag_pos);
			float TD = -dot(ray_direction, sphereCenter);
			float c = dot(sphereCenter, sphereCenter) - point_size * point_size;
			float arg = TD * TD - c;
			if (arg < 0.0)
				discard;
			float t = -c / (TD - sqrt(arg));
			vec3 frag_position_eye = ray_direction * t ;
			vec4 pos = projection_matrix * vec4(frag_position_eye, 1.0);
			gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;
			vec3 N = normalize(frag_position_eye - sphereCenter);
			vec3 L = normalize (light_position - frag_position_eye);
			float lambertTerm = dot(N, L);
			vec4 result = vec4(color_f.rgb * lambertTerm, color_f.a);
			result += vec4(ambiant.rgb, 0.0);
			frag_out = result.rgba;
		}
	)";

	load3_bind(vertex_shader_source, fragment_shader_source, geometry_shader_source, "vertex_position", "vertex_color",
			   "vertex_size", "clipping_position");
	get_uniforms("ambiant", "light_position", "plane_clip", "plane_clip2");
}

void ShaderParamPointSpriteColorSize::set_uniforms()
{
	shader_->set_uniforms_values(ambiant_color_, light_position_, plane_clip_, plane_clip2_);
}

} // namespace rendering

} // namespace cgogn
