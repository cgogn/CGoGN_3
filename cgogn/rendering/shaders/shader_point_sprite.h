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

#ifndef CGOGN_RENDERING_SHADER_POINT_SPRITE_H_
#define CGOGN_RENDERING_SHADER_POINT_SPRITE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(PointSprite)

class CGOGN_RENDERING_EXPORT ShaderParamPointSprite : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(color_, ambiant_color_, light_pos_, size_, plane_clip_, plane_clip2_);
	}

public:

	GLColor color_;
	GLColor ambiant_color_;
	GLVec3 light_pos_;
	float32 size_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using LocalShader = ShaderPointSprite;

	ShaderParamPointSprite(LocalShader* sh) :
		ShaderParam(sh),
		color_(color_point_default),
		ambiant_color_(color_ambiant_default),
		light_pos_(10, 100, 1000),
		size_(2),
		plane_clip_(0, 0, 0, 0),
		plane_clip2_(0, 0, 0, 0)
	{}

	inline ~ShaderParamPointSprite() override {}

	inline void set_vbos(VBO* vbo_pos)
	{
		bind_vao();
		associate_vbos(vbo_pos);
		release_vao();
	}
};

DECLARE_SHADER_CLASS(PointSpriteColor)

class CGOGN_RENDERING_EXPORT ShaderParamPointSpriteColor : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(ambiant_color_, light_pos_, size_, plane_clip_, plane_clip2_);
	}

public:

	GLColor ambiant_color_;
	GLVec3 light_pos_;
	float32 size_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using LocalShader = ShaderPointSpriteColor;

	ShaderParamPointSpriteColor(LocalShader* sh) :
		ShaderParam(sh),
		ambiant_color_(color_ambiant_default),
		light_pos_(10, 100, 1000),
		size_(2),
		plane_clip_(0, 0, 0, 0),
		plane_clip2_(0, 0, 0, 0)
	{}

	inline ~ShaderParamPointSpriteColor() override {}

	inline void set_vbos(VBO* vbo_pos, VBO* vbo_col)
	{
		bind_vao();
		associate_vbos(vbo_pos,vbo_col);
		release_vao();
	}
};

DECLARE_SHADER_CLASS(PointSpriteSize)

class CGOGN_RENDERING_EXPORT ShaderParamPointSpriteSize : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(color_, ambiant_color_, light_pos_, plane_clip_, plane_clip2_);
	}

public:

	GLColor color_;
	GLColor ambiant_color_;
	GLVec3 light_pos_;
	float32 size_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using LocalShader = ShaderPointSpriteSize;

	ShaderParamPointSpriteSize(LocalShader* sh) :
		ShaderParam(sh),
		color_(color_point_default),
		ambiant_color_(color_ambiant_default),
		light_pos_(10, 100, 1000),
		plane_clip_(0, 0, 0, 0),
		plane_clip2_(0, 0, 0, 0)
	{}

	inline ~ShaderParamPointSpriteSize() override {}

	inline void set_vbos(VBO* vbo_pos, VBO* vbo_size)
	{
		bind_vao();
		associate_vbos(vbo_pos,vbo_size);
		release_vao();
	}
};

DECLARE_SHADER_CLASS(PointSpriteColorSize)

class CGOGN_RENDERING_EXPORT ShaderParamPointSpriteColorSize : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(ambiant_color_, light_pos_, plane_clip_, plane_clip2_);
	}

public:

	GLColor ambiant_color_;
	GLVec3 light_pos_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using LocalShader = ShaderPointSpriteColorSize;

	ShaderParamPointSpriteColorSize(LocalShader* sh) :
		ShaderParam(sh),
		ambiant_color_(color_ambiant_default),
		light_pos_(10, 100, 1000),
		plane_clip_(0, 0, 0, 0),
		plane_clip2_(0, 0, 0, 0)
	{}

	inline ~ShaderParamPointSpriteColorSize() override {}

	inline void set_vbos(VBO* vbo_pos, VBO* vbo_col, VBO* vbo_size)
	{
		bind_vao();
		associate_vbos(vbo_pos,vbo_col,vbo_size);
		release_vao();
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_POINT_SPRITE_H_
