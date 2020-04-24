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

#ifndef CGOGN_RENDERING_SHADERS_POINT_SPRITE_H_
#define CGOGN_RENDERING_SHADERS_POINT_SPRITE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders_new/shader_program.h>

namespace cgogn
{

namespace rendering
{

namespace shader
{

DECLARE_SHADER_CLASS(PointSprite, CGOGN_STR(PointSprite))

class CGOGN_RENDERING_EXPORT PointSpriteParam : public ShaderParam
{
	void set_uniforms() override;

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	GLColor color_;
	float32 size_;

	template <typename... Args>
	void fill(Args&&... args)
	{
		auto a = std::forward_as_tuple(args...);
		color_ = std::get<0>(a);
		ambiant_color_ = std::get<1>(a);
		light_position_ = std::get<2>(a);
	}

	using LocalShader = PointSprite;

	PointSpriteParam(LocalShader* sh)
		: ShaderParam(sh), color_(color_point_default), ambiant_color_(color_ambiant_default),
		  light_position_(10, 100, 1000), size_(2)
	{
	}

	inline ~PointSpriteParam() override
	{
	}
};

DECLARE_SHADER_CLASS(PointSpriteColor, CGOGN_STR(PointSpriteColor))

class CGOGN_RENDERING_EXPORT PointSpriteColorParam : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(ambiant_color_, light_position_, size_);
	}

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	float32 size_;

	using LocalShader = PointSpriteColor;

	PointSpriteColorParam(LocalShader* sh)
		: ShaderParam(sh), ambiant_color_(color_ambiant_default), light_position_(10, 100, 1000), size_(2)
	{
	}

	inline ~PointSpriteColorParam() override
	{
	}
};

DECLARE_SHADER_CLASS(PointSpriteSize, CGOGN_STR(PointSpriteSize))

class CGOGN_RENDERING_EXPORT PointSpriteSizeParam : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(color_, ambiant_color_, light_position_);
	}

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	GLColor color_;
	float32 size_;

	using LocalShader = PointSpriteSize;

	PointSpriteSizeParam(LocalShader* sh)
		: ShaderParam(sh), color_(color_point_default), ambiant_color_(color_ambiant_default),
		  light_position_(10, 100, 1000)
	{
	}

	inline ~PointSpriteSizeParam() override
	{
	}
};

DECLARE_SHADER_CLASS(PointSpriteColorSize, CGOGN_STR(PointSpriteColorSize))

class CGOGN_RENDERING_EXPORT PointSpriteColorSizeParam : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(ambiant_color_, light_position_);
	}

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;

	using LocalShader = PointSpriteColorSize;

	PointSpriteColorSizeParam(LocalShader* sh)
		: ShaderParam(sh), ambiant_color_(color_ambiant_default), light_position_(10, 100, 1000)
	{
	}

	inline ~PointSpriteColorSizeParam() override
	{
	}
};

} // namespace shader

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_POINT_SPRITE_H_
