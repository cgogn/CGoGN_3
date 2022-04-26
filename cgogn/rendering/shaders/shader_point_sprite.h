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

#ifndef CGOGN_RENDERING_SHADER_POINT_SPRITE_H_
#define CGOGN_RENDERING_SHADER_POINT_SPRITE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(PointSprite, false, CGOGN_STR(PointSprite))

class CGOGN_RENDERING_EXPORT ShaderParamPointSprite : public ShaderParam
{
	void set_uniforms() override;

public:
	GLColor color_;
	GLColor ambiant_color_;
	GLVec3 light_position_;
	float32 point_size_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using ShaderType = ShaderPointSprite;

	ShaderParamPointSprite(ShaderType* sh)
		: ShaderParam(sh,true), color_(1, 1, 1, 1), ambiant_color_(0.05f, 0.05f, 0, 1), light_position_(10, 100, 1000),
		  point_size_(2), plane_clip_(0, 0, 0, 0), plane_clip2_(0, 0, 0, 0)
	{
	}

	inline ~ShaderParamPointSprite() override
	{
	}
};

DECLARE_SHADER_CLASS(PointSpriteColor, false, CGOGN_STR(PointSpriteColor))

class CGOGN_RENDERING_EXPORT ShaderParamPointSpriteColor : public ShaderParam
{
	void set_uniforms() override;

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	float32 point_size_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using ShaderType = ShaderPointSpriteColor;

	ShaderParamPointSpriteColor(ShaderType* sh)
		: ShaderParam(sh,true), ambiant_color_(0.05f, 0.05f, 0, 1), light_position_(10, 100, 1000), point_size_(2),
		  plane_clip_(0, 0, 0, 0), plane_clip2_(0, 0, 0, 0)
	{
	}

	inline ~ShaderParamPointSpriteColor() override
	{
	}
};

DECLARE_SHADER_CLASS(PointSpriteSize, false, CGOGN_STR(PointSpriteSize))

class CGOGN_RENDERING_EXPORT ShaderParamPointSpriteSize : public ShaderParam
{
	void set_uniforms() override;

public:
	GLColor color_;
	GLColor ambiant_color_;
	GLVec3 light_position_;
	float32 point_size_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using ShaderType = ShaderPointSpriteSize;

	ShaderParamPointSpriteSize(ShaderType* sh)
		: ShaderParam(sh,true), color_(1, 1, 1, 1), ambiant_color_(0.05f, 0.05f, 0, 1), light_position_(10, 100, 1000),
		  point_size_(1.0f), plane_clip_(0, 0, 0, 0), plane_clip2_(0, 0, 0, 0)
	{
	}

	inline ~ShaderParamPointSpriteSize() override
	{
	}
};

DECLARE_SHADER_CLASS(PointSpriteColorSize, false, CGOGN_STR(PointSpriteColorSize))

class CGOGN_RENDERING_EXPORT ShaderParamPointSpriteColorSize : public ShaderParam
{
	void set_uniforms() override;

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using ShaderType = ShaderPointSpriteColorSize;

	ShaderParamPointSpriteColorSize(ShaderType* sh)
		: ShaderParam(sh,true), ambiant_color_(0.05f, 0.05f, 0, 1), light_position_(10, 100, 1000), plane_clip_(0, 0, 0, 0),
		  plane_clip2_(0, 0, 0, 0)
	{
	}

	inline ~ShaderParamPointSpriteColorSize() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_POINT_SPRITE_H_
