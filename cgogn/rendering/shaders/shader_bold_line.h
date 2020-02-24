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

#ifndef CGOGN_RENDERING_SHADERS_BOLDLINE_H_
#define CGOGN_RENDERING_SHADERS_BOLDLINE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{
DECLARE_SHADER_CLASS(BoldLine, CGOGN_STR(BoldLine))

class CGOGN_RENDERING_EXPORT ShaderParamBoldLine : public ShaderParam
{
	inline void set_uniforms() override
	{
		int viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		GLVec2 wd(width_ / float32(viewport[2]), width_ / float32(viewport[3]));
		shader_->set_uniforms_values(color_, wd, plane_clip_, plane_clip2_);
	}

public:
	GLColor color_;
	float32 width_;
	bool blending_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	inline void pick_parameters(const PossibleParameters& pp) override
	{
		color_ = pp.color_;
		width_ = pp.width_;
		blending_ = pp.blending_;
		plane_clip_ = pp.plane_clip_;
		plane_clip2_ = pp.plane_clip2_;
	}

	using LocalShader = ShaderBoldLine;

	ShaderParamBoldLine(LocalShader* sh)
		: ShaderParam(sh), color_(color_line_default), width_(2), blending_(true), plane_clip_(0, 0, 0, 0),
		  plane_clip2_(0, 0, 0, 0)
	{
	}

	inline ~ShaderParamBoldLine() override
	{
	}
};

} // namespace rendering
} // namespace cgogn

#endif
