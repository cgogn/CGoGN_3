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

<<<<<<< HEAD:cgogn/rendering/shaders_new/bold_line.h
#ifndef CGOGN_RENDERING_SHADERS_BOLDLINE_H_
#define CGOGN_RENDERING_SHADERS_BOLDLINE_H_
=======
#ifndef CGOGN_RENDERING_SHADERS_NoIllum_COLOR_PER_FACE_H_
#define CGOGN_RENDERING_SHADERS_NoIllum_COLOR_PER_FACE_H_
>>>>>>> 3f10aeda972891997bbe3ed8dacd36fbc4392f95:cgogn/rendering/shaders/shader_no_illum_color_per_face.h

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{
DECLARE_SHADER_CLASS(BoldLine,CGOGN_STR(BoldLine))

<<<<<<< HEAD:cgogn/rendering/shaders_new/bold_line.h
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

	template<typename ...Args>
	void fill(Args&&... args)
	{
		auto a = std::forward_as_tuple(args...);
		color_ = std::get<0>(a);
		width_ = std::get<1>(a);
		blending_ = std::get<2>(a);
	}

	using LocalShader = ShaderBoldLine;

	ShaderParamBoldLine(LocalShader* sh)
		: ShaderParam(sh), color_(color_line_default), width_(2),blending_(true),
		plane_clip_(0, 0, 0, 0), plane_clip2_(0, 0, 0, 0)
	{
	}

	inline ~ShaderParamBoldLine() override
	{
	}


=======
DECLARE_SHADER_CLASS(NoIllumColorPerFace, true, CGOGN_STR(NoIllumColorPerFace))

class CGOGN_RENDERING_EXPORT ShaderParamNoIllumColorPerFace : public ShaderParam
{
	void set_uniforms() override;

public:
	std::array<VBO*, 2> vbos_;
	bool double_side_;

	inline void pick_parameters(const PossibleParameters& pp) override
	{
		double_side_ = pp.double_side_;
	}

	using LocalShader = ShaderNoIllumColorPerFace;

	ShaderParamNoIllumColorPerFace(LocalShader* sh) : ShaderParam(sh), double_side_(true)
	{
		for (auto& v : vbos_)
			v = nullptr;
	}

	inline ~ShaderParamNoIllumColorPerFace() override
	{
	}

	inline VBO** vbo_tb(uint32 i) override
	{
		return &vbos_[i];
	}
>>>>>>> 3f10aeda972891997bbe3ed8dacd36fbc4392f95:cgogn/rendering/shaders/shader_no_illum_color_per_face.h
};

} // namespace rendering
} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_NoIllum_H_
