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

#ifndef CGOGN_RENDERING_SHADERS_NO_ILLUM_H_
#define CGOGN_RENDERING_SHADERS_NO_ILLUM_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(NoIllum, false, CGOGN_STR(NoIllum))

class CGOGN_RENDERING_EXPORT ShaderParamNoIllum : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	GLColor color_;
	bool double_side_;

	inline void pick_parameters(const PossibleParameters& pp) override
	{
		color_ = pp.color_;
		double_side_ = pp.double_side_;
	}

	using ShaderType = ShaderNoIllum;

	ShaderParamNoIllum(ShaderType* sh) : ShaderParam(sh), color_(color_line_default)
	{
	}

	inline ~ShaderParamNoIllum() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif
