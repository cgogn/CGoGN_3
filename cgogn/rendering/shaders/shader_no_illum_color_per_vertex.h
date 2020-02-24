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

#ifndef CGOGN_RENDERING_SHADERS_NO_ILLUM_COLOR_PER_VERTEX_H_
#define CGOGN_RENDERING_SHADERS_NO_ILLUM_COLOR_PER_VERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(NoIllumColorPerVertex, CGOGN_STR(NoIllumColorPerVertex))

class CGOGN_RENDERING_EXPORT ShaderParamNoIllumColorPerVertex : public ShaderParam
{
	void set_uniforms() override;

public:
	bool double_side_;

	inline void pick_parameters(const PossibleParameters& pp) override
	{
		double_side_ = pp.double_side_;
	}

	using LocalShader = ShaderNoIllumColorPerVertex;

	ShaderParamNoIllumColorPerVertex(LocalShader* sh) : ShaderParam(sh), double_side_(true)
	{
	}

	inline ~ShaderParamNoIllumColorPerVertex() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif
