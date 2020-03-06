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

#ifndef CGOGN_RENDERING_FRAME_MANIP_DRAWER_H_
#define CGOGN_RENDERING_FRAME_MANIP_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(Axis, false, CGOGN_STR(Axis))

class CGOGN_RENDERING_EXPORT ShaderParamAxis : public ShaderParam
{
	void set_uniforms() override;

public:
	int selected_;

	inline void pick_parameters(const PossibleParameters&) override
	{
	}

	using LocalShader = ShaderAxis;

	ShaderParamAxis(LocalShader* sh) : ShaderParam(sh), selected_(4)
	{
	}

	inline ~ShaderParamAxis() override
	{
	}
};

DECLARE_SHADER_CLASS(Rings, false, CGOGN_STR(Rings))

class CGOGN_RENDERING_EXPORT ShaderParamRings : public ShaderParam
{
	void set_uniforms() override;

public:
	int selected_;

	inline void pick_parameters(const PossibleParameters&) override
	{
	}

	using LocalShader = ShaderRings;

	ShaderParamRings(LocalShader* sh) : ShaderParam(sh), selected_(4)
	{
	}

	inline ~ShaderParamRings() override
	{
	}
};

class FrameManipDrawer
{
	static FrameManipDrawer* instance_;
	std::unique_ptr<ShaderRings::Param> param_rings_;

	FrameManipDrawer();

public:
	inline static FrameManipDrawer* generate()
	{
		if (instance_ == nullptr)
			instance_ = new FrameManipDrawer();
		return instance_;
	}

	inline void set_selected(int32 s)
	{
		param_rings_->selected_ = s;
	}

	~FrameManipDrawer();

	void draw(const GLMat4& projection, const GLMat4& view, const GLMat4& frame);
};

} // namespace rendering
} // namespace cgogn

#endif
