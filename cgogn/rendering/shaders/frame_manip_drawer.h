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

DECLARE_SHADER_CLASS(XYGrid, false, CGOGN_STR(XYGrid))

class CGOGN_RENDERING_EXPORT ShaderParamXYGrid : public ShaderParam
{
	void set_uniforms() override;

public:
	float sc_;
	int32 nb_;
	GLColor color_;

	inline void pick_parameters(const PossibleParameters&) override
	{
	}

	using LocalShader = ShaderXYGrid;

	ShaderParamXYGrid(LocalShader* sh) : ShaderParam(sh), sc_(1), nb_(4), color_(1, 1, 1, 1)
	{
	}

	inline ~ShaderParamXYGrid() override
	{
	}
};

class FrameManipDrawer
{
	static FrameManipDrawer* instance_;
	std::unique_ptr<ShaderRings::Param> param_rings_;
	std::unique_ptr<ShaderAxis::Param> param_axis_;
	std::unique_ptr<ShaderXYGrid::Param> param_grid_;

	FrameManipDrawer();

public:
	inline static FrameManipDrawer* generate()
	{
		if (instance_ == nullptr)
			instance_ = new FrameManipDrawer();
		return instance_;
	}

	inline void set_axis_selected(int32 s)
	{
		param_axis_->selected_ = s;
	}
	inline void set_ring_selected(int32 s)
	{
		param_rings_->selected_ = s;
	}
	~FrameManipDrawer();

	void draw_transla(const GLMat4& projection, const GLMat4& view, const GLMat4& frame);
	void draw_rota(const GLMat4& projection, const GLMat4& view, const GLMat4& frame);
	void draw_grid(const GLMat4& projection, const GLMat4& view, const GLMat4& frame, float32 scale = 5);

	inline void draw(const GLMat4& projection, const GLMat4& view, const GLMat4& frame)
	{
		draw_transla(projection, view, frame);
		draw_rota(projection, view, frame);
		draw_grid(projection, view, frame);
	}
};

} // namespace rendering
} // namespace cgogn

#endif
