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

#ifndef CGOGN_RENDERING_SHADERS_COMPUTER_VOLUME_CENTERS_H_
#define CGOGN_RENDERING_SHADERS_COMPUTER_VOLUME_CENTERS_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/ebo.h>
#include <cgogn/rendering/fbo.h>
#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/shaders/shader_program.h>
#include <cgogn/rendering/shaders/transform_feedback.h>
#include <cgogn/rendering/texture.h>
#include <cgogn/rendering/vbo.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(ComputeCenter1, true, CGOGN_STR(ComputeCenter1))
DECLARE_SHADER_CLASS(ComputeCenter2, false, CGOGN_STR(ComputeCenter2))

class CGOGN_RENDERING_EXPORT ShaderParamComputeCenter1 : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	VBO* vbo_pos_;
	int32 height_tex_;

	using ShaderType = ShaderComputeCenter1;

	ShaderParamComputeCenter1(ShaderType* sh) : ShaderParam(sh), vbo_pos_(nullptr), height_tex_(0)
	{
	}
};

class CGOGN_RENDERING_EXPORT ShaderParamComputeCenter2 : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	Texture2D* tex_;

	using ShaderType = ShaderComputeCenter2;

	ShaderParamComputeCenter2(ShaderType* sh) : ShaderParam(sh)
	{
	}
};

using TFB_ComputeCenter = TransformFeedback<ShaderComputeCenter2>;

class ComputeCenterEngine
{
	Texture2D* tex_;
	FBO* fbo_;
	std::unique_ptr<ShaderComputeCenter1::Param> param1_;
	std::unique_ptr<ShaderComputeCenter2::Param> param2_;
	TFB_ComputeCenter* tfb_;

public:
	ComputeCenterEngine();

	~ComputeCenterEngine();

	template <typename MESH, typename ATT>
	inline void check_primitives(rendering::MeshRender* render, const MESH& m, const ATT* vert_pos)
	{
		if (!render->is_primitive_uptodate(rendering::VOLUMES_VERTICES))
			render->init_primitives(m, rendering::VOLUMES_VERTICES, vert_pos);
	}

	void compute(VBO* pos, MeshRender* renderer, VBO* centers);
};

} // namespace rendering
} // namespace cgogn

#endif
