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

#ifndef CGOGN_RENDERING_SHADERS_COMPUTE_NORMALS_H_
#define CGOGN_RENDERING_SHADERS_COMPUTE_NORMALS_H_

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

DECLARE_SHADER_CLASS(ComputeNormal1, CGOGN_STR(ComputeNormal1))
DECLARE_SHADER_CLASS(ComputeNormal2, CGOGN_STR(ComputeNormal2))

class CGOGN_RENDERING_EXPORT ShaderParamComputeNormal1 : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	VBO* vbo_pos_;
	int32 height_tex_;

	using LocalShader = ShaderComputeNormal1;

	ShaderParamComputeNormal1(LocalShader* sh) : ShaderParam(sh), vbo_pos_(nullptr), height_tex_(0)
	{
	}
};

class CGOGN_RENDERING_EXPORT ShaderParamComputeNormal2 : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	Texture2D* tex_;

	using LocalShader = ShaderComputeNormal2;

	ShaderParamComputeNormal2(LocalShader* sh) : ShaderParam(sh)
	{
	}
};

using TFB_ComputeNormal = TransformFeedback<ShaderComputeNormal2>;

class ComputeNormalEngine
{
	static ComputeNormalEngine* instance_;
	Texture2D* tex_;
	FBO* fbo_;
	std::unique_ptr<ShaderComputeNormal1::Param> param1_;
	std::unique_ptr<ShaderComputeNormal2::Param> param2_;
	TFB_ComputeNormal* tfb_;

	ComputeNormalEngine();

public:
	inline static ComputeNormalEngine* generate()
	{
		if (instance_ == nullptr)
			instance_ = new ComputeNormalEngine();
		return instance_;
	}

	~ComputeNormalEngine();
	void compute(VBO* pos, MeshRender* renderer, VBO* normals);
};

} // namespace rendering
} // namespace cgogn

#endif
