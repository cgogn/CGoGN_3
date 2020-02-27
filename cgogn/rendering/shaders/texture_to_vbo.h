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

#ifndef CGOGN_RENDERING_SHADERS_TEX2VBO_H_
#define CGOGN_RENDERING_SHADERS_TEX2VBO_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/texture.h>
#include <cgogn/rendering/vbo.h>

#include <cgogn/rendering/shaders/shader_program.h>
#include <cgogn/rendering/shaders/transform_feedback.h>
namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(TEX2VBO, CGOGN_STR(TEX2VBO))

class CGOGN_RENDERING_EXPORT ShaderParamTEX2VBO : public ShaderParam
{
protected:
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(TUin_);
	}

	inline void set_normalization(bool n)
	{
		shader_->set_uniform_value(1, n);
	}

public:
	using LocalShader = ShaderTEX2VBO;

	ShaderParamTEX2VBO(LocalShader* sh) : ShaderParam(sh), TUin_(-1)
	{
	}
	//
	int32_t TUin_;
};

using TFB_TEX2VBO = TransformFeedback<ShaderTEX2VBO>;

inline void texture_to_vbo(TFB_TEX2VBO* tfb, VBO* vbo)
{
	tfb->start(GL_POINTS, {vbo});
}

} // namespace rendering
} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_FLAT_H_
