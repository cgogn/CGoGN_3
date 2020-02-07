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

#ifndef CGOGN_RENDERING_SHADERS_VECTORPERVERTEX_H_
#define CGOGN_RENDERING_SHADERS_VECTORPERVERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(VectorPerVertex,CGOGN_STR(VectorPerVertex))

class CGOGN_RENDERING_EXPORT ShaderParamVectorPerVertex : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(color_, length_);
	}

public:
	GLColor color_;
	float32 length_;

	using LocalShader = ShaderVectorPerVertex;

	ShaderParamVectorPerVertex(LocalShader* sh) : ShaderParam(sh), color_(color_line_default), length_(1.0f)
	{
	}

	inline ~ShaderParamVectorPerVertex() override
	{
	}


};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_VECTORPERVERTEX_H_
