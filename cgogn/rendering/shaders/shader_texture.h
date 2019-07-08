/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_RENDERING_SHADERS_TEXTURE_H_
#define CGOGN_RENDERING_SHADERS_TEXTURE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>
#include <cgogn/rendering/texture.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(Texture)

class CGOGN_RENDERING_EXPORT ShaderParamTexture : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(texture_->bind(unit_));
	}

public:

	Texture2D* texture_;
	GLuint unit_;

	using LocalShader = ShaderTexture;

	ShaderParamTexture(LocalShader* sh) :
		ShaderParam(sh),
		unit_(0)
	{}

	inline ~ShaderParamTexture() override {}

	inline void set_vbos(VBO* vbo_pos, VBO* vbo_tc)
	{
		bind_vao();
		associate_vbos(vbo_pos, vbo_tc);
		release_vao();
	}
};

} // namespace rendering

} // namespace cgogn

#endif
