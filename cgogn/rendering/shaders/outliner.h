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

#ifndef CGOGN_RENDERING_SHADERS_OUTLINER_H_
#define CGOGN_RENDERING_SHADERS_OUTLINER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/fbo.h>
#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/shaders/shader_program.h>
#include <cgogn/rendering/texture.h>

namespace cgogn
{

namespace rendering
{

namespace outline_shaders
{

DECLARE_SHADER_CLASS(Mask, false, CGOGN_STR(Mask))

class CGOGN_RENDERING_EXPORT ShaderParamMask : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	using ShaderType = ShaderMask;

	ShaderParamMask(ShaderType* sh) : ShaderParam(sh)
	{
	}

	inline ~ShaderParamMask() override
	{
	}
};

DECLARE_SHADER_CLASS(Sobel, false, CGOGN_STR(Sobel))

class CGOGN_RENDERING_EXPORT ShaderParamSobel : public ShaderParam
{
	void set_uniforms() override;

public:
	Texture2D* texture_;

	using ShaderType = ShaderSobel;

	ShaderParamSobel(ShaderType* sh) : ShaderParam(sh)
	{
	}

	inline ~ShaderParamSobel() override
	{
	}
};

DECLARE_SHADER_CLASS(Blur, false, CGOGN_STR(Blur))

class CGOGN_RENDERING_EXPORT ShaderParamBlur : public ShaderParam
{
	void set_uniforms() override;

public:
	Texture2D* texture_;
	GLuint pass_;

	using ShaderType = ShaderBlur;

	ShaderParamBlur(ShaderType* sh) : ShaderParam(sh), pass_(0)
	{
	}

	inline ~ShaderParamBlur() override
	{
	}
};

DECLARE_SHADER_CLASS(Colorize, false, CGOGN_STR(Colorize))

class CGOGN_RENDERING_EXPORT ShaderParamColorize : public ShaderParam
{
	void set_uniforms() override;

public:
	Texture2D* texture_blur_;
	Texture2D* texture_mask_;
	GLColor color_;

	using ShaderType = ShaderColorize;

	ShaderParamColorize(ShaderType* sh) : ShaderParam(sh)
	{
	}

	inline ~ShaderParamColorize() override
	{
	}

	inline void draw()
	{
		bind();
		glDrawArrays(GL_TRIANGLES, 0, 3);
		release();
	}
};

} // namespace outline_shaders

class OutLiner
{
	static OutLiner* instance_;
	Texture2D* tex_;
	FBO* fbo_mask_;
	FBO* fbo_blur1_;
	FBO* fbo_blur2_;
	std::unique_ptr<outline_shaders::ShaderMask::Param> param_mask_;
	std::unique_ptr<outline_shaders::ShaderSobel::Param> param_sobel_;
	std::unique_ptr<outline_shaders::ShaderBlur::Param> param_blur_;
	std::unique_ptr<outline_shaders::ShaderColorize::Param> param_colorize_;

	OutLiner();

public:
	inline static OutLiner* generate()
	{
		if (instance_ == nullptr)
			instance_ = new OutLiner();
		return instance_;
	}

	~OutLiner();

	void draw(VBO* position, MeshRender* renderer, const rendering::GLMat4& projection_matrix,
			  const rendering::GLMat4& view_matrix, const GLColor& color);
};

} // namespace rendering

} // namespace cgogn

#endif
