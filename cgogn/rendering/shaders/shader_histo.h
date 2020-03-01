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

#ifndef CGOGN_RENDERING_SHADERS_HISTO_H_
#define CGOGN_RENDERING_SHADERS_HISTO_H_
#include <GL/gl3w.h>
#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/fbo.h>
#include <cgogn/rendering/shaders/shader_program.h>
#include <cgogn/rendering/texture.h>

namespace cgogn
{

namespace rendering
{
DECLARE_SHADER_CLASS(Histo, false, CGOGN_STR(Histo))

class CGOGN_RENDERING_EXPORT ShaderParamHisto : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(texture_->bind(0), texture_->width(), 1.0f);
	}

public:
	Texture2D* tex_fbo_;
	FBO* fbo_;
	Texture2D* texture_;

	using LocalShader = ShaderHisto;

	inline ShaderParamHisto(LocalShader* sh) : ShaderParam(sh)
	{
		tex_fbo_ = new Texture2D();
		tex_fbo_->alloc(1, 1, GL_R32F, GL_RED, nullptr, GL_FLOAT);
		fbo_ = new FBO({tex_fbo_}, false, nullptr);
	}

	inline ~ShaderParamHisto() override
	{
	}

	inline void draw(int nbb, std::vector<float>& histogram)
	{
		histogram.resize(nbb, 5.5f);

		fbo_->resize(nbb, 1);

		bind();
		shader_->set_uniform_value(2, 1.0f - 0.5f / nbb);
		fbo_->bind();

		GLenum idbuf = GL_COLOR_ATTACHMENT0;
		glDrawBuffers(1, &idbuf);
		glClear(GL_COLOR_BUFFER_BIT);
		glClearColor(0.0, 0, 0, 0);
		glViewport(0, 0, nbb, 1);
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE);
		glDrawArrays(GL_POINTS, 0, texture_->width() * texture_->height());
		glDisable(GL_BLEND);
		fbo_->release();
		release();

		fbo_->bind();
		glReadBuffer(GL_COLOR_ATTACHMENT0);
		glReadPixels(0, 0, nbb, 1, GL_RED, GL_FLOAT, histogram.data());
		fbo_->release();

		for (float h : histogram)
			std::cout << "|" << h;
		float tot = 0;
		for (float h : histogram)
			tot += h;
		std::cout << "| => " << tot << std::endl;
	}
};

} // namespace rendering
} // namespace cgogn

#endif
