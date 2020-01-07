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

#ifndef CGOGN_RENDERING_FLAT_TR_DR_H_
#define CGOGN_RENDERING_FLAT_TR_DR_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/transparency_shaders/shader_copy_depth.h>
#include <cgogn/rendering/transparency_shaders/shader_transparent_flat.h>
#include <cgogn/rendering/transparency_shaders/shader_transparent_phong.h>
#include <cgogn/rendering/transparency_shaders/shader_transparent_quad.h>
#include <cgogn/rendering/transparency_shaders/shader_transparent_volumes.h>

#include <cgogn/rendering/shaders/vbo.h>

#include <GLColor>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT SurfaceTransparencyDrawer
{
	int max_nb_layers_;

	/// rendering shader
	std::unique_ptr<ShaderFlatTransp::Param> param_flat_;

	std::unique_ptr<ShaderPhongTransp::Param> param_phong_;

	/// shader for quad blending with opaque scene
	std::unique_ptr<ShaderTranspQuad::Param> param_trq_;

	std::unique_ptr<ShaderCopyDepth::Param> param_copy_depth_;

	/// FBO
	std::unique_ptr<QOpenGLFramebufferObject> fbo_layer_;

	/// Occlusion query
	GLuint oq_transp_;

	QOpenGLFunctions_3_3_Core* ogl33_;

	int width_;

	int height_;

	GLuint depthTexture_;

public:
	~SurfaceTransparencyDrawer();

	/**
	 * @brief create and init

	 */
	SurfaceTransparencyDrawer();

	/**
	 * @brief resize call_back need to be called when resize windows
	 * @param w width of GL widget (do not forget to multiply by devicePixelRatio())
	 * @param h height GL widget (do not forget to multiply by devicePixelRatio())
	 */
	void resize(int w, int h);

	/**
	 * @brief draw the transparent object (can draw only one mesh)
	 * @param param ShaderFlatTransp::Param or ShaderPhongTransp::Param
	 * @param proj projection matrix
	 * @param view modelview matrix
	 * @param draw_func the func/lambda that draw transparent objects
	 */
	// template<typename PARAM, typename TFUNC>
	// void draw(PARAM& param, const QMatrix4x4& proj, const QMatrix4x4& view, const TFUNC& draw_func);

	/**
	 * @brief draw the transparent objects (can draw several meshes)
	 * @param draw_func the func/lambda that draw transparent objects (must bind/unbind ShaderFlatTransp::Param or
	 * ShaderPhongTransp::Param)
	 */
	template <typename TFUNC>
	void draw(const TFUNC& draw_func);

	/**
	 * @brief draw the transparent object with local ShaderFlatTransp::Param (can draw only one mesh)
	 * @param proj projection matrix
	 * @param view modelview matrix
	 * @param draw_func the func/lambda that draw transparent objects
	 */
	template <typename TFUNC>
	void draw_flat(const QMatrix4x4& proj, const QMatrix4x4& view, const TFUNC& draw_func)
	{
		// draw(*param_flat_, proj, view, draw_func);
		draw([&]() -> void {
			param_flat_->bind(proj, view);
			draw_func();
			param_flat_->release();
		});
	}

	/**
	 * @brief draw the transparent object with local ShaderPhongTransp::Param (can draw only one mesh)
	 * @param proj projection matrix
	 * @param view modelview matrix
	 * @param draw_func the func/lambda that draw transparent objects
	 */
	template <typename TFUNC>
	void draw_phong(const QMatrix4x4& proj, const QMatrix4x4& view, const TFUNC& draw_func)
	{
		//		draw(*param_phong_, proj, view, draw_func);
		draw([&]() -> void {
			param_phong_->bind(proj, view);
			draw_func();
			param_phong_->release();
		});
	}

	/**
	 * @brief set the max number of layers (eq drawing passes)
	 * @param nbl
	 */
	inline void set_max_nb_layers(int nbl)
	{
		max_nb_layers_ = nbl;
	}

	inline void set_light_position(const QVector3D& l)
	{
		param_flat_->light_pos_ = l;
		param_phong_->light_pos_ = l;
	}

	inline void set_front_color(const GLColor& rgba)
	{
		param_flat_->front_color_ = rgba;
		param_phong_->front_color_ = rgba;
	}

	inline void set_back_color(const GLColor& rgba)
	{
		param_flat_->back_color_ = rgba;
		param_phong_->back_color_ = rgba;
	}

	inline void set_ambiant_color(const GLColor& rgba)
	{
		param_flat_->ambiant_color_ = rgba;
		param_phong_->ambiant_color_ = rgba;
	}

	inline void set_position_vbo(cgogn::rendering::VBO* vbo_pos)
	{
		param_flat_->set_position_vbo(vbo_pos);
		param_phong_->set_position_vbo(vbo_pos);
	}

	inline void set_normal_vbo(cgogn::rendering::VBO* vbo_norm)
	{
		param_phong_->set_normal_vbo(vbo_norm);
	}

	inline void set_back_face_culling(bool cull)
	{
		param_flat_->bf_culling_ = cull;
		param_phong_->bf_culling_ = cull;
	}

	inline void set_lighted(bool lighted)
	{
		param_flat_->lighted_ = lighted;
	}
};

template <typename TFUNC>
void SurfaceTransparencyDrawer::draw(const TFUNC& draw_func)
{
	if (ogl33_ == nullptr)
		return;

	ShaderFlatTransp* sh_flat = ShaderFlatTransp::get_instance();
	ShaderPhongTransp* sh_phong = ShaderPhongTransp::get_instance();
	ShaderTransparentVolumes* sh_vol = ShaderTransparentVolumes::get_instance();
	QMatrix4x4 fake_mat;

	GLfloat bkColor[4];
	ogl33_->glGetFloatv(GL_COLOR_CLEAR_VALUE, bkColor);

	ogl33_->glEnable(GL_TEXTURE_2D);
	ogl33_->glBindTexture(GL_TEXTURE_2D, depthTexture_);
	ogl33_->glReadBuffer(GL_BACK);
	ogl33_->glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, width_, height_);

	QVector<GLuint> textures = fbo_layer_->textures();
	GLenum buffs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT2};

	sh_flat->bind();
	sh_flat->set_rgba_sampler(0);
	sh_flat->set_depth_sampler(1);
	sh_flat->release();
	sh_phong->bind();
	sh_phong->set_rgba_sampler(0);
	sh_phong->set_depth_sampler(1);
	sh_phong->release();
	sh_vol->bind();
	sh_vol->set_rgba_sampler(0);
	sh_vol->set_depth_sampler(1);
	sh_vol->release();

	fbo_layer_->bind();

	GLenum clear_buff[1] = {GL_COLOR_ATTACHMENT3};
	ogl33_->glDrawBuffers(1, clear_buff);
	ogl33_->glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	ogl33_->glClear(GL_COLOR_BUFFER_BIT);

	ogl33_->glDrawBuffers(1, buffs);
	ogl33_->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	GLenum opaq_buff[1] = {GL_COLOR_ATTACHMENT5};

	for (int p = 0; p < max_nb_layers_; ++p)
	{
		ogl33_->glClear(GL_DEPTH_BUFFER_BIT);

		if (p > 0)
		{
			ogl33_->glDrawBuffers(1, opaq_buff);
			ogl33_->glActiveTexture(GL_TEXTURE0);
			ogl33_->glBindTexture(GL_TEXTURE_2D, depthTexture_);
			param_copy_depth_->depth_texture_sampler_ = 0;
			param_copy_depth_->bind(fake_mat, fake_mat);
			ogl33_->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
			param_copy_depth_->release();
		}

		ogl33_->glDrawBuffers(2, buffs);
		sh_flat->bind();
		sh_flat->set_layer(p);
		sh_flat->release();
		sh_phong->bind();
		sh_phong->set_layer(p);
		sh_phong->release();
		sh_vol->bind();
		sh_vol->set_layer(p);
		sh_vol->release();

		ogl33_->glActiveTexture(GL_TEXTURE0);
		ogl33_->glBindTexture(GL_TEXTURE_2D, textures[3]);
		ogl33_->glActiveTexture(GL_TEXTURE1);
		ogl33_->glBindTexture(GL_TEXTURE_2D, textures[1]);

		ogl33_->glBeginQuery(GL_SAMPLES_PASSED, oq_transp_);
		draw_func();
		ogl33_->glEndQuery(GL_SAMPLES_PASSED);

		GLuint nb_samples;
		ogl33_->glGetQueryObjectuiv(oq_transp_, GL_QUERY_RESULT, &nb_samples);

		if (nb_samples == 0) // finished ?
			p = max_nb_layers_;
		else
		{
			ogl33_->glReadBuffer(GL_COLOR_ATTACHMENT2);
			ogl33_->glBindTexture(GL_TEXTURE_2D, textures[1]);
			ogl33_->glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 0, 0, width_, height_, 0);

			if (p == 0)
			{
				ogl33_->glBindTexture(GL_TEXTURE_2D, textures[4]);
				ogl33_->glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 0, 0, width_, height_, 0);
			}

			ogl33_->glReadBuffer(GL_COLOR_ATTACHMENT0);
			ogl33_->glBindTexture(GL_TEXTURE_2D, textures[3]);
			ogl33_->glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 0, 0, width_, height_, 0);
		}
	}

	fbo_layer_->release();

	// real draw with blending with opaque object

	param_trq_->rgba_texture_sampler_ = 0;
	param_trq_->depth_texture_sampler_ = 1;

	ogl33_->glActiveTexture(GL_TEXTURE0);
	ogl33_->glBindTexture(GL_TEXTURE_2D, textures[3]);

	ogl33_->glActiveTexture(GL_TEXTURE1);
	ogl33_->glBindTexture(GL_TEXTURE_2D, textures[4]);

	ogl33_->glEnable(GL_BLEND);
	ogl33_->glBlendFunc(GL_ONE, GL_SRC_ALPHA);
	param_trq_->bind(fake_mat, fake_mat);
	ogl33_->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	param_trq_->release();
	ogl33_->glDisable(GL_BLEND);

	ogl33_->glClearColor(bkColor[0], bkColor[1], bkColor[2], bkColor[3]);
}

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_FLAT_TR_DR_H_
