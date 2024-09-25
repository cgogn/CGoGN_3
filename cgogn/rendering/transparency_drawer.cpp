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

#include <cgogn/rendering/transparency_drawer.h>

namespace cgogn
{

namespace rendering
{

SurfaceTransparencyDrawer::~SurfaceTransparencyDrawer()
{
	param_flat_.reset();
	param_trq_.reset();
	fbo_layer_.reset();
	if (ogl33_)
		ogl33_->glDeleteQueries(1, &oq_transp_);
}

SurfaceTransparencyDrawer::SurfaceTransparencyDrawer()
	: max_nb_layers_(8), param_flat_(nullptr), param_trq_(nullptr), fbo_layer_(nullptr), oq_transp_(0u),
	  ogl33_(nullptr), depthTexture_(0)
{
	param_flat_ = cgogn::rendering::ShaderFlatTransp::generate_param();
	param_flat_->front_color_ = QColor(0, 250, 0, 120);
	param_flat_->back_color_ = QColor(0, 0, 250, 120);
	param_flat_->ambiant_color_ = QColor(0, 0, 0, 0);

	param_phong_ = cgogn::rendering::ShaderPhongTransp::generate_param();
	param_phong_->front_color_ = QColor(0, 250, 0, 120);
	param_phong_->back_color_ = QColor(0, 0, 250, 120);
	param_phong_->ambiant_color_ = QColor(0, 0, 0, 0);
	param_phong_->specular_color_ = QColor(255, 255, 255, 0);
	param_phong_->specular_coef_ = 100.0f;

	param_trq_ = cgogn::rendering::ShaderTranspQuad::generate_param();

	param_copy_depth_ = ShaderCopyDepth::generate_param();
}

void SurfaceTransparencyDrawer::resize(int w, int h)
{
	QOpenGLFunctions_3_3_Core* ogl33 = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();

	width_ = w;
	height_ = h;
	ogl33_ = ogl33;

	fbo_layer_ = cgogn::make_unique<QOpenGLFramebufferObject>(width_, height_, QOpenGLFramebufferObject::Depth,
															  GL_TEXTURE_2D, /*GL_RGBA8*/ GL_RGBA32F);
	fbo_layer_->addColorAttachment(width_, height_, GL_R32F);
	fbo_layer_->addColorAttachment(width_, height_, GL_R32F);
	fbo_layer_->addColorAttachment(width_, height_);
	fbo_layer_->addColorAttachment(width_, height_, GL_R32F); // first depth
	fbo_layer_->addColorAttachment(width_, height_);

	if (!ogl33->glIsQuery(oq_transp_))
		ogl33->glGenQueries(1, &oq_transp_);

	if (ogl33->glIsTexture(depthTexture_))
		ogl33->glDeleteTextures(1, &depthTexture_);

	ogl33_->glGenTextures(1, &depthTexture_);
	ogl33_->glBindTexture(GL_TEXTURE_2D, depthTexture_);
	ogl33_->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	ogl33_->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	ogl33_->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	ogl33_->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	ogl33_->glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width_, height_, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
}

} // namespace rendering

} // namespace cgogn
