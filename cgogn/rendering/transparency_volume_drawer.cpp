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

#define CGOGN_RENDERING_TRANSP_VOLUME_RENDER_CPP_

#include <cgogn/rendering/transparency_volume_drawer.h>

#include <QColor>
#include <QImage>
#include <QOpenGLFunctions>
#include <iostream>

namespace cgogn
{

namespace rendering
{

VolumeTransparencyDrawer::VolumeTransparencyDrawer() : vbo_pos_(nullptr), face_color_(0, 150, 0), shrink_v_(0.6f)
{
	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
}

VolumeTransparencyDrawer::Renderer::~Renderer()
{
	param_transp_vol_.reset();
}

VolumeTransparencyDrawer::Renderer::Renderer(VolumeTransparencyDrawer* vr)
	: param_transp_vol_(nullptr), volume_drawer_data_(vr)
{
	param_transp_vol_ = ShaderTransparentVolumes::generate_param();
	param_transp_vol_->set_position_vbo(vr->vbo_pos_.get());
	param_transp_vol_->explode_factor_ = vr->shrink_v_;
	param_transp_vol_->color_ = vr->face_color_;
}

void VolumeTransparencyDrawer::Renderer::draw_faces(const QMatrix4x4& projection, const QMatrix4x4& modelview)
{
	QOpenGLFunctions_3_3_Core* ogl33 = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();
	param_transp_vol_->bind(projection, modelview);
	ogl33->glDrawArrays(GL_LINES_ADJACENCY, 0, volume_drawer_data_->vbo_pos_->size());
	param_transp_vol_->release();
}

void VolumeTransparencyDrawer::Renderer::set_explode_volume(float32 x)
{
	if (param_transp_vol_)
		param_transp_vol_->explode_factor_ = x;
}

void VolumeTransparencyDrawer::Renderer::set_color(const QColor& rgb)
{
	if (param_transp_vol_)
		param_transp_vol_->color_ = rgb;
}

void VolumeTransparencyDrawer::Renderer::set_clipping_plane(const QVector4D& pl)
{
	if (param_transp_vol_)
		param_transp_vol_->plane_clip_ = pl;
}

void VolumeTransparencyDrawer::Renderer::set_clipping_plane2(const QVector4D& pl)
{
	if (param_transp_vol_)
		param_transp_vol_->plane_clip2_ = pl;
}

void VolumeTransparencyDrawer::Renderer::set_thick_clipping_plane(const QVector4D& p, float32 th)
{
	QVector4D p1 = p;
	p1[3] -= th / 2.0f;
	set_clipping_plane(p1);

	QVector4D p2 = -p;
	p2[3] -= th / 2.0f;
	set_clipping_plane2(p2);
}

void VolumeTransparencyDrawer::Renderer::set_back_face_culling(bool cull)
{
	param_transp_vol_->bf_culling_ = cull;
}

void VolumeTransparencyDrawer::Renderer::set_lighted(bool lighted)
{
	param_transp_vol_->lighted_ = lighted;
}

} // namespace rendering

} // namespace cgogn
