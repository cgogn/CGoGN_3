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

#define CGOGN_RENDERING_VOLUME_RENDER_CPP_

#include <cgogn/rendering/volume_drawer.h>

namespace cgogn
{

namespace rendering
{

VolumeDrawerGen::VolumeDrawerGen(bool with_color_per_face) :
	vbo_pos_(nullptr),
	vbo_col_(nullptr),
	face_color_(0,0.8f,0,1),
	vbo_pos2_(nullptr),
	edge_color_(0,0,0,1),
	shrink_v_(0.6f)
{
	vbo_pos_ = cgogn::make_unique<VBO>(3);
	vbo_pos2_ = cgogn::make_unique<VBO>(3);

	if (with_color_per_face)
		vbo_col_ = cgogn::make_unique<VBO>(3);
}

VolumeDrawerGen::~VolumeDrawerGen()
{}

VolumeDrawerGen::Renderer::Renderer(VolumeDrawerGen* vr) :
	param_expl_vol_(nullptr),
	param_expl_vol_col_(nullptr),
	param_expl_vol_line_(nullptr),
	volume_drawer_data_(vr)
{
	if (vr->vbo_col_)
	{
		param_expl_vol_col_ = ShaderExplodeVolumesColor::generate_param();
		param_expl_vol_col_->explode_vol_ = vr->shrink_v_;
		param_expl_vol_col_->set_vbos(vr->vbo_pos_.get(), vr->vbo_col_.get());
	}
	else
	{
		param_expl_vol_ = ShaderExplodeVolumes::generate_param();
		param_expl_vol_->set_vbos(vr->vbo_pos_.get());
		param_expl_vol_->explode_vol_ = vr->shrink_v_;
		param_expl_vol_->color_ = vr->face_color_;
	}

	if (vr->vbo_pos2_)
	{
		param_expl_vol_line_ = ShaderExplodeVolumesLine::generate_param();
		param_expl_vol_line_->set_vbos(vr->vbo_pos2_.get());
		param_expl_vol_line_->explode_factor_ = vr->shrink_v_;
		param_expl_vol_line_->color_ = vr->edge_color_;
	}
}

VolumeDrawerGen::Renderer::~Renderer()
{}

void VolumeDrawerGen::Renderer::draw_faces(const GLMat4& projection, const GLMat4& modelview)
{
	if (param_expl_vol_col_)
	{
		param_expl_vol_col_->bind(projection, modelview);
		glDrawArrays(GL_LINES_ADJACENCY, 0, volume_drawer_data_->vbo_pos_->size());
		param_expl_vol_col_->release();
	}
	else
	{ 
		param_expl_vol_->bind(projection, modelview);
		glDrawArrays(GL_LINES_ADJACENCY, 0, volume_drawer_data_->vbo_pos_->size());
		param_expl_vol_->release();
	}
}

void VolumeDrawerGen::Renderer::draw_edges(const GLMat4& projection, const GLMat4& modelview)
{
	param_expl_vol_line_->bind(projection,modelview);
	glDrawArrays(GL_TRIANGLES, 0, volume_drawer_data_->vbo_pos2_->size());
	param_expl_vol_line_->release();
}

void VolumeDrawerGen::Renderer::set_explode_volume(float32 x)
{
	if (param_expl_vol_)
		param_expl_vol_->explode_vol_ = x;
	if (param_expl_vol_col_)
		param_expl_vol_col_->explode_vol_ = x;
	if (param_expl_vol_line_)
		param_expl_vol_line_->explode_factor_ = x;
}

void VolumeDrawerGen::Renderer::set_face_color(const GLColor& rgba)
{
	if (param_expl_vol_)
		param_expl_vol_->color_ = rgba;
}

void VolumeDrawerGen::Renderer::set_edge_color(const GLColor& rgba)
{
	if (param_expl_vol_line_)
		param_expl_vol_line_->color_=rgba;
}

void VolumeDrawerGen::Renderer::set_clipping_plane(const GLVec4& pl)
{
	if (param_expl_vol_)
		param_expl_vol_->plane_clip_ = pl;
	if (param_expl_vol_col_)
		param_expl_vol_col_->plane_clip_ = pl;
	if (param_expl_vol_line_)
		param_expl_vol_line_->plane_clip_ = pl;
}

void VolumeDrawerGen::Renderer::set_clipping_plane2(const GLVec4& pl)
{
	if (param_expl_vol_)
		param_expl_vol_->plane_clip2_ = pl;
	if (param_expl_vol_col_)
		param_expl_vol_col_->plane_clip2_ = pl;
	if (param_expl_vol_line_)
		param_expl_vol_line_->plane_clip2_ = pl;
}

void VolumeDrawerGen::Renderer::set_thick_clipping_plane(const GLVec4& p, float32 th)
{
	GLVec4 p1 = p;
	p1[3] -= th/2.0f;
	set_clipping_plane(p1);

	GLVec4 p2 = -p;
	p2[3] -= th/2.0f;
	set_clipping_plane2(p2);
}

VolumeDrawerTpl<false>::~VolumeDrawerTpl()
{}

VolumeDrawerTpl<true>::~VolumeDrawerTpl()
{}

} // namespace rendering

} // namespace cgogn
