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

#include <cgogn/rendering/wall_paper.h>

#include <iostream>

namespace cgogn
{

namespace rendering
{

void WallPaper::init(const GLImage& img)
{
	if (vbo_pos_ == nullptr)
	{
		vbo_pos_ = std::make_unique<VBO>(3);
		vbo_pos_->allocate(4, 3);
	}
	if (vbo_tc_ == nullptr)
	{
		vbo_tc_ = std::make_unique<VBO>(2);
		vbo_tc_->allocate(4, 2);
		float32* ptr_tc = vbo_tc_->lock_pointer();
		*ptr_tc++ = 0.0f;
		*ptr_tc++ = 0.0f;
		*ptr_tc++ = 1.0f;
		*ptr_tc++ = 0.0f;
		*ptr_tc++ = 1.0f;
		*ptr_tc++ = 1.0f;
		*ptr_tc++ = 0.0f;
		*ptr_tc++ = 1.0f;
		vbo_tc_->release_pointer();
	}
	texture_ = std::make_unique<Texture2D>();
	texture_->load(img);
	set_full_screen(false);
}

WallPaper::WallPaper(const GLImage& img) : vbo_pos_(nullptr), vbo_tc_(nullptr), texture_(nullptr)
{
	init(img);
}

WallPaper::WallPaper(const GLColor& col) : vbo_pos_(nullptr), vbo_tc_(nullptr), texture_(nullptr)
{
	GLImage img(1, 1, 3);
	img.set_pixel(0, 0, col);
	init(img);
}

WallPaper::WallPaper(const GLColor& col_tl, const GLColor& col_tr, const GLColor& col_bl, const GLColor& col_br)
	: vbo_pos_(nullptr), vbo_tc_(nullptr), texture_(nullptr)
{
	GLImage img(2, 2, 3);
	img.set_pixel(0, 1, col_bl);
	img.set_pixel(1, 1, col_br);
	img.set_pixel(1, 0, col_tr);
	img.set_pixel(0, 0, col_tl);
	init(img);
}

WallPaper::~WallPaper()
{
}

void WallPaper::set_full_screen(bool front)
{
	float32 depth = 0.9999999f;
	if (front)
		depth = 0.0f;

	float32* ptr_pos = vbo_pos_->lock_pointer();
	*ptr_pos++ = -1.0f;
	*ptr_pos++ = 1.0f;
	*ptr_pos++ = depth;
	*ptr_pos++ = 1.0f;
	*ptr_pos++ = 1.0f;
	*ptr_pos++ = depth;
	*ptr_pos++ = 1.0f;
	*ptr_pos++ = -1.0f;
	*ptr_pos++ = depth;
	*ptr_pos++ = -1.0f;
	*ptr_pos++ = -1.0f;
	*ptr_pos++ = depth;
	vbo_pos_->release_pointer();
}

void WallPaper::set_local_position(uint32 win_w, uint32 win_h, uint32 x, uint32 y, uint32 w, uint32 h, bool front)
{
	float32 depth = 0.0f;
	if (!front)
		depth = 0.9999999f;

	float32 xmin = -1.0f + float32(2 * x) / float32(win_w);
	float32 xmax = xmin + float32(2 * w) / float32(win_w);

	float32 ymin = 1.0f - float32(2 * y) / float32(win_h);
	float32 ymax = ymin - float32(2 * h) / float32(win_h);

	float32* ptr_pos = vbo_pos_->lock_pointer();
	*ptr_pos++ = xmin;
	*ptr_pos++ = ymin;
	*ptr_pos++ = depth;
	*ptr_pos++ = xmax;
	*ptr_pos++ = ymin;
	*ptr_pos++ = depth;
	*ptr_pos++ = xmax;
	*ptr_pos++ = ymax;
	*ptr_pos++ = depth;
	*ptr_pos++ = xmin;
	*ptr_pos++ = ymax;
	*ptr_pos++ = depth;
	vbo_pos_->release_pointer();
}

void WallPaper::set_local_position(float x, float y, float w, float h, bool front)
{
	float32 depth = 0.0f;
	if (!front)
		depth = 0.9999999f;

	float32 xmin = -1.0f + 2 * x;
	float32 xmax = xmin + 2 * w;

	float32 ymin = 1.0f - 2 * y;
	float32 ymax = ymin - 2 * h;

	float32* ptr_pos = vbo_pos_->lock_pointer();
	*ptr_pos++ = xmin;
	*ptr_pos++ = ymin;
	*ptr_pos++ = depth;
	*ptr_pos++ = xmax;
	*ptr_pos++ = ymin;
	*ptr_pos++ = depth;
	*ptr_pos++ = xmax;
	*ptr_pos++ = ymax;
	*ptr_pos++ = depth;
	*ptr_pos++ = xmin;
	*ptr_pos++ = ymax;
	*ptr_pos++ = depth;
	vbo_pos_->release_pointer();
}

WallPaper::Renderer::Renderer(WallPaper* wp) : wall_paper_data_(wp)
{
	param_texture_ = ShaderTexture::generate_param();
	param_texture_->set_vbos(wp->vbo_pos_.get(), wp->vbo_tc_.get());
	param_texture_->texture_ = wp->texture_.get();
}

WallPaper::Renderer::~Renderer()
{
}

void WallPaper::Renderer::draw()
{
	param_texture_->bind(GLMat4::Identity(), GLMat4::Identity());
	glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
	param_texture_->release();
}

} // namespace rendering

} // namespace cgogn
