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

#define CGOGN_RENDERING_TEXT_RENDER_CPP_

#include <cgogn/rendering/text_drawer.h>
#include <iostream>

namespace cgogn
{

namespace rendering
{

TextDrawer::End TextDrawer::end = TextDrawer::End();

Texture2D* TextDrawer::texture_ = nullptr;


TextDrawer::TextDrawer() :
	vbo_pos_(nullptr),
	vbo_char_(nullptr),
	vbo_colsz_(nullptr)
{
//	Q_INIT_RESOURCE(fonte);
	vbo_pos_ = cgogn::make_unique<VBO>(4);
	vbo_char_ = cgogn::make_unique<VBO>(1);
	vbo_colsz_ = cgogn::make_unique<VBO>(4);
	GLImage img("fonte4064.png");
	if (texture_ == nullptr)
	{
		texture_ = new Texture2D();
		texture_->load(img);
	}

}

TextDrawer::~TextDrawer()
{}


TextDrawer& TextDrawer::operator << (const GLColor& col)
{
	current_color_ = col;
	return *this;
}

TextDrawer& TextDrawer::operator << (float32 sz)
{
	current_size_ = sz;
	return *this;
}

TextDrawer& TextDrawer::operator << (const std::string& str)
{
	if (next_pos_)
	{
		strings_.push_back(str);
		positions_.push_back(current_pos_);
		colors_.push_back(current_color_);
		sizes_.push_back(current_size_);
		next_pos_ = false;
	}
	else
		strings_.back() += str;
	
	return *this;
}

TextDrawer& TextDrawer::operator << (TextDrawer::End)
{
	std::size_t nb = 0;
	for (const auto& s : strings_)
		nb += s.size();

	std::vector<Vec4f> pos;
	pos.reserve(nb);
	std::vector<float32> chars;
	chars.reserve(nb);
	std::vector<Vec4f> colsz;
	colsz.reserve(nb);

	auto it = positions_.begin();
	auto jt = colors_.begin();
	auto kt = sizes_.begin();
	for (const auto& s : strings_)
	{
		Vec4f p{ (*it)[0], (*it)[1], (*it)[2], 0 };
		Vec4f cs{ jt->x(), jt->y(), jt->z(), *kt++ };
		it++;
		jt++;
		for (auto c : s)
		{
			chars.push_back(float32(c - ' ') / 95.0f);
			pos.push_back(p);
			p[3] += 1;
			colsz.push_back(cs);
		}

	}

	vbo_pos_->allocate(nb, 4);
	vbo_pos_->bind();
	vbo_pos_->copy_data(0, nb * sizeof(Vec4f), pos.data());
	vbo_pos_->release();

	vbo_char_->allocate(nb, 1);
	vbo_char_->bind();
	vbo_char_->copy_data(0, nb * sizeof(float32), chars.data());
	vbo_char_->release();

	vbo_colsz_->allocate(nb, 4);
	vbo_colsz_->bind();
	vbo_colsz_->copy_data(0, nb * sizeof(Vec4f), colsz.data());
	vbo_colsz_->release();

	positions_.clear();
	strings_.clear();
	colors_.clear();
	sizes_.clear();
	positions_.shrink_to_fit();
	strings_.shrink_to_fit();
	colors_.shrink_to_fit();
	sizes_.shrink_to_fit();
	next_pos_ = false;

	return *this;
}



void TextDrawer::update_text(std::size_t pos, const std::string& str)
{
	float32* chars = vbo_char_->lock_pointer() + pos;

	for (auto c : str)
		*chars++ = float32(c - ' ') / 95.0f;

	vbo_char_->release_pointer();
}

void TextDrawer::scale_text(float sc)
{
	float32* sz = vbo_colsz_->lock_pointer() + 3;

	for (uint32 i = 0; i != vbo_colsz_->size(); ++i)
	{
		*sz *= sc;
		sz += 4;
	}

	vbo_colsz_->release_pointer();
}




TextDrawer::Renderer::Renderer(TextDrawer* tr) :
	param_text_(nullptr),
	text_drawer_data_(tr)
{
	param_text_ = ShaderText::generate_param();
	param_text_->texture_ = text_drawer_data_->texture_;
	param_text_->italic_ = 0;
	param_text_->set_vbos(text_drawer_data_->vbo_pos_.get(), text_drawer_data_->vbo_char_.get(), text_drawer_data_->vbo_colsz_.get());
}

TextDrawer::Renderer::~Renderer()
{}


void TextDrawer::Renderer::draw(const GLMat4& projection, const GLMat4& modelview)
{
	param_text_->bind(projection, modelview);
	glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, text_drawer_data_->vbo_pos_->size());
	param_text_->release();
}


} // namespace rendering

} // namespace cgogn
