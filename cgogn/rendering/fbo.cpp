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

#include <cgogn/rendering/fbo.h>

namespace cgogn
{

namespace rendering
{

FBO::FBO(const std::vector<Texture2D*>& textures, bool add_depth, FBO* from)
{
	glGenFramebuffers(1, &id_);
	glBindFramebuffer(GL_FRAMEBUFFER, id_);
	GLenum att = GL_COLOR_ATTACHMENT0;
	for (auto* t : textures)
	{
		tex_.clear();
		tex_.push_back(t);
		glFramebufferTexture2D(GL_FRAMEBUFFER, att++, GL_TEXTURE_2D, t->id(), 0);
	}

	//	if (add_depth)
	//	{
	//		if (from)
	//		{
	//			depth_render_buffer_ = from->depth_render_buffer_;
	//			glBindRenderbuffer(GL_RENDERBUFFER, depth_render_buffer_);
	//			glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_render_buffer_);
	//			glBindRenderbuffer(GL_RENDERBUFFER, 0);
	//		}
	//		else
	//		{
	//			glGenRenderbuffers(1, &depth_render_buffer_);
	//			glBindRenderbuffer(GL_RENDERBUFFER, depth_render_buffer_);
	//			glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, tex_[0]->width(), tex_[0]->height());
	//			glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_render_buffer_);
	//			glBindRenderbuffer(GL_RENDERBUFFER, 0);
	//		}
	//	}

	if (add_depth)
	{
		if (from)
			depth_tex_ = from->depth_tex_;
		else
		{
			depth_tex_ = new Texture2D({
				{GL_TEXTURE_MIN_FILTER, GL_NEAREST},
				{GL_TEXTURE_MAG_FILTER, GL_NEAREST},
				{GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE},
				{GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE},
			});
			depth_tex_->allocate(0, 0, GL_DEPTH_COMPONENT32F, GL_DEPTH_COMPONENT);
		}
		glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_tex_->id(), 0);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	}
	else
		depth_tex_ = nullptr;
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void FBO::resize(int w, int h)
{
	for (auto* t : tex_)
		t->resize(w, h);

	if (depth_tex_ != nullptr)
		depth_tex_->resize(w, h);
}

} // namespace rendering

} // namespace cgogn
