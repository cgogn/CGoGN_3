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

#ifndef CGOGN_RENDERING_FBO_H_
#define CGOGN_RENDERING_FBO_H_

#include <GL/gl3w.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_fullscreen_texture.h>
#include <cgogn/rendering/texture.h>

#include <vector>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT FBO
{
	GLint prev_id_;
	GLint prev_viewport[4];

public:
	FBO(const std::vector<Texture2D*>& textures, bool add_depth, FBO* from);

	inline void bind()
	{
		glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &prev_id_);
		glGetIntegerv(GL_VIEWPORT, prev_viewport);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, id_);
		glViewport(0, 0, tex_[0]->width(), tex_[0]->height());
		std::cout << "FBO BIND "<< id_ <<" VP 0, 0, "<< tex_[0]->width() <<","<< tex_[0]->height() <<std::endl;
	}

	inline void release()
	{
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, prev_id_);
		glViewport(prev_viewport[0], prev_viewport[1], prev_viewport[2], prev_viewport[3]);GL_ASSERT()

		std::cout << "BACK TO FBO "<< prev_id_ <<" VP "<<prev_viewport[0] <<","<< prev_viewport[1] <<","<<  prev_viewport[2]  <<","<<  prev_viewport[3]<<std::endl;
	}


	inline void bind_read()
	{
		glBindFramebuffer(GL_READ_FRAMEBUFFER, id_);GL_ASSERT()
	}

	inline void release_read()
	{
		glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);GL_ASSERT()
	}

	void resize(int w, int h);

	inline Texture2D* texture(std::size_t i)
	{
		return tex_[i];
	}

	inline std::size_t nb_textures()
	{
		return tex_.size();
	}

	inline GLint width() const
	{
		return tex_.front()->width();
	}
	inline GLint height() const
	{
		return tex_.front()->height();
	}

protected:
	GLuint id_;
	GLuint depth_render_buffer_;
	std::vector<Texture2D*> tex_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_FBO_H_
