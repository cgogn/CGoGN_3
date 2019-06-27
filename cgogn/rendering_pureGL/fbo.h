
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

#ifndef CGOGN_RENDERING_FBO_H_
#define CGOGN_RENDERING_FBO_H_

#include <GL/gl3w.h>
#include <vector>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>
#include <cgogn/rendering_pureGL/texture.h>
#include <cgogn/rendering_pureGL/shaders/shader_fullscreen_texture.h>


namespace cgogn
{

namespace rendering_pgl
{

class CGOGN_RENDERING_PUREGL_EXPORT FBO
{
	GLint initial_viewport_[4];
public:
	FBO(const std::vector<Texture2D*>& textures, bool add_depth, FBO* from );

	inline void bind()
	{
		glGetIntegerv(GL_VIEWPORT, initial_viewport_);
		glBindFramebuffer(GL_FRAMEBUFFER, id_);
		glViewport(0,0,tex_[0]->width(),tex_[0]->height());
	}

	inline void release()
	{
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(initial_viewport_[0],initial_viewport_[1],initial_viewport_[2],initial_viewport_[3]);
	}

	void resize(int w, int h);

	inline Texture2D* texture(std::size_t i) { return tex_[i]; }

	inline std::size_t nb_textures() { return tex_.size(); }

	inline GLint width() const { return tex_.front()->width(); }
	inline GLint height() const { return tex_.front()->height(); }

protected:
	GLuint id_;
	GLuint depth_render_buffer_;
	std::vector<Texture2D*> tex_;
};

}
}
#endif
