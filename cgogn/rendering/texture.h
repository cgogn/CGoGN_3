
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

#ifndef CGOGN_RENDERING_TEXTURE_H_
#define CGOGN_RENDERING_TEXTURE_H_

#include <GL/gl3w.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/types.h>

#include <string>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT Texture2D
{
public:
	inline Texture2D() : internal_(0), external_(0), data_type_(0), width_(0), height_(0), depth_(false)
	{
		glGenTextures(1, &id_);
		glBindTexture(GL_TEXTURE_2D, id_);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		depth_ = false;
	}

	inline ~Texture2D()
	{
		glDeleteTextures(1, &id_);
	}

	inline GLuint id()
	{
		return id_;
	}

	inline void alloc(GLsizei w, GLsizei h, GLint internal, GLenum external, const uint8* ptr = nullptr,
					  GLenum data_type = GL_UNSIGNED_BYTE)
	{
		internal_ = internal;
		external_ = external;
		width_ = w;
		height_ = h;
		data_type_ = data_type;
		depth_ = ((internal == GL_DEPTH_COMPONENT32F) || (internal == GL_DEPTH_COMPONENT24));
		if (w * h > 0)
		{
			bind();
			glTexImage2D(GL_TEXTURE_2D, 0, internal, w, h, 0, external, data_type, ptr);
			release();
		}
	}

	inline void load(const GLImage& img)
	{
		switch (img.depth())
		{
		case 1:
			alloc(img.width(), img.height(), GL_R8, GL_RED, img.data());
			break;
		case 3:
			alloc(img.width(), img.height(), GL_RGB8, GL_RGB, img.data());
			break;
		case 4:
			alloc(img.width(), img.height(), GL_RGBA8, GL_RGBA, img.data());
			break;
		}
		
	}

	inline void resize(GLsizei w, GLsizei h)
	{
		bind();
		glTexImage2D(GL_TEXTURE_2D, 0, internal_, w, h, 0, external_, data_type_, nullptr);
		width_ = w;
		height_ = h;
		release();
		
	}

	inline GLsizei width() const
	{
		return width_;
	}

	inline GLsizei height() const
	{
		return height_;
	}

	inline void bind()
	{
		glBindTexture(GL_TEXTURE_2D, id_);
		
	}

	inline GLint bind(GLint unit)
	{
		glActiveTexture(GL_TEXTURE0 + GLuint(unit));
		glBindTexture(GL_TEXTURE_2D, id_);
		
		return unit;
	}

	inline static void release()
	{
		glBindTexture(GL_TEXTURE_2D, 0);
		
	}

protected:
	GLuint id_;
	GLint internal_;
	GLenum external_;
	GLenum data_type_;
	GLsizei width_;
	GLsizei height_;
	bool depth_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_TEXTURE_H_
