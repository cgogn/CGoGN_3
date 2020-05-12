
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

#ifndef CGOGN_RENDERING_VBO_H_
#define CGOGN_RENDERING_VBO_H_

#include <GL/gl3w.h>

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/types.h>

#include <cgogn/core/utils/numerics.h>

#include <iostream>
#include <string>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT VBO
{
protected:
	GLuint id_;
	GLuint id_texture_buffer_;
	std::size_t nb_vectors_;
	int32 vector_dimension_;
	std::string name_;

public:
	inline VBO(int32 vec_dim = 3) : id_texture_buffer_(0), nb_vectors_(0), vector_dimension_(vec_dim)
	{
		glGenBuffers(1, &id_);
	}

	inline ~VBO()
	{
		glDeleteBuffers(1, &id_);
		if (id_texture_buffer_ != 0)
			glDeleteTextures(1, &id_texture_buffer_);
	}

	/**
	 * @brief dimension of vectors stored in buffer
	 */
	inline int32 vector_dimension() const
	{
		return vector_dimension_;
	}

	inline uint32 size() const
	{
		return uint32(nb_vectors_);
	}

	inline GLuint id() const
	{
		return id_;
	}

	inline GLuint id_texture_buffer() const
	{
		return id_texture_buffer_;
	}

	inline void set_name(const std::string& name)
	{
		name_ = name;
		gl_debug_name(GL_BUFFER, id_, "VBO_" + name_);
	}

	inline const std::string& name() const
	{
		return name_;
	}

	inline void bind()
	{
		glBindBuffer(GL_ARRAY_BUFFER, id_);
	}

	inline void release()
	{
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	inline void bind_texture_buffer(GLint unit)
	{
		glActiveTexture(GL_TEXTURE0 + unit);
		if (id_texture_buffer_ == 0)
		{
			static GLenum internals[] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
			glGenTextures(1, &id_texture_buffer_);
			glBindTexture(GL_TEXTURE_BUFFER, id_texture_buffer_);
			glTexBuffer(GL_TEXTURE_BUFFER, internals[vector_dimension_ - 1], id_);
		}
		else
			glBindTexture(GL_TEXTURE_BUFFER, id_texture_buffer_);
	}

	inline static void release_texture_buffer(GLint unit)
	{
		glActiveTexture(GL_TEXTURE0 + unit);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
	}

	/**
	 * @brief allocate VBO memory
	 * @param nb_vectors number of vectors
	 * @param vector_dimension number of component of each vector
	 */
	inline void allocate(std::size_t nb_vectors, int32 vector_dimension)
	{
		std::size_t total = nb_vectors * uint32(vector_dimension);
		if (total != nb_vectors_ * uint64(vector_dimension)) // only allocate when > ?
		{
			glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(total * 4), nullptr, GL_DYNAMIC_DRAW);
			nb_vectors_ = nb_vectors;
			vector_dimension_ = vector_dimension;
		}
	}

	/**
	 * @brief get and lock pointer on buffer memory (muste be bind)
	 * @return the pointer
	 */
	inline float32* lock_pointer()
	{
		float32* ptr = reinterpret_cast<float32*>(glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE));
		return ptr;
	}

	inline float32* lock_pointer_read()
	{
		float32* ptr = reinterpret_cast<float32*>(glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY));
		return ptr;
	}

	/**
	 * @brief release_pointer
	 */
	inline void release_pointer()
	{
		glUnmapBuffer(GL_ARRAY_BUFFER);
	}

	/**
	 * @brief copy data
	 * @param offset offset in bytes in the bufffer
	 * @param nb number of bytes to copy
	 * @param src source pointer
	 */
	inline void copy_data(uint32 offset, std::size_t nb, const void* src)
	{
		glBufferSubData(GL_ARRAY_BUFFER, offset, GLsizeiptr(nb), src);
	}

	inline void associate(GLuint attrib, int32 stride = 0, uint32 first = 0)
	{
		bind();
		glEnableVertexAttribArray(attrib);
		glVertexAttribPointer(attrib, vector_dimension(), GL_FLOAT, GL_FALSE, stride * vector_dimension() * 4,
							  reinterpret_cast<GLvoid*>(first * uint64(vector_dimension()) * 4u));
		release();
	}
};

inline std::ostream& operator<<(std::ostream& out, VBO& vbo)
{
	const int NB = 5;

	auto vd = vbo.vector_dimension();

	std::cout << "DEBUG VBO " << vbo.id() << " : " << vbo.name();
	std::cout << " ( " << vbo.size() << " vec of dim " << vd << " ) " << std::endl;
	vbo.bind();
	float* ff = vbo.lock_pointer_read();
	float* f = ff;
	for (int i = 0; i < NB; ++i)
	{
		for (int j = 0; j < vd; ++j)
			std::cout << *f++ << ", ";
		std::cout << std::endl;
	}

	std::cout << " . . . . " << std::endl;

	f = ff + vbo.size() - int64(NB) * vd;

	for (int i = 0; i < NB; ++i)
	{
		for (int j = 0; j < vd; ++j)
			std::cout << *f++ << ", ";
		std::cout << std::endl;
	}

	std::cout << "----------end vbo debug-----------" << std::endl;
	vbo.release_pointer();
	vbo.release();
	return out;
}

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_VBO_H_
