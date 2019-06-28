
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

#ifndef CGOGN_RENDERING_VBO_H_
#define CGOGN_RENDERING_VBO_H_

#include <GL/gl3w.h>

#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>
#include <cgogn/core/utils/numerics.h>

#include <string>

namespace cgogn
{

namespace rendering_pgl
{

class CGOGN_RENDERING_PUREGL_EXPORT VBO
{
protected:

	GLuint id_;
	std::size_t nb_vectors_;
	int32 vector_dimension_;
	std::string name_;

public:

	inline VBO(int32 vec_dim = 3) :
		nb_vectors_(0),
		vector_dimension_(vec_dim)
	{
		glGenBuffers(1, &id_);
	}

	inline ~VBO()
	{
		glDeleteBuffers(1,&id_);
		id_ = 0;
	}

	inline void set_name(const std::string& name)
	{
		name_ = name;
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

	/**
	 * @brief allocate VBO memory
	 * @param nb_vectors number of vectors
	 * @param vector_dimension_ number of component of each vector
	 */
	inline void allocate(std::size_t nb_vectors, int32 vector_dim)
	{
		std::size_t total = nb_vectors * uint32(vector_dim);
//		if (total != nb_vectors_ * uint64(vector_dim)) // only allocate when > ?
		{
			glBindBuffer(GL_ARRAY_BUFFER, id_);
			glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(total*4), nullptr, GL_STATIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
		}
		nb_vectors_ = nb_vectors;
		vector_dimension_ = vector_dim;
	}

	/**
	 * @brief get and lock pointer on buffer memory
	 * @return  the pointer
	 */
	inline float32* lock_pointer()
	{
		this->bind();
		float32* ptr = reinterpret_cast<float32*>(glMapBuffer(GL_ARRAY_BUFFER,GL_READ_WRITE));
		this->release();
		return ptr;
	}

	/**
	 * @brief release_pointer
	 */
	inline void release_pointer()
	{
		this->bind();
		glUnmapBuffer(GL_ARRAY_BUFFER);
		this->release();
	}

	/**
	 * @brief copy data
	 * @param offset offset in bytes in the bufffer
	 * @param nb number of bytes to copy
	 * @param src source pointer
	 */
	inline void copy_data(uint32 offset, std::size_t nb, const void* src)
	{
		glBindBuffer(GL_ARRAY_BUFFER, id_);
		glBufferSubData(GL_ARRAY_BUFFER, offset, GLsizeiptr(nb), src);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
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

	inline void associate(GLuint attrib, int32 stride = 0, uint32 first = 0)
	{
		bind();
		glEnableVertexAttribArray(attrib);
		glVertexAttribPointer(attrib, vector_dimension(), GL_FLOAT, GL_FALSE, stride * vector_dimension() * 4, reinterpret_cast<GLvoid *>(first * uint64(vector_dimension()) * 4u));
		release();
	}
};

} // namespace rendering_pgl

} // namespace cgogn

#endif // CGOGN_RENDERING_VBO_H_
