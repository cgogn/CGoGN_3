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

#ifndef CGOGN_RENDERING_EBO_H_
#define CGOGN_RENDERING_EBO_H_

#include <GL/gl3w.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/types.h>
#include <iostream>
#include <string>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT EBO
{
protected:
	GLuint id_;
	GLuint id_tb_;
	std::size_t nb_;
	std::string name_;

public:
	inline EBO() : id_(0), id_tb_(0), nb_(0)
	{
	}

	inline void create()
	{
		glGenBuffers(1, &id_);
	}

	inline bool is_created()
	{
		return id_ != 0;
	}

	inline ~EBO()
	{
		glDeleteBuffers(1, &id_);
		id_ = 0;
	}

	inline void bind()
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id_);
	}

	inline void release()
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	inline GLint bind_tb(GLint unit)
	{
		if (id_tb_ == 0)
		{
			glGenTextures(1, &id_tb_);
			glBindTexture(GL_TEXTURE_BUFFER, id_tb_);
			glTexBuffer(GL_TEXTURE_BUFFER, GL_R32UI, id_);
			glBindTexture(GL_TEXTURE_BUFFER, 0);
		}

		glActiveTexture(GL_TEXTURE0 + unit);
		glBindTexture(GL_TEXTURE_BUFFER, id_tb_);
		return unit;
	}

	inline static void release_tb()
	{
		glBindTexture(GL_TEXTURE_BUFFER, 0);
	}

	inline void allocate(std::size_t nb_ind)
	{
		if (nb_ind != nb_) // only allocate when > ?
		{
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id_);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(nb_ind * sizeof(GLuint)), nullptr, GL_STATIC_DRAW);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
			nb_ = nb_ind;
		}
	}

	inline void allocate(const GLuint* indices, std::size_t nb_ind)
	{
		if (nb_ind != nb_) // only allocate when > ?
		{
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id_);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(nb_ind * sizeof(GLuint)), indices, GL_STATIC_DRAW);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
			nb_ = nb_ind;
		}
		else
		{
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id_);
			glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, GLsizeiptr(nb_ind * sizeof(GLuint)), indices);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	}

	/**
	 * @brief get and lock pointer on buffer memory, you must bind before
	 * @return  the pointer
	 */
	inline GLuint* lock_pointer()
	{
		return reinterpret_cast<GLuint*>(glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_READ_WRITE));
	}

	inline GLuint* lock_pointer_read()
	{
		return reinterpret_cast<GLuint*>(glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_READ_ONLY));
	}

	/**
	 * @brief release_pointer (must be binded)
	 */
	inline void release_pointer()
	{
		glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
	}

	/**
	 * @brief copy data
	 * @param offset offset in GLuint in the bufffer
	 * @param nb number of GLuint to copy
	 * @param src source pointer
	 */
	inline void copy_data(uint32 offset, std::size_t nb, const void* src)
	{
		this->bind();
		glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, offset * sizeof(GLuint), GLsizeiptr(nb) * sizeof(GLuint), src);
		this->release();
	}

	uint32 size() const
	{
		return uint32(nb_);
	}

	GLuint id() const
	{
		return id_;
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
};

inline std::ostream& operator<<(std::ostream& out, EBO& ebo)
{
	const int NB = 10;
	std::cout << "Debug EBO " << ebo.id() << " : " << ebo.name();
	std::cout << " ( " << ebo.size() << " indices) " << std::endl;
	ebo.bind();
	uint32* f = ebo.lock_pointer_read();
	for (int i = 0; i < NB; ++i)
	{
		std::cout << *f++ << " / ";
	}

	std::cout << std::endl << " . . . . " << std::endl;

	f += ebo.size() - 2 * NB;

	for (int i = 0; i < NB; ++i)
	{
		std::cout << *f++ << " / ";
	}

	std::cout << std::endl;
	ebo.release_pointer();
	ebo.release();
	std::cout << "-----end Debug EBO------" << std::endl;
	return out;
}

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_EBO_H_
