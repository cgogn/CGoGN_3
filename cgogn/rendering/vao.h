
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

#ifndef CGOGN_RENDERING_VAO_H_
#define CGOGN_RENDERING_VAO_H_

#include <GL/gl3w.h>

#include <cgogn/rendering/vbo.h>

#include <string>
#include <tuple>
#include <vector>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT VAO
{
protected:
	GLuint id_;
	std::size_t nb_;
	std::string name_;

public:
	inline VAO() : id_(0), nb_(0)
	{
		nb_ = 0;
		glGenVertexArrays(1, &id_);
	}

	inline ~VAO()
	{
		if (id_ != 0)
			glDeleteVertexArrays(1, &id_);
	}

	inline void set_name(const std::string& name)
	{
		name_ = name;
		gl_debug_name(GL_VERTEX_ARRAY, id_, "VAO_" + name_);
	}

	inline GLuint id() const
	{
		return id_;
	}

	inline const std::string& name() const
	{
		return name_;
	}

	inline void bind()
	{
		glBindVertexArray(id_);
	}

	static inline void release()
	{
		glBindVertexArray(0);
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_VAO_H_
