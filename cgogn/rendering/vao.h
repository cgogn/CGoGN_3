
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
	}

	inline bool is_created()
	{
		return id_ != 0;
	}

	inline void create()
	{
		nb_ = 0;
		if (is_created())
			glDeleteVertexArrays(1, &id_);
		glGenVertexArrays(1, &id_);
	}

	inline VAO(const std::vector<std::tuple<GLint, VBO, GLint>>& params)
	{
		nb_ = std::get<1>(params[0]).size();
		glGenVertexArrays(1, &id_);
		glBindVertexArray(id_);
		for (const auto& p : params)
		{
			glBindBuffer(GL_ARRAY_BUFFER, std::get<1>(p).id());
			GLuint vid = GLuint(std::get<0>(p));
			glEnableVertexAttribArray(vid);
			glVertexAttribPointer(vid, GLint(std::get<1>(p).vector_dimension()), GL_FLOAT, GL_FALSE, 0, nullptr);
		}
		glBindVertexArray(0);
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
