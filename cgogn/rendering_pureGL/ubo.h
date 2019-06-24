
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

#ifndef CGOGN_RENDERING_UBO_H_
#define CGOGN_RENDERING_UBO_H_

#include <GL/gl3w.h>
#include <string>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

class UBO
{
protected:
	GLuint id_;
	GLuint block_index_;
	GLuint binding_point_index_;

public:
	template<typename T>
	inline UBO(const T& data, const ShaderProgram& program, const std::string& uniform_block_name, GLuint binding_point)
	{
		glGenBuffers(1, &id_);
		glBindBuffer(GL_UNIFORM_BUFFER, id_);
		glBufferData(GL_UNIFORM_BUFFER, sizeof(T), &data, GL_DYNAMIC_COPY);
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
		block_index_ = glGetUniformBlockIndex(program.id(), uniform_block_name.c_str());
		binding_point_index_ = binding_point;
		glBindBufferBase(GL_UNIFORM_BUFFER, binding_point_index_, id_);
		glUniformBlockBinding(program.id(), block_index_, binding_point_index_);
	}

	inline ~UBO()
	{
		glDeleteBuffers(1,&id_);
	}

	template<typename T>
	inline void update_data(const T& data)
	{
		glBindBuffer(GL_UNIFORM_BUFFER, id_);
		glBufferSubData(GL_UNIFORM_BUFFER, 0, GLsizeiptr(sizeof(T)), &data);
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
	}

	template<typename T, typename S>
	inline void update_data(const T& data, const S& sub)
	{
		glBindBuffer(GL_UNIFORM_BUFFER, id_);
		glBufferSubData(GL_UNIFORM_BUFFER, &sub - &data, GLsizeiptr(sizeof(S)), &sub);
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
	}


	GLuint id() const
	{
		return id_;
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_VBO_H_
