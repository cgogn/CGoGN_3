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

#ifndef CGOGN_RENDERING_SHADERS_COLORPERVERTEX_H_
#define CGOGN_RENDERING_SHADERS_COLORPERVERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

// forward
class ShaderParamColorPerVertex;

class CGOGN_RENDERING_EXPORT ShaderColorPerVertex : public ShaderProgram
{
	friend class ShaderParamColorPerVertex;

protected:

	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

public:

	using Self = ShaderColorPerVertex;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderColorPerVertex);

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_COLOR
	};

	using Param = ShaderParamColorPerVertex;
	static std::unique_ptr<Param> generate_param();

protected:

	ShaderColorPerVertex();
	static ShaderColorPerVertex* instance_;
};

class CGOGN_RENDERING_EXPORT ShaderParamColorPerVertex : public ShaderParam
{
protected:

	inline void set_uniforms() override
	{}

public:
	using ShaderType = ShaderColorPerVertex;

	ShaderParamColorPerVertex(ShaderColorPerVertex* prg) : ShaderParam(prg)
	{}

	/**
	 * @brief set a vbo configuration
	 * @param vbo_pos pointer on position vbo (XYZ)
	 * @param vbo_col pointer on color vbo (RGB)
	 */
	void set_all_vbos(VBO* vbo_pos, VBO* vbo_color)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderColorPerVertex::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderColorPerVertex::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		// color vbo
		vbo_color->bind();
		ogl->glEnableVertexAttribArray(ShaderColorPerVertex::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderColorPerVertex::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_color->release();
		vao_->release();
		shader_->release();
	}

	void set_position_vbo(VBO* vbo_pos)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderColorPerVertex::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderColorPerVertex::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}

	void set_color_vbo(VBO* vbo_color)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_color->bind();
		ogl->glEnableVertexAttribArray(ShaderColorPerVertex::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderColorPerVertex::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_color->release();
		vao_->release();
		shader_->release();
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_COLORPERVERTEX_H_
