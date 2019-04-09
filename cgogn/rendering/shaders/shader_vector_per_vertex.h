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

#ifndef CGOGN_RENDERING_SHADERS_VECTORPERVERTEX_H_
#define CGOGN_RENDERING_SHADERS_VECTORPERVERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <QOpenGLFunctions>
#include <QColor>

namespace cgogn
{

namespace rendering
{

class ShaderParamVectorPerVertex;

class CGOGN_RENDERING_EXPORT ShaderVectorPerVertex : public ShaderProgram
{
	friend class ShaderParamVectorPerVertex;

	static const char* vertex_shader_source_;
	static const char* geometry_shader_source_;
	static const char* fragment_shader_source_;

	// uniform ids
	GLint unif_color_;
	GLint unif_length_;

public:

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_NORMAL
	};

	using Param = ShaderParamVectorPerVertex;
	static std::unique_ptr<Param> generate_param();

	/**
	 * @brief set current color
	 * @param rgb
	 */
	void set_color(const QColor& rgb);

	/**
	 * @brief set length of normal
	 * @param l length
	 */
	void set_length(float32 l);

protected:

	ShaderVectorPerVertex();
	static ShaderVectorPerVertex* instance_;
};

class CGOGN_RENDERING_EXPORT ShaderParamVectorPerVertex : public ShaderParam
{
protected:

	inline void set_uniforms() override
	{
		ShaderVectorPerVertex* sh = static_cast<ShaderVectorPerVertex*>(this->shader_);
		sh->set_color(color_);
		sh->set_length(length_);
	}

public:

	using ShaderType = ShaderVectorPerVertex;

	QColor color_;
	float32 length_;

	ShaderParamVectorPerVertex(ShaderVectorPerVertex* sh) :
		ShaderParam(sh),
		color_(255, 255, 255),
		length_(1.0)
	{}

	void set_all_vbos(VBO* vbo_pos, VBO* vbo_vect)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderVectorPerVertex::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderVectorPerVertex::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		// vector vbo
		vbo_vect->bind();
		ogl->glEnableVertexAttribArray(ShaderVectorPerVertex::ATTRIB_NORMAL);
		ogl->glVertexAttribPointer(ShaderVectorPerVertex::ATTRIB_NORMAL, vbo_vect->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_vect->release();
		vao_->release();
		shader_->release();
	}

	void set_position_vbo(VBO* vbo_pos)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderVectorPerVertex::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderVectorPerVertex::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}

	void set_vector_vbo(VBO* vbo_vect)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_vect->bind();
		ogl->glEnableVertexAttribArray(ShaderVectorPerVertex::ATTRIB_NORMAL);
		ogl->glVertexAttribPointer(ShaderVectorPerVertex::ATTRIB_NORMAL, vbo_vect->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_vect->release();
		vao_->release();
		shader_->release();
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_VECTORPERVERTEX_H_
