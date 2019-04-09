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

#ifndef CGOGN_RENDERING_SHADERS_SIMPLECOLOR_H_
#define CGOGN_RENDERING_SHADERS_SIMPLECOLOR_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <QOpenGLFunctions>
#include <QColor>

namespace cgogn
{

namespace rendering
{

class ShaderParamSimpleColor;

class CGOGN_RENDERING_EXPORT ShaderSimpleColor : public ShaderProgram
{
	friend class ShaderParamSimpleColor;

protected:

	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

	// uniform ids
	GLint unif_color_;

public:

	enum
	{
		ATTRIB_POS = 0
	};

	using Param = ShaderParamSimpleColor;
	static std::unique_ptr<Param> generate_param();

	/**
	 * @brief set current color
	 * @param rgb
	 */
	void set_color(const QColor& rgb);

protected:

	ShaderSimpleColor();
	static ShaderSimpleColor* instance_;
};

class CGOGN_RENDERING_EXPORT ShaderParamSimpleColor : public ShaderParam
{
protected:

	inline void set_uniforms() override
	{
		ShaderSimpleColor* sh = static_cast<ShaderSimpleColor*>(this->shader_);
		sh->set_color(color_);
	}

public:

	using ShaderType = ShaderSimpleColor;

	QColor color_;

	ShaderParamSimpleColor(ShaderSimpleColor* sh) :
		ShaderParam(sh),
		color_(255, 255, 255)
	{}

	inline void set_position_vbo(VBO* vbo_pos, uint32 stride = 0, uint32 first = 0)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderSimpleColor::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderSimpleColor::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, stride * vbo_pos->vector_dimension() * 4, void_ptr(first * vbo_pos->vector_dimension() * 4));
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_SIMPLECOLOR_H_
