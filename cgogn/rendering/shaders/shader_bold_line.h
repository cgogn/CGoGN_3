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

#ifndef CGOGN_RENDERING_SHADERS_BOLDLINE_H_
#define CGOGN_RENDERING_SHADERS_BOLDLINE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <cgogn/core/utils/unique_ptr.h>

#include <QColor>
#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

// forward
template <bool CPV>
class ShaderParamBoldLine : public ShaderParam
{};

class CGOGN_RENDERING_EXPORT ShaderBoldLineGen : public ShaderProgram
{
	template <bool CPV> friend class ShaderParamBoldLine;

protected:

	static const char* vertex_shader_source_;
	static const char* geometry_shader_source_;
	static const char* fragment_shader_source_;

	static const char* vertex_shader_source2_;
	static const char* geometry_shader_source2_;
	static const char* fragment_shader_source2_;

	// uniform ids
	GLint unif_color_;
	GLint unif_width_;
	GLint unif_plane_clip_;
	GLint unif_plane_clip2_;

public:

	using Self = ShaderBoldLineGen;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderBoldLineGen);

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_COLOR
	};

	/**
	 * @brief set current color
	 * @param rgb
	 */
	void set_color(const QColor& rgb);

	/**
	 * @brief set the width of lines (call before each draw)
	 * @param w width in pixel
	 */
	void set_width(float32 w);

	/**
	 * @brief set_plane_clip
	 * @param plane
	 */
	void set_plane_clip(const QVector4D& plane);

	/**
	 * @brief set_plane_clip2
	 * @param plane
	 */
	void set_plane_clip2(const QVector4D& plane);


protected:

	ShaderBoldLineGen(bool color_per_vertex);
};

template <bool CPV>
class ShaderBoldLineTpl : public ShaderBoldLineGen
{
public:

	using Param = ShaderParamBoldLine<CPV>;
	static std::unique_ptr<Param> generate_param();

private:

	ShaderBoldLineTpl() : ShaderBoldLineGen(CPV) {}
	static ShaderBoldLineTpl* instance_;
};

template <bool CPV>
ShaderBoldLineTpl<CPV>* ShaderBoldLineTpl<CPV>::instance_ = nullptr;


// COLOR UNIFORM VERSION
template <>
class ShaderParamBoldLine<false> : public ShaderParam
{
protected:

	void set_uniforms() override
	{
		ShaderBoldLineGen* sh = static_cast<ShaderBoldLineGen*>(this->shader_);
		sh->set_color(color_);
		sh->set_width(width_);
		sh->set_plane_clip(plane_clip_);
		sh->set_plane_clip2(plane_clip2_);
	}

public:
	using ShaderType = ShaderBoldLineTpl<false>;

	QColor color_;
	float32 width_;
	QVector4D plane_clip_;
	QVector4D plane_clip2_;

	ShaderParamBoldLine(ShaderBoldLineTpl<false>* sh) :
		ShaderParam(sh),
		color_(255, 255, 255),
		width_(2.0f),
		plane_clip_(0,0,0,0),
		plane_clip2_(0,0,0,0)
	{}

	void set_position_vbo(VBO* vbo_pos)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderBoldLineGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderBoldLineGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}
};

// COLOR PER VERTEX VERSION
template <>
class ShaderParamBoldLine<true> : public ShaderParam
{
protected:

	void set_uniforms() override
	{
		ShaderBoldLineGen* sh = static_cast<ShaderBoldLineGen*>(this->shader_);
		sh->set_width(width_);
		sh->set_plane_clip(plane_clip_);
		sh->set_plane_clip2(plane_clip2_);
	}

public:

	using ShaderType = ShaderBoldLineTpl<true>;
	using Self = ShaderParamBoldLine<true>;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderParamBoldLine);

	float32 width_;
	QVector4D plane_clip_;
	QVector4D plane_clip2_;

	ShaderParamBoldLine(ShaderBoldLineTpl<true>* sh) :
		ShaderParam(sh),
		width_(2.0f),
		plane_clip_(0,0,0,0),
		plane_clip2_(0,0,0,0)
	{}

	void set_all_vbos(VBO* vbo_pos, VBO* vbo_color)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderBoldLineGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderBoldLineGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		// color vbo
		vbo_color->bind();
		ogl->glEnableVertexAttribArray(ShaderBoldLineGen::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderBoldLineGen::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
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
		ogl->glEnableVertexAttribArray(ShaderBoldLineGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderBoldLineGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
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
		ogl->glEnableVertexAttribArray(ShaderBoldLineGen::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderBoldLineGen::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_color->release();
		vao_->release();
		shader_->release();
	}
};


template <bool CPV>
std::unique_ptr<typename ShaderBoldLineTpl<CPV>::Param> ShaderBoldLineTpl<CPV>::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderBoldLineTpl<CPV>;
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

using ShaderBoldLine = ShaderBoldLineTpl<false>;
using ShaderBoldLineColor = ShaderBoldLineTpl<true>;

#if defined(CGOGN_USE_EXTERNAL_TEMPLATES) && !defined(CGOGN_RENDER_SHADERS_BOLD_LINE_CPP_)
extern template class CGOGN_RENDERING_EXPORT ShaderBoldLineTpl<false>;
extern template class CGOGN_RENDERING_EXPORT ShaderBoldLineTpl<true>;
extern template class CGOGN_RENDERING_EXPORT ShaderParamBoldLine<false>;
extern template class CGOGN_RENDERING_EXPORT ShaderParamBoldLine<true>;
#endif

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_BOLDLINE_H_
