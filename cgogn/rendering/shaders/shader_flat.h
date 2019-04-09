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

#ifndef CGOGN_RENDERING_SHADERS_FLAT_H_
#define CGOGN_RENDERING_SHADERS_FLAT_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <cgogn/core/utils/unique_ptr.h>

#include <QOpenGLFunctions>
#include <QColor>

namespace cgogn
{

namespace rendering
{

// forward
template <bool CPV>
class ShaderParamFlat: public ShaderParam
{};

class CGOGN_RENDERING_EXPORT ShaderFlatGen : public ShaderProgram
{
	template <bool CPV> friend class ShaderParamFlat;

protected:

	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

	static const char* vertex_shader_source2_;
	static const char* fragment_shader_source2_;

	// uniform ids
	GLint unif_front_color_;
	GLint unif_back_color_;
	GLint unif_ambiant_color_;
	GLint unif_light_position_;
	GLint unif_bf_culling_;

public:

	using Self = ShaderFlatGen;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderFlatGen);

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_COLOR
	};

	/**
	 * @brief set current front color
	 * @param rgb
	 */
	void set_front_color(const QColor& rgb);

	/**
	 * @brief set current front color
	 * @param rgb
	 */
	void set_back_color(const QColor& rgb);

	/**
	 * @brief set current ambiant color
	 * @param rgb
	 */
	void set_ambiant_color(const QColor& rgb);

	/**
	 * @brief set light position relative to screen
	 * @param l light position
	 */
	void set_light_position(const QVector3D& l);

	/**
	 * @brief set light position relative to world
	 * @param l light position
	 * @param view_matrix
	 */
	void set_local_light_position(const QVector3D& l, const QMatrix4x4& view_matrix);

	void set_bf_culling(bool cull);

protected:

	ShaderFlatGen(bool color_per_vertex);
};

template <bool CPV>
class ShaderFlatTpl : public ShaderFlatGen
{
public:

	using Param = ShaderParamFlat<CPV>;
	static std::unique_ptr<Param> generate_param();

private:

	ShaderFlatTpl() : ShaderFlatGen(CPV) {}
	static ShaderFlatTpl* instance_;
};

template <bool CPV>
ShaderFlatTpl<CPV>* ShaderFlatTpl<CPV>::instance_ = nullptr;


// COLOR UNIFORM PARAM
template <>
class ShaderParamFlat<false> : public ShaderParam
{
protected:

	void set_uniforms() override
	{
		ShaderFlatGen* sh = static_cast<ShaderFlatGen*>(this->shader_);
		sh->set_front_color(front_color_);
		sh->set_back_color(back_color_);
		sh->set_ambiant_color(ambiant_color_);
		sh->set_light_position(light_pos_);
		sh->set_bf_culling(bf_culling_);
	}

public:

	using ShaderType = ShaderFlatTpl<false>;

	QColor front_color_;
	QColor back_color_;
	QColor ambiant_color_;
	QVector3D light_pos_;
	bool bf_culling_;

	ShaderParamFlat(ShaderFlatTpl<false>* sh) :
		ShaderParam(sh),
		front_color_(250, 0, 0),
		back_color_(0, 250, 0),
		ambiant_color_(5, 5, 5),
		light_pos_(10, 100, 1000),
		bf_culling_(false)
	{}

	void set_position_vbo(VBO* vbo_pos)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderFlatGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderFlatGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}
};


// COLOR PER VERTEX PARAM
template <>
class ShaderParamFlat<true> : public ShaderParam
{
protected:

	void set_uniforms() override
	{
		ShaderFlatGen* sh = static_cast<ShaderFlatGen*>(this->shader_);
		sh->set_ambiant_color(ambiant_color_);
		sh->set_light_position(light_pos_);
		sh->set_bf_culling(bf_culling_);
	}

public:

	using ShaderType = ShaderFlatTpl<true>;

	QColor ambiant_color_;
	QVector3D light_pos_;
	bool bf_culling_;

	ShaderParamFlat(ShaderFlatTpl<true>* sh) :
		ShaderParam(sh),
		ambiant_color_(5, 5, 5),
		light_pos_(10, 100, 1000),
		bf_culling_(false)
	{}

	void set_all_vbos(VBO* vbo_pos, VBO* vbo_color)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderFlatGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderFlatGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		// color
		vbo_color->bind();
		ogl->glEnableVertexAttribArray(ShaderFlatGen::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderFlatGen::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
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
		ogl->glEnableVertexAttribArray(ShaderFlatGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderFlatGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
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
		ogl->glEnableVertexAttribArray(ShaderFlatGen::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderFlatGen::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_color->release();
		vao_->release();
		shader_->release();
	}
};

template <bool CPV>
std::unique_ptr<typename ShaderFlatTpl<CPV>::Param> ShaderFlatTpl<CPV>::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderFlatTpl<CPV>();
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}

using ShaderFlat = ShaderFlatTpl<false>;
using ShaderFlatColor = ShaderFlatTpl<true>;

#if defined(CGOGN_USE_EXTERNAL_TEMPLATES) && !defined(CGOGN_RENDER_SHADERS_FLAT_CPP_)
extern template class CGOGN_RENDERING_EXPORT ShaderFlatTpl<false>;
extern template class CGOGN_RENDERING_EXPORT ShaderFlatTpl<true>;
extern template class CGOGN_RENDERING_EXPORT ShaderParamFlat<false>;
extern template class CGOGN_RENDERING_EXPORT ShaderParamFlat<true>;
#endif

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_FLAT_H_
