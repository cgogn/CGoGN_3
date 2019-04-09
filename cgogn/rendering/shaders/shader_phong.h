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

#ifndef CGOGN_RENDERING_SHADERS_PHONG_H_
#define CGOGN_RENDERING_SHADERS_PHONG_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <cgogn/core/utils/unique_ptr.h>

#include <QColor>
#include <QVector3D>
#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

// forward
template <bool CPV>
class ShaderParamPhong: public ShaderParam
{};

class CGOGN_RENDERING_EXPORT ShaderPhongGen : public ShaderProgram
{
	template <bool CPV> friend class ShaderParamPhong;

protected:

	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

	static const char* vertex_shader_source_2_;
	static const char* fragment_shader_source_2_;

	// uniform ids
	GLint unif_front_color_;
	GLint unif_back_color_;
	GLint unif_ambiant_color_;
	GLint unif_spec_color_;
	GLint unif_spec_coef_;
	GLint unif_double_side_;
	GLint unif_light_position_;

public:

	using Self = ShaderPhongGen;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderPhongGen);

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_NORM,
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
	 * @brief set current specular color
	 * @param rgb
	 */
	void set_specular_color(const QColor& rgb);

	/**
	 * @brief set current specular coefficient
	 * @param rgb
	 */
	void set_specular_coef(float32 coef);

	/**
	 * @brief set double side option
	 * @param ts
	 */
	void set_double_side(bool ts);

	/**
	 * @brief set_light_position
	 * @param l
	 */
	void set_light_position(const QVector3D& l);

	/**
	 * @brief set light position relative to world
	 * @param l light position
	 * @param view_matrix
	 */
	void set_local_light_position(const QVector3D& l, const QMatrix4x4& view_matrix);

protected:

	ShaderPhongGen(bool color_per_vertex);
};


template <bool CPV>
class ShaderPhongTpl : public ShaderPhongGen
{
public:

	using Param = ShaderParamPhong<CPV>;
	static std::unique_ptr<Param> generate_param();

private:

	ShaderPhongTpl() : ShaderPhongGen(CPV) {}
	static ShaderPhongTpl* instance_;
};

template <bool CPV>
ShaderPhongTpl<CPV>* ShaderPhongTpl<CPV>::instance_ = nullptr;


// COLOR UNIFORM PARAM
template <>
class ShaderParamPhong<false> : public ShaderParam
{
protected:

	void set_uniforms() override
	{
		ShaderPhongGen* sh = static_cast<ShaderPhongGen*>(this->shader_);
		sh->set_front_color(front_color_);
		sh->set_back_color(back_color_);
		sh->set_ambiant_color(ambiant_color_);
		sh->set_specular_color(specular_color_);
		sh->set_specular_coef(specular_coef_);
		sh->set_double_side(double_side_);
		sh->set_light_position(light_position_);
	}

public:
	using ShaderType = ShaderPhongTpl<false>;

	QVector3D light_position_;
	QColor front_color_;
	QColor back_color_;
	QColor ambiant_color_;
	QColor specular_color_;
	float32 specular_coef_;
	bool double_side_;

	ShaderParamPhong(ShaderPhongTpl<false>* sh) :
		ShaderParam(sh),
		light_position_(10.0f, 100.0f, 1000.0f),
		front_color_(250, 0, 0),
		back_color_(0, 250, 5),
		ambiant_color_(5, 5, 5),
		specular_color_(100, 100, 100),
		specular_coef_(50.0f),
		double_side_(true)
	{}

	void set_all_vbos(VBO* vbo_pos, VBO* vbo_norm)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		// normal vbo
		vbo_norm->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_NORM);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_NORM, vbo_norm->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_norm->release();
		vao_->release();
		shader_->release();
	}

	void set_position_vbo(VBO* vbo_pos)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}

	void set_normal_vbo(VBO* vbo_norm)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_norm->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_NORM);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_NORM, vbo_norm->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_norm->release();
		vao_->release();
		shader_->release();
	}
};

// COLOR PER VERTEX PARAM
template <>
class ShaderParamPhong<true> : public ShaderParam
{
protected:

	void set_uniforms() override
	{
		ShaderPhongGen* sh = static_cast<ShaderPhongGen*>(this->shader_);
		sh->set_ambiant_color(ambiant_color_);
		sh->set_specular_color(specular_color_);
		sh->set_specular_coef(specular_coef_);
		sh->set_double_side(double_side_);
		sh->set_light_position(light_position_);
	}

public:
	using ShaderType = ShaderPhongTpl<true>;

	QVector3D light_position_;
	QColor ambiant_color_;
	QColor specular_color_;
	float32 specular_coef_;
	bool double_side_;

	ShaderParamPhong(ShaderPhongTpl<true>* sh) :
		ShaderParam(sh),
		light_position_(10.0f, 100.0f, 1000.0f),
		ambiant_color_(5, 5, 5),
		specular_color_(100, 100, 100),
		specular_coef_(50.0f),
		double_side_(true)
	{}

	void set_all_vbos(VBO* vbo_pos, VBO* vbo_norm, VBO* vbo_color)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		// position vbo
		vbo_pos->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		// normal vbo
		vbo_norm->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_NORM);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_NORM, vbo_norm->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_norm->release();
		// color  vbo
		vbo_color->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
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
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_POS);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_POS, vbo_pos->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_pos->release();
		vao_->release();
		shader_->release();
	}

	void set_normal_vbo(VBO* vbo_norm)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_norm->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_NORM);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_NORM, vbo_norm->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_norm->release();
		vao_->release();
		shader_->release();
	}

	void set_color_vbo(VBO* vbo_color)
	{
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();
		shader_->bind();
		vao_->bind();
		vbo_color->bind();
		ogl->glEnableVertexAttribArray(ShaderPhongGen::ATTRIB_COLOR);
		ogl->glVertexAttribPointer(ShaderPhongGen::ATTRIB_COLOR, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, 0);
		vbo_color->release();
		vao_->release();
		shader_->release();
	}
};


template <bool CPV>
std::unique_ptr<typename ShaderPhongTpl<CPV>::Param> ShaderPhongTpl<CPV>::generate_param()
{
	if (!instance_)
	{
		instance_ = new ShaderPhongTpl<CPV>;
		ShaderProgram::register_instance(instance_);
	}
	return cgogn::make_unique<Param>(instance_);
}


using ShaderPhong = ShaderPhongTpl<false>;
using ShaderPhongColor = ShaderPhongTpl<true>;


#if defined(CGOGN_USE_EXTERNAL_TEMPLATES) && !defined(CGOGN_RENDER_SHADERS_PHONG_CPP_)
extern template class CGOGN_RENDERING_EXPORT ShaderPhongTpl<false>;
extern template class CGOGN_RENDERING_EXPORT ShaderPhongTpl<true>;
extern template class CGOGN_RENDERING_EXPORT ShaderParamPhong<false>;
extern template class CGOGN_RENDERING_EXPORT ShaderParamPhong<true>;
#endif
} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_PHONG_H_
