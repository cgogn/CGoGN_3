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

#ifndef CGOGN_RENDERING_SHADERS_CIRCLE_H_
#define CGOGN_RENDERING_SHADERS_CIRCLE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
//#include <cgogn/rendering/vbo.h>

namespace cgogn
{

namespace rendering
{

// forward
class ShaderParamCircle;

class CGOGN_RENDERING_EXPORT ShaderCircle : public ShaderProgram
{
	friend class ShaderParamCircle;

protected:
	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;
	GLint unif_color_;
	void set_locations();

public:
	using Self = ShaderCircle;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderCircle);

	using Param = ShaderParamCircle;
	static std::unique_ptr<Param> generate_param();

	void set_color(const GLColor& rgba);

private:
	ShaderCircle();
	static ShaderCircle* instance_;
};

// COLOR UNIFORM PARAM
class CGOGN_RENDERING_EXPORT ShaderParamCircle : public ShaderParam
{
protected:
	inline void set_uniforms() override
	{
		ShaderCircle* sh = static_cast<ShaderCircle*>(this->shader_);
		sh->set_color(color_);
	}

public:
	inline void set_color_vbo(VBO* vbo_color)
	{
		shader_->bind();
		vao_->bind();
		vbo_color->bind();
		glEnableVertexAttribArray(4);
		glVertexAttribPointer(4, vbo_color->vector_dimension(), GL_FLOAT, GL_FALSE, 0, nullptr);
		vbo_color->release();
		vao_->release();
		shader_->release();
	}

	using ShaderType = ShaderCircle;
	GLColor color_;

	inline ShaderParamCircle(ShaderCircle* sh) : ShaderParam(sh), color_(1.0, 1.0, 1.0, 1.0)
	{
	}

	inline virtual ~ShaderParamCircle() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_FLAT_H_
