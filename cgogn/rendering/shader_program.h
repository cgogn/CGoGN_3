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

#ifndef CGOGN_RENDERING_SHADERS_SHADERPROGRAM_H_
#define CGOGN_RENDERING_SHADERS_SHADERPROGRAM_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/core/utils/numerics.h>

#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLFunctions_3_3_Core>

#include <iostream>
#include <cassert>
#include <memory>

namespace cgogn
{

namespace rendering
{

// convenient conversion function
inline void* void_ptr(uint32 x)
{
	return reinterpret_cast<void*>(uint64_t(x));
}

// forward
class ShaderProgram;

class CGOGN_RENDERING_EXPORT ShaderParam
{
protected:

	ShaderProgram* shader_;
	std::unique_ptr<QOpenGLVertexArrayObject> vao_;
	QOpenGLFunctions_3_3_Core* ogl33_;

	virtual void set_uniforms() = 0;

public:

	ShaderParam(ShaderProgram* prg);

	virtual ~ShaderParam();

	inline ShaderProgram* get_shader()
	{
		return shader_;
	}

	/**
	 * @brief bind vao (and set uniform)
	 * @param with_uniforms ask to set uniforms
	 */
	void bind_vao_only(bool with_uniforms = true);

	/**
	 * @brief release vao
	 */
	void release_vao_only();

	/**
	 * @brief bind the shader set uniforms & matrices, bind vao
	 * @param proj projectiob matrix
	 * @param mv model-view matrix
	 */
	void bind(const QMatrix4x4& proj, const QMatrix4x4& mv);

	/**
	 * @brief release vao and shader
	 */
	void release();
};

class CGOGN_RENDERING_EXPORT ShaderProgram : protected QOpenGLFunctions_3_3_Core
{
protected:

	QOpenGLShaderProgram prg_;

	GLint unif_mv_matrix_;
	GLint unif_projection_matrix_;
	GLint unif_normal_matrix_;

	static std::vector<ShaderProgram*>* instances_;

public:

	static void register_instance(ShaderProgram* sh);

	static void clean_all();

	~ShaderProgram();

	/**
	 * @brief get the matrices uniforms (call after link)
	 */
	void get_matrices_uniforms();

	/**
	 * @brief set matrices (normal matrix computed if needed)
	 * @param proj projection matrix
	 * @param mv model view matrix
	 */
	void set_matrices(const QMatrix4x4& proj, const QMatrix4x4& mv);

	/**
	 * @brief set model-view matrice (normal matrix computed if needed)
	 * @param mv model view matrix
	 */
	void set_view_matrix(const QMatrix4x4& mv);

	/**
	 * @brief bind the shader
	 */
	inline void bind() { prg_.bind(); }

	/**
	 * @brief release the shader
	 */
	inline void release() { prg_.release(); }
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_SHADERPROGRAM_H_
