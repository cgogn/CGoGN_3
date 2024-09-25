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

#ifndef CGOGN_RENDERING_SHADER_COPY_DEPTH_H_
#define CGOGN_RENDERING_SHADER_COPY_DEPTH_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/shaders/vbo.h>

#include <QOpenGLFunctions>
#include <QOpenGLTexture>

namespace cgogn
{

namespace rendering
{

// forward
class ShaderCopyDepth;

class CGOGN_RENDERING_EXPORT ShaderParamCopyDepth : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	using ShaderType = ShaderCopyDepth;
	QOpenGLTexture* texture_;
	GLuint depth_texture_sampler_;
	ShaderParamCopyDepth(ShaderCopyDepth* sh);
};

class CGOGN_RENDERING_EXPORT ShaderCopyDepth : public ShaderProgram
{
	friend class ShaderParamCopyDepth;

protected:
	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

	// uniform ids
	GLint unif_depth_texture_sampler_;

public:
	using Self = ShaderCopyDepth;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderCopyDepth);

	void set_depth_sampler(GLuint depth_samp);

	using Param = ShaderParamCopyDepth;

	static std::unique_ptr<Param> generate_param();

private:
	ShaderCopyDepth();
	static ShaderCopyDepth* instance_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_COPY_DEPTH_H_
