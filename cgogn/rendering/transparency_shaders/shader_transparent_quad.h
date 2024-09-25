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

#ifndef CGOGN_RENDERING_SHADER_TRANSP_QUAD_H_
#define CGOGN_RENDERING_SHADER_TRANSP_QUAD_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/shaders/vbo.h>

#include <GLColor>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

// forward
class ShaderTranspQuad;

class CGOGN_RENDERING_EXPORT ShaderParamTranspQuad : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	using ShaderType = ShaderTranspQuad;
	GLuint rgba_texture_sampler_;
	GLuint depth_texture_sampler_;
	ShaderParamTranspQuad(ShaderTranspQuad* sh);
};

class CGOGN_RENDERING_EXPORT ShaderTranspQuad : public ShaderProgram
{
	friend class ShaderParamTranspQuad;

protected:
	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

	// uniform ids
	GLint unif_depth_texture_sampler_;
	GLint unif_rgba_texture_sampler_;

public:
	using Self = ShaderTranspQuad;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderTranspQuad);

	void set_rgba_sampler(GLuint rgba_samp);

	void set_depth_sampler(GLuint depth_samp);

	using Param = ShaderParamTranspQuad;

	static std::unique_ptr<Param> generate_param();

private:
	ShaderTranspQuad();
	static ShaderTranspQuad* instance_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_TRANSP_QUAD_H_
