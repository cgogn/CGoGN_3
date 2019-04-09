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

#ifndef CGOGN_RENDERING_SHADERS_TEXTURE_H_
#define CGOGN_RENDERING_SHADERS_TEXTURE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo.h>

#include <cgogn/core/utils/unique_ptr.h>

#include <QOpenGLTexture>

class QOpenGLTexture;

namespace cgogn
{

namespace rendering
{

class ShaderTexture;

class CGOGN_RENDERING_EXPORT ShaderParamTexture : public ShaderParam
{
protected:

	void set_uniforms();

public:

	using ShaderType = ShaderTexture;

	QOpenGLTexture* texture_;

	ShaderParamTexture(ShaderTexture* sh);

	void set_vbo(VBO* vbo_pos, VBO* vbo_tc);
};

class CGOGN_RENDERING_EXPORT ShaderTexture : public ShaderProgram
{
	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

public:

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_TC
	};

	using Param = ShaderParamTexture;

	/**
	 * @brief generate shader parameter object
	 * @return pointer
	 */
	inline static std::unique_ptr<Param> generate_param()
	{
		if (!instance_)
		{
			instance_ = new ShaderTexture();
			ShaderProgram::register_instance(instance_);
		}
		return cgogn::make_unique<Param>(instance_);
	}

protected:

	ShaderTexture();
	static ShaderTexture* instance_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_TEXTURE_H_
