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

#ifndef CGOGN_RENDERING_SHADERS_XXXX_H_
#define CGOGN_RENDERING_SHADERS_XXXX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

// forward
class ShaderParamXXXX;

class CGOGN_RENDERING_EXPORT ShaderXXXX : public ShaderProgram
{
public:
	using Self = ShaderXXXX;
	using Param = ShaderParamXXXX;
	friend Param;

protected:
	ShaderXXXX();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderXXXX);

	void set_locations() override;
	static Self* instance_;

public:
	inline static std::unique_ptr<Param> generate_param()
	{
		if (!instance_)
		{
			instance_ = new Self();
			ShaderProgram::register_instance(instance_);
		}
		return cgogn::make_unique<Param>(instance_);
	}
};

class CGOGN_RENDERING_EXPORT ShaderParamXXXX : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(front_color_, back_color_, ambiant_color_, light_pos_, bf_culling_);
	}

public:
	GLColor front_color_;
	GLColor back_color_;
	GLColor ambiant_color_;
	GLVec3 light_pos_;
	bool bf_culling_;

	using LocalShader = ShaderXXXX;

	ShaderParamXXXX(LocalShader* sh)
		: ShaderParam(sh),

		  light_pos_(10, 100, 1000), bf_culling_(false)
	{
	}

	inline ~ShaderParamXXXX() override
	{
	}

};

} // namespace rendering

} // namespace cgogn

#endif
