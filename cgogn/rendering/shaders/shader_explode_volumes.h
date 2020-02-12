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

#ifndef CGOGN_RENDERING_SHADERS_EXPLODE_VOLUMES_H_
#define CGOGN_RENDERING_SHADERS_EXPLODE_VOLUMES_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(ExplodeVolumes,CGOGN_STR(ExplodeVolumes))

class CGOGN_RENDERING_EXPORT ShaderParamExplodeVolumes : public ShaderParam
{
	void set_uniforms() override;

public:
	GLColor color_;
	GLVec3 light_pos_;
	float32 explode_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;
//	std::shared_ptr<VBO> vbo_pos_;
//	std::shared_ptr<VBO> vbo_center_;
	VBO* vbo_pos_;
	VBO* vbo_center_;


	using LocalShader = ShaderExplodeVolumes;

	inline ShaderParamExplodeVolumes(LocalShader* sh)
		: ShaderParam(sh), color_(color_front_default), light_pos_(10, 100, 1000), explode_(0.9f),
		  plane_clip_(0, 0, 0, 0), plane_clip2_(0, 0, 0, 0),
		  vbo_pos_(nullptr), vbo_center_(nullptr)
	{
	}

	inline ~ShaderParamExplodeVolumes() override
	{
	}

	inline void set_vbos(const std::vector<VBO*>& vbos) override
	{
		vbo_pos_ = vbos[0];
		vbo_center_ = vbos[1];
	}

};

} // namespace rendering
} // namespace cgogn

#endif
