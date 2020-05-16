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
#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(ExplodeVolumes, true, CGOGN_STR(ExplodeVolumes))

class CGOGN_RENDERING_EXPORT ShaderParamExplodeVolumes : public ShaderParam
{
	void set_uniforms() override;

	std::array<VBO*, 2> vbos_;
	inline void set_texture_buffer_vbo(uint32 i, VBO* vbo) override
	{
		vbos_[i] = vbo;
	}
	void bind_texture_buffers() override;
	void release_texture_buffers() override;

	enum VBOName : int32
	{
		VERTEX_POSITION = 0,
		VOLUME_CENTER
	};

public:
	GLColor color_;
	GLVec3 light_position_;
	float32 explode_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;

	using ShaderType = ShaderExplodeVolumes;

	inline ShaderParamExplodeVolumes(ShaderType* sh)
		: ShaderParam(sh), color_(0.9f, 0, 0, 1), light_position_(10, 100, 1000), explode_(0.9f),
		  plane_clip_(0, 0, 0, 0), plane_clip2_(0, 0, 0, 0)
	{
		for (auto& v : vbos_)
			v = nullptr;
	}

	inline ~ShaderParamExplodeVolumes() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif
