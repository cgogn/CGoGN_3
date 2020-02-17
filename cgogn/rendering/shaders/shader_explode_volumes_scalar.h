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

#ifndef CGOGN_RENDERING_SHADERS_EXPLODE_VOLUMES_SCALAR_H_
#define CGOGN_RENDERING_SHADERS_EXPLODE_VOLUMES_SCALAR_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

enum ColorMap : int32
{
	BWR = 0,
	CWR,
	BCGYR,
	BGR
};


DECLARE_SHADER_CLASS(ExplodeVolumesScalar,CGOGN_STR(ExplodeVolumesScalar))

class CGOGN_RENDERING_EXPORT ShaderParamExplodeVolumesScalar : public ShaderParam
{
	void set_uniforms() override;

public:
	GLVec3 light_pos_;
	float32 explode_;
	VBO* vbo_pos_;
	VBO* vbo_center_;
	VBO* vbo_scalar_vol_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;
	int32 color_map_;
	int32 expansion_;
	float32 min_value_;
	float32 max_value_;

	using LocalShader = ShaderExplodeVolumesScalar;

	ShaderParamExplodeVolumesScalar(LocalShader* sh)
		: ShaderParam(sh), light_pos_(10, 100, 1000), explode_(0.8f),
		  vbo_pos_(nullptr),vbo_center_(nullptr),vbo_scalar_vol_(nullptr),
		  plane_clip_(0, 0, 0, 0),
		  plane_clip2_(0, 0, 0, 0),
		  color_map_(BWR), expansion_(0), min_value_(.0f), max_value_(1.0f)
	{
	}

	inline ~ShaderParamExplodeVolumesScalar() override
	{
	}

	inline void set_vbos(const std::vector<VBO*>& vbos) override
	{
		vbo_pos_ = vbos[0];
		vbo_center_ = vbos[1];
		vbo_scalar_vol_ = vbos[2];
	}
};

} // namespace rendering
} // namespace cgogn

#endif
