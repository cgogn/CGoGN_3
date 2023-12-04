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

#ifndef CGOGN_RENDERING_SHADERS_MESH_2D_EDGES_H_
#define CGOGN_RENDERING_SHADERS_MESH_2D_EDGES_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/texture.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(Mesh2DEdges, true, CGOGN_STR(Mesh2DEdges))

class CGOGN_RENDERING_EXPORT ShaderParamMesh2DEdges : public ShaderParam
{
	void set_uniforms() override;

	std::array<VBO*, 1> vbos_;

	inline void set_texture_buffer_vbo(uint32 i, VBO* vbo) override
	{
		vbos_[i] = vbo;
	}
	void bind_texture_buffers() override;
	void release_texture_buffers() override;

	enum VBOName : uint32
	{
		VERTEX_TC = 0,
	};

public:
	using ShaderType = ShaderMesh2DEdges;
	GLVec2 ratio_;
	GLColor color_;

	ShaderParamMesh2DEdges(ShaderType* sh) : ShaderParam(sh)
	{
		for (auto& v : vbos_)
			v = nullptr;
	}

	inline ~ShaderParamMesh2DEdges() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_NO_ILLUM_2D_H_
