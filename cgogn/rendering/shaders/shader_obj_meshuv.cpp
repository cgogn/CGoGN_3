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

#define CGOGN_RENDER_SHADERS_OBJ_NORMAL_CPP_

#include <iostream>

#include <cgogn/rendering/shaders/shader_obj_meshuv.h>

namespace cgogn
{

namespace rendering
{

ShaderObjMeshUV* ShaderObjMeshUV::instance_ = nullptr;

ShaderObjMeshUV::ShaderObjMeshUV()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform vec2 ratio;
		uniform usamplerBuffer tex_coord_ind;
		uniform samplerBuffer vertex_tc;

		out vec2 tc;
		out vec2 normal;

		void main()
		{
			int ind = 3 * gl_InstanceID + gl_VertexID;
			int ind_tc = int(texelFetch(tex_coord_ind, ind).r);
			tc = texelFetch(vertex_tc, ind_tc).rg;
			vec2 P2 = vec2(2.0*tc-1.0)*ratio;
			gl_Position = vec4(P2.xy,0.5,1);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330

		out vec3 frag_out;

		void main()
		{
			frag_out = vec3(1,1,0);
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("tex_coord_ind","vertex_tc","ratio");
}

void ShaderParamObjMeshUV::set_uniforms()
{

	shader_->set_uniforms_values(10, 11,ratio_);
}

void ShaderParamObjMeshUV::bind_texture_buffers()
{
	vbos_[VERTEX_TC]->bind_texture_buffer(11);
}

void ShaderParamObjMeshUV::release_texture_buffers()
{
	vbos_[VERTEX_TC]->release_texture_buffer(11);
}

} // namespace rendering

} // namespace cgogn
