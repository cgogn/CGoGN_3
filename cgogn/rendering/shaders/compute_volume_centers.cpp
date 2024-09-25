/********************************************************************************
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

#include <cgogn/rendering/shaders/compute_volume_centers.h>

namespace cgogn
{

namespace rendering
{

namespace compute_center_shaders
{

/*****************************************************************************/
// ShaderComputeCenter1
/*****************************************************************************/

ShaderComputeCenter1* ShaderComputeCenter1::instance_ = nullptr;

ShaderComputeCenter1::ShaderComputeCenter1()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform usamplerBuffer vertex_ind;
		uniform samplerBuffer vertex_position;
		uniform float inv_h;

		out vec4 position;
		
		void main()
		{
			int ind_v = int(texelFetch(vertex_ind, 2 * gl_VertexID).r);
			uint ind_vol = texelFetch(vertex_ind, 2 * gl_VertexID + 1).r;

			vec2 d = vec2(1.0 / 1024.0, inv_h);
			vec2 coord_c = (-1.0 + d) + 2.0 * d * vec2(float(ind_vol % 1024u), float(ind_vol / 1024u));
			
			gl_Position = vec4(coord_c, 0, 1);
			position = vec4(texelFetch(vertex_position, ind_v).rgb, 1.0);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		in vec4 position;

		out vec4 frag_out;
		
		void main()
		{
			frag_out = position;
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("vertex_ind", "vertex_position", "inv_h");

	nb_attributes_ = 1;
}

void ShaderParamComputeCenter1::set_uniforms()
{
	shader_->set_uniforms_values(10, 11, 1.0f / float(tex_height_));
}

void ShaderParamComputeCenter1::bind_texture_buffers()
{
	vbo_position_->bind_texture_buffer(11);
}

void ShaderParamComputeCenter1::release_texture_buffers()
{
	vbo_position_->release_texture_buffer(11);
}

/*****************************************************************************/
// ShaderComputeCenter2
/*****************************************************************************/

ShaderComputeCenter2* ShaderComputeCenter2::instance_ = nullptr;

ShaderComputeCenter2::ShaderComputeCenter2()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform sampler2D tex_centers;

		out vec3 vbo_out;
		
		void main()
		{
			ivec2 icoord = ivec2(gl_VertexID % 1024, gl_VertexID / 1024);
			vec4 position4 = texelFetch(tex_centers, icoord, 0);
			vbo_out = position4.xyz / position4.w;
		}
	)";

	load_tfb1_bind(vertex_shader_source, {"vbo_out"});
	get_uniforms("tex_centers");
}

void ShaderParamComputeCenter2::set_uniforms()
{
	shader_->set_uniforms_values(tex_->bind(0));
}

} // namespace compute_center_shaders

ComputeVolumeCenterEngine::ComputeVolumeCenterEngine()
{
	param1_ = compute_center_shaders::ShaderComputeCenter1::generate_param();
	param2_ = compute_center_shaders::ShaderComputeCenter2::generate_param();
	param2_->tex_ = new Texture2D();
	param2_->tex_->allocate(0, 0, GL_RGBA32F, GL_RGBA, nullptr, GL_FLOAT);
	fbo_ = new FBO({param2_->tex_}, false, nullptr);
	tfb_ = new TFB_ComputeCenter(*(param2_.get()));
}

ComputeVolumeCenterEngine::~ComputeVolumeCenterEngine()
{
	delete fbo_;
	delete tfb_;
	delete param2_->tex_;
}

void ComputeVolumeCenterEngine::compute(VBO* vertex_position, MeshRender* renderer, VBO* volume_center)
{
	int32 h = (volume_center->size() + 1023) / 1024;
	fbo_->resize(1024, h);

	glDisable(GL_DEPTH_TEST);

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);

	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);

	glPointSize(1.001f);

	fbo_->bind();
	param1_->vbo_position_ = vertex_position;
	param1_->tex_height_ = h;
	param1_->bind();
	renderer->draw(VOLUMES_VERTICES);
	param1_->release();
	fbo_->release();

	glDisable(GL_BLEND);

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);

	tfb_->start(GL_POINTS, {volume_center});
	glDrawArrays(GL_POINTS, 0, volume_center->size());
	tfb_->stop();

	glEnable(GL_DEPTH_TEST);
}

} // namespace rendering

} // namespace cgogn
