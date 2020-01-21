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

#include <cgogn/rendering/shaders/compute_volume_centers.h>

namespace cgogn
{

namespace rendering
{
static const char* vertex_shader_source1 =
		R"(
		#version 330
		uniform usamplerBuffer vertex_volume;
		uniform samplerBuffer pos_vertex;
		uniform float inv_h;
		out vec4 P;
		void main()
		{
			vec2 d = vec2(1.0/1024,inv_h)
			int vid = gl_VertexID;
			int ind_v = int(texelFetch(vertex_volume, 2*gl_VertexID).r);
			uint ind_c = texelFetch(vertex_volume, 2*gl_VertexID+1).r;
			vec2 coord_c = (-1.0+d) + 2.0 * d * vec2(float(ind_c%1024),float(ind_c/1024),1);
			gl_Position = vec4(coord_c,0,1);
			P.xyz = texelFetch(pos_vertex, ind_v).rgb;
			P.w = 1.0;
		}
		)";

static const char* fragment_shader_source1 =
		R"(
		#version 330
		out vec4 fragOout;
		in vec4 P;
		void main()
		{
			fragOut = P;
		};
		)";

ShaderComputeCenter1* ShaderComputeCenter1::instance_ = nullptr;

ShaderComputeCenter1::ShaderComputeCenter1()
{
	load2_bind(vertex_shader_source1, fragment_shader_source1);
	add_uniforms("vertex_volume","pos_vertex","inv_h");
}

void ShaderParamComputeCenter1::set_uniforms()
{
	shader_->set_uniforms_values(1.0f/float(height_tex_), vbo_pos_->bind_tb(20));
}


static const char* vertex_shader_source2 =
		R"(
		#version 330
		uniform sampler2D TU;
		out vec3 vbo_out;

		const int w = 1024;

		void main()
		{
			ivec2 icoord = ivec2(gl_VertexID%w,gl_VertexID/w);
			vec4 P4 = texelFetch(TU,icoord,0);
			vbo_out = P4.xyz / P4.w;
		}
		)";


ShaderComputeCenter2* ShaderComputeCenter2::instance_ = nullptr;

ShaderComputeCenter2::ShaderComputeCenter2()
{
	load_tfb1_bind(vertex_shader_source2, {"vbo_out"});
	add_uniforms("TU");
}
void ShaderParamComputeCenter2::set_uniforms()
{
	shader_->set_uniforms_values(tex_->bind(0));
}



ComputeCenterEngine::ComputeCenterEngine()
{
	param1_ = ShaderComputeCenter1::generate_param();
	param2_ = ShaderComputeCenter2::generate_param();
	param2_->tex_ = new Texture2D();
	param2_->tex_->alloc(0,0,GL_RGB32F,GL_RGB,nullptr,GL_FLOAT);
	fbo_ = new FBO(std::vector<Texture2D*>{param2_->tex_},false, nullptr);
	tfb_ = new TFB_ComputeCenter(*param2_);
}

ComputeCenterEngine::~ComputeCenterEngine()
{
	delete tfb_;
	delete fbo_;
	delete param2_->tex_;
}

void ComputeCenterEngine::compute(VBO* pos, MeshRender* renderer, VBO* centers)
{
	int32 h = (centers->size()+1023)/1024;
	fbo_->resize(1024,h);
	fbo_->bind();
	glDisable(GL_DEPTH_TEST);
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE,GL_ONE);
	glPointSize(1.0f);
	param1_->vbo_pos_ = pos;
	param1_->height_tex_ = h;
	param1_->bind();
	renderer->draw(BUFFER_VOLUMES_VERTICES);
	param1_->release();
	fbo_->release();

	tfb_->start(GL_POINTS,{centers});
	glDrawArrays(GL_POINTS,0,centers->size());
	tfb_->stop();
}



} // namespace rendering

} // namespace cgogn
