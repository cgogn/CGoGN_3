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

#include <cgogn/rendering/shaders/outliner.h>

namespace cgogn
{

namespace rendering
{

namespace outline_shaders
{

/*****************************************************************************/
// ShaderMask
/*****************************************************************************/

ShaderMask* ShaderMask::instance_ = nullptr;

ShaderMask::ShaderMask()
{
	const char* vertex_shader_source = R"(
		#version 330
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		in vec3 vertex_position;

		void main()
		{
			gl_Position = projection_matrix * model_view_matrix * vec4(vertex_position, 1.0);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		out float frag_out;
		
		void main()
		{
			frag_out = 1.0;
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position");
}

void ShaderParamMask::set_uniforms()
{
}

/*****************************************************************************/
// ShaderSobel
/*****************************************************************************/

ShaderSobel* ShaderSobel::instance_ = nullptr;

ShaderSobel::ShaderSobel()
{
	const char* vertex_shader_source = R"(
		#version 330
		out vec2 tc;
		
		void main()
		{
			vec2 p = 2.0 * vec2(gl_VertexID % 2, gl_VertexID / 2);
			tc = p;
			p = 2.0 * p - 1.0;
			gl_Position = vec4(p, 0.0, 1.0);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform sampler2D TU;
		uniform vec2 texel_size;
		
		in vec2 tc;
		out float frag_out;
		
		void main()
		{
			vec2 N = tc;
			N.y -= texel_size.y;
			vec2 S = tc;
			S.y += texel_size.y;
			float v = abs(texture(TU, N).r - texture(TU, S).r);

			vec2 W = tc;
			W.x -= texel_size.x;
			vec2 E = tc;
			E.x += texel_size.x;
			float h = abs(texture(TU, E).r - texture(TU, W).r);
			frag_out = step(0.5, v + h);
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("TU", "texel_size");
}

void ShaderParamSobel::set_uniforms()
{
	shader_->set_uniforms_values(texture_->bind(0), GLVec2(1.0f / texture_->width(), 1.0f / texture_->height()));
}

/*****************************************************************************/
// ShaderBlur
/*****************************************************************************/

ShaderBlur* ShaderBlur::instance_ = nullptr;

ShaderBlur::ShaderBlur()
{
	const char* vertex_shader_source = R"(
		#version 330
		out vec2 tc;
		
		void main()
		{
			vec2 p = 2.0 * vec2(gl_VertexID % 2, gl_VertexID / 2);
			tc = p;
			p = 2.0 * p - 1.0;
			gl_Position = vec4(p, 0.0, 1.0);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform sampler2D TU;
		uniform vec2 texel_size_dir;
		
		in vec2 tc;
		out float frag_out;
		
		void main()
		{
			float v = 0.38774 * texture(TU, tc).r;
			vec2 p = tc - texel_size_dir;
			vec2 n = tc + texel_size_dir;
			v += 0.24477 * (texture(TU, p).r + texture(TU, n).r);
			p -= texel_size_dir;
			n+= texel_size_dir;
			frag_out = v + 0.06136 * (texture(TU, p).r + texture(TU, n).r);
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("TU", "texel_size_dir");
}

void ShaderParamBlur::set_uniforms()
{
	GLVec2 td = (pass_ % 2) ? GLVec2(1.0f / texture_->width(), 0.0) : GLVec2(0.0, 1.0f / texture_->height());
	shader_->set_uniforms_values(texture_->bind(0), td);
}

/*****************************************************************************/
// ShaderColorize
/*****************************************************************************/

ShaderColorize* ShaderColorize::instance_ = nullptr;

ShaderColorize::ShaderColorize()
{
	const char* vertex_shader_source = R"(
		#version 330
		out vec2 tc;
		
		void main()
		{
			vec2 p = 2.0 * vec2(gl_VertexID % 2, gl_VertexID / 2);
			tc = p;
			p = 2.0 * p - 1.0;
			gl_Position = vec4(p, 0.0, 1.0);
		}
	)";

	const char* fragment_shader_source = R"(
		#version 330
		uniform sampler2D TU_blur;
		uniform sampler2D TU_mask;
		uniform vec4 color;
		
		in vec2 tc;
		out vec3 frag_out;
		
		void main()
		{
			frag_out = (1.0 - 0.000001 * texture(TU_mask, tc).r) * texture(TU_blur, tc).r * color.rgb;
		}
	)";

	load(vertex_shader_source, fragment_shader_source);
	get_uniforms("TU_blur", "TU_mask", "color");
}

void ShaderParamColorize::set_uniforms()
{
	shader_->set_uniforms_values(texture_blur_->bind(0), texture_mask_->bind(1), color_);
}

} // namespace outline_shaders

Outliner* Outliner::instance_ = nullptr;

Outliner::Outliner()
{
	param_mask_ = outline_shaders::ShaderMask::generate_param();
	param_sobel_ = outline_shaders::ShaderSobel::generate_param();
	param_blur_ = outline_shaders::ShaderBlur::generate_param();
	param_colorize_ = outline_shaders::ShaderColorize::generate_param();

	Texture2D* t1 = new Texture2D();
	t1->allocate(0, 0, GL_R8, GL_RED, nullptr, GL_UNSIGNED_BYTE);
	fbo_mask_ = new FBO({t1}, false, nullptr);

	Texture2D* t2 = new Texture2D();
	t2->allocate(0, 0, GL_R8, GL_RED, nullptr, GL_UNSIGNED_BYTE);
	fbo_blur1_ = new FBO({t2}, false, nullptr);

	Texture2D* t3 = new Texture2D();
	t3->allocate(0, 0, GL_R8, GL_RED, nullptr, GL_UNSIGNED_BYTE);
	fbo_blur2_ = new FBO({t3}, false, nullptr);
}

Outliner::~Outliner()
{
	delete fbo_mask_->texture(0);
	delete fbo_mask_;
	delete fbo_blur1_->texture(0);
	delete fbo_blur1_;
	delete fbo_blur2_->texture(0);
	delete fbo_blur2_;
}

void Outliner::draw(VBO* position, MeshRender* renderer, const rendering::GLMat4& projection_matrix,
					const rendering::GLMat4& view_matrix, const GLColor& color)
{
	GLint prev_viewport[4];
	glGetIntegerv(GL_VIEWPORT, prev_viewport);
	prev_viewport[2] /= 4;
	prev_viewport[3] /= 4;

	if (fbo_mask_->width() != prev_viewport[2] || fbo_mask_->height() != prev_viewport[3])
	{
		fbo_mask_->resize(prev_viewport[2], prev_viewport[3]);
		fbo_blur1_->resize(prev_viewport[2], prev_viewport[3]);
		fbo_blur2_->resize(prev_viewport[2], prev_viewport[3]);
	}

	fbo_mask_->bind();
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	param_mask_->set_vbos({position});
	param_mask_->bind(projection_matrix, view_matrix);
	renderer->draw(TRIANGLES);
	param_mask_->release();

	fbo_blur1_->bind_no_release();
	glDisable(GL_DEPTH_TEST);
	param_sobel_->texture_ = fbo_mask_->texture(0);
	param_sobel_->bind();
	glDrawArrays(GL_TRIANGLES, 0, 3);
	param_sobel_->release();

	param_blur_->pass_ = 0;
	for (int i = 0; i < 8; ++i)
	{
		fbo_blur2_->bind_no_release();
		param_blur_->texture_ = fbo_blur1_->texture(0);
		param_blur_->bind();
		glDrawArrays(GL_TRIANGLES, 0, 3);
		param_blur_->pass_++;

		fbo_blur1_->bind_no_release();
		param_blur_->texture_ = fbo_blur2_->texture(0);
		param_blur_->bind();
		glDrawArrays(GL_TRIANGLES, 0, 3);
		param_blur_->pass_++;
	}

	fbo_mask_->release();

	param_colorize_->texture_mask_ = fbo_mask_->texture(0);
	param_colorize_->texture_blur_ = fbo_blur1_->texture(0);
	param_colorize_->color_ = color;
	param_colorize_->bind();
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
}

} // namespace rendering

} // namespace cgogn
