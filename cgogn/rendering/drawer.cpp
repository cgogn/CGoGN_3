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

#include <cgogn/rendering/drawer.h>

namespace cgogn
{

namespace rendering
{

DisplayListDrawer::DisplayListDrawer() : current_size_(1.0f), current_aa_(true), current_ball_(true)
{
	vbo_pos_ = std::make_unique<VBO>(3);
	vbo_col_ = std::make_unique<VBO>(3);
}

DisplayListDrawer::~DisplayListDrawer()
{
}

void DisplayListDrawer::new_list()
{
	data_pos_.clear();
	data_col_.clear();
	begins_point_.clear();
	begins_round_point_.clear();
	begins_balls_.clear();
	begins_line_.clear();
	begins_bold_line_.clear();
	begins_face_.clear();
}

void DisplayListDrawer::begin(GLenum mode)
{
	switch (mode)
	{
	case GL_POINTS:
		if (current_ball_)
		{
			begins_balls_.push_back(PrimParam(data_pos_.size(), mode, current_size_, false));
			current_begin_ = &begins_balls_;
		}
		else if (current_size_ > 2.0f)
		{
			begins_round_point_.push_back(PrimParam(data_pos_.size(), mode, current_size_, current_aa_));
			current_begin_ = &begins_round_point_;
		}
		else
		{
			begins_point_.push_back(PrimParam(data_pos_.size(), mode, current_size_, false));
			current_begin_ = &begins_point_;
		}
		break;
	case GL_LINES:
	case GL_LINE_STRIP:
	case GL_LINE_LOOP:
		if (current_size_ > 1.0f)
		{
			begins_bold_line_.push_back(PrimParam(data_pos_.size(), mode, current_size_, current_aa_));
			current_begin_ = &begins_bold_line_;
		}
		else
		{
			begins_line_.push_back(PrimParam(data_pos_.size(), mode, 1.0f, current_aa_));
			current_begin_ = &begins_line_;
		}
		break;
	default:
		begins_face_.push_back(PrimParam(data_pos_.size(), mode, 1.0f, false));
		current_begin_ = &begins_face_;
		break;
	}
}

void DisplayListDrawer::end()
{
	current_begin_->back().nb = uint32(data_pos_.size() - current_begin_->back().begin);
}

void DisplayListDrawer::vertex3ff(float32 x, float32 y, float32 z)
{
	if (data_pos_.size() == data_col_.size())
	{
		if (data_col_.empty())
			data_col_.push_back(GLVec3{1.0f, 1.0f, 1.0f});
		else
			data_col_.push_back(data_col_.back());
	}
	data_pos_.push_back(GLVec3{x, y, z});
}

void DisplayListDrawer::color3ff(float32 r, float32 g, float32 b)
{
	if (data_pos_.size() == data_col_.size())
		data_col_.push_back(GLVec3{r, g, b});
	else
		data_col_.back() = GLVec3{r, g, b};
}

void DisplayListDrawer::end_list()
{
	uint32 nb_elts = uint32(data_pos_.size());

	if (nb_elts == 0)
		return;

	vbo_pos_->allocate(nb_elts, 3);
	float32* ptr = vbo_pos_->lock_pointer();
	std::memcpy(ptr, data_pos_[0].data(), nb_elts * 12);
	vbo_pos_->release_pointer();

	vbo_col_->allocate(nb_elts, 3);
	ptr = vbo_col_->lock_pointer();
	std::memcpy(ptr, data_col_[0].data(), nb_elts * 12);
	vbo_col_->release_pointer();

	// free memory
	data_pos_.clear();
	data_pos_.shrink_to_fit();
	data_col_.clear();
	data_col_.shrink_to_fit();
}

// void DisplayListDrawer::call_list(const QMatrix4x4& projection, const QMatrix4x4& modelview,
// QOpenGLFunctions_3_3_Core* ogl33)
//{

//	//classic rendering
//	if (!begins_point_.empty() || !begins_line_.empty() || !begins_face_.empty())
//	{
//		param_cpv_->bind(projection, modelview);

//		for (auto& pp : begins_point_)
//		{
//			glPointSize(pp.width);
//			glDrawArrays(pp.mode, pp.begin, pp.nb);
//		}

//		for (auto& pp : begins_line_)
//			glDrawArrays(pp.mode, pp.begin, pp.nb);

//		for (auto& pp : begins_face_)
//			glDrawArrays(pp.mode, pp.begin, pp.nb);

//		param_cpv_->release();
//	}

//	// balls
//	if (!begins_balls_.empty())
//	{
//		param_ps_->bind(projection,modelview);

//		for (auto& pp : begins_balls_)
//		{
//			ShaderPointSpriteColor* shader_ps_ = static_cast<ShaderPointSpriteColor*>(param_ps_->get_shader());
//			shader_ps_->set_size(pp.width);
//			glDrawArrays(pp.mode, pp.begin, pp.nb);
//		}
//		param_ps_->release();
//	}

//	// round points
//	if (!begins_round_point_.empty())
//	{
//		param_rp_->bind(projection, modelview);

//		for (auto& pp : begins_round_point_)
//		{
//			if (pp.aa)
//			{
//				glEnable(GL_BLEND);
//				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//			}
//			ShaderRoundPointColor* shader_rp_ = static_cast<ShaderRoundPointColor*>(param_rp_->get_shader());
//			shader_rp_->set_size(pp.width);
//			glDrawArrays(pp.mode, pp.begin, pp.nb);

//			if (pp.aa)
//				glDisable(GL_BLEND);
//		}
//		param_rp_->release();
//	}

//	// bold lines
//	if (!begins_bold_line_.empty())
//	{
//		param_bl_->bind(projection, modelview);

//		for (auto& pp : begins_bold_line_)
//		{
//			ShaderBoldLineColor* shader_bl_ = static_cast<ShaderBoldLineColor*>(param_bl_->get_shader());
//			shader_bl_->set_width(pp.width);
//			shader_bl_->set_color(QColor(255, 255, 0));

//			if (pp.aa)
//			{
//				glEnable(GL_BLEND);
//				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//			}

//			glDrawArrays(pp.mode, pp.begin, pp.nb);

//			if (pp.aa)
//				glDisable(GL_BLEND);
//		}

//		param_bl_->release();
//	}
//}

DisplayListDrawer::Renderer::Renderer(DisplayListDrawer* dr) : drawer_data_(dr)
{
	param_cpv_ = ShaderColorPerVertex::generate_param();
	param_bl_ = ShaderBoldLineColor::generate_param();
	param_rp_ = ShaderRoundPointColor::generate_param();
	param_ps_ = ShaderPointSpriteColor::generate_param();

	param_cpv_->set_vbos(dr->vbo_pos_.get(), dr->vbo_col_.get());
	param_bl_->set_vbos(dr->vbo_pos_.get(), dr->vbo_col_.get());
	param_rp_->set_vbos(dr->vbo_pos_.get(), dr->vbo_col_.get());
	param_ps_->set_vbos(dr->vbo_pos_.get(), dr->vbo_col_.get());
}

DisplayListDrawer::Renderer::~Renderer()
{
	param_cpv_.reset();
	param_bl_.reset();
	param_rp_.reset();
	param_ps_.reset();
}

void DisplayListDrawer::Renderer::draw(const GLMat4& projection, const GLMat4& modelview)
{

	// classic rendering
	if (!drawer_data_->begins_point_.empty() || !drawer_data_->begins_line_.empty() ||
		!drawer_data_->begins_face_.empty())
	{
		param_cpv_->bind(projection, modelview);

		for (const auto& pp : drawer_data_->begins_point_)
		{
			glPointSize(pp.width);
			glDrawArrays(pp.mode, pp.begin, pp.nb);
		}

		for (const auto& pp : drawer_data_->begins_line_)
			glDrawArrays(pp.mode, pp.begin, pp.nb);

		for (const auto& pp : drawer_data_->begins_face_)
			glDrawArrays(pp.mode, pp.begin, pp.nb);

		param_cpv_->release();
	}

	// balls
	if (!drawer_data_->begins_balls_.empty())
	{
		param_ps_->bind(projection, modelview);

		for (const auto& pp : drawer_data_->begins_balls_)
		{
			// get direct access to the shader to modify parameters while keeping the original param binded
			ShaderPointSpriteColor* shader_ps_ = static_cast<ShaderPointSpriteColor*>(param_ps_->get_shader());
			//			shader_ps_->set_size(pp.width);
			shader_ps_->set_uniform_value(2, pp.width);

			glDrawArrays(pp.mode, pp.begin, pp.nb);
		}

		param_ps_->release();
	}

	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// round points
	if (!drawer_data_->begins_round_point_.empty())
	{
		param_rp_->bind(projection, modelview);
		for (const auto& pp : drawer_data_->begins_round_point_)
		{
			// get direct access to the shader to modify parameters while keeping the original param binded
			ShaderRoundPointColor* shader_rp_ = static_cast<ShaderRoundPointColor*>(param_rp_->get_shader());
			// shader_rp_->set_size(pp.width);
			GLVec2 wd(pp.width / float32(viewport[2]), pp.width / float32(viewport[3]));
			shader_rp_->set_uniform_value(0, wd);

			if (pp.aa)
			{
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			}

			glDrawArrays(pp.mode, pp.begin, pp.nb);

			if (pp.aa)
				glDisable(GL_BLEND);
		}

		param_rp_->release();
	}

	// bold lines
	if (!drawer_data_->begins_bold_line_.empty())
	{
		param_bl_->bind(projection, modelview);

		for (const auto& pp : drawer_data_->begins_bold_line_)
		{
			// get direct access to the shader to modify parameters while keeping the original param binded
			ShaderBoldLineColor* shader_bl_ = static_cast<ShaderBoldLineColor*>(param_bl_->get_shader());
			GLVec2 wd(pp.width / float32(viewport[2]), pp.width / float32(viewport[3]));
			shader_bl_->set_uniform_value(0, wd);

			if (pp.aa)
			{
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			}

			glDrawArrays(pp.mode, pp.begin, pp.nb);

			if (pp.aa)
				glDisable(GL_BLEND);
		}

		param_bl_->release();
	}
}

} // namespace rendering

} // namespace cgogn
