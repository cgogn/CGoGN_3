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

#include <cgogn/ui/gl_viewer.h>

#include <GLFW/glfw3.h>

namespace cgogn
{

namespace ui
{

GLViewer::GLViewer(Inputs* inputs) : viewport_width_(0), viewport_height_(0), inputs_(inputs), need_redraw_(true)
{
	current_frame_ = &camera_;
}

GLViewer::~GLViewer()
{
}

void GLViewer::set_manipulated_frame(MovingFrame* frame)
{
	if (frame != nullptr)
	{
		current_frame_ = frame;
		inv_camera_ = camera_.frame_.inverse();
	}
	else
		current_frame_ = &camera_;
}

void GLViewer::resize_event(int32 viewport_width, int32 viewport_height)
{
	viewport_width_ = viewport_width;
	viewport_height_ = viewport_height;
	camera_.set_aspect_ratio(double(viewport_width_) / viewport_height_);
	need_redraw_ = true;
}

void GLViewer::close_event()
{
}

void GLViewer::mouse_press_event(int32 button, int32, int32)
{
	if (button == 0)
	{
		current_frame_->is_moving_ = false;
		spinning_speed_ = 0;
		current_frame_->spin_ = rendering::Transfo3d::Identity();
		need_redraw_ = true;
	}
}

void GLViewer::mouse_release_event(int32 button, int32, int32)
{
	if (button == 0)
	{
		current_frame_->is_moving_ = (spinning_speed_ > inputs_->spin_sensitivity_);
		if (!current_frame_->is_moving_)
		{
			spinning_speed_ = 0;
			current_frame_->spin_ = rendering::Transfo3d::Identity();
		}
		need_redraw_ = true;
	}
}

void GLViewer::mouse_dbl_click_event(int32 buttons, int32 x, int32 y)
{
	if (inputs_->shift_pressed_ && mouse_button_pressed(GLFW_MOUSE_BUTTON_LEFT))
	{
		rendering::GLVec3d P;
		if (pixel_scene_position(x, y, P))
		{
			set_scene_pivot(P);
			need_redraw_ = true;
		}
	}
}

void GLViewer::mouse_move_event(int32 x, int32 y)
{
	float64 dx = float64(x - inputs_->previous_mouse_x_);
	float64 dy = float64(y - inputs_->previous_mouse_y_);

	if (mouse_button_pressed(GLFW_MOUSE_BUTTON_LEFT) && ((std::abs(dx) + std::abs(dy)) > 0.0))
	{
		rendering::GLVec3d axis(dy, dx, 0.0);
		spinning_speed_ = axis.norm();
		axis /= spinning_speed_;
		spinning_speed_ *= inputs_->mouse_sensitivity_;
		if (obj_mode())
		{
			rendering::Transfo3d sm(Eigen::AngleAxisd(2.0 * spinning_speed_, axis));
			current_frame_->spin_ = inv_camera_ * sm * camera_.frame_;
			auto tr = current_frame_->frame_.translation().eval();
			current_frame_->frame_.translation().setZero();
			current_frame_->frame_ = current_frame_->spin_ * current_frame_->frame_;
			current_frame_->frame_.translation() = tr;
		}
		else
		{
			current_frame_->spin_ = Eigen::AngleAxisd(0.2 * spinning_speed_, axis);
			auto tr = current_frame_->frame_.translation().eval();
			current_frame_->frame_.translation().setZero();
			current_frame_->frame_ = Eigen::AngleAxisd(spinning_speed_, axis) * current_frame_->frame_;
			current_frame_->frame_.translation() = tr;
		}
		need_redraw_ = true;
	}

	if (mouse_button_pressed(GLFW_MOUSE_BUTTON_RIGHT))
	{
		float64 zcam = 1.0 / std::tan(camera_.field_of_view() / 2.0);
		float64 a = camera_.scene_radius() - camera_.frame_.translation().z() / zcam;
		if (obj_mode())
		{
			float64 tx = dx / viewport_width_ * camera_.width() * a;
			float64 ty = -dy / viewport_height_ * camera_.height() * a;
			rendering::Transfo3d ntr =
				inv_camera_ * Eigen::Translation3d(rendering::GLVec3d(tx, ty, 0.0)) * camera_.frame_;
			current_frame_->frame_ = ntr * current_frame_->frame_;
		}
		else
		{
			float64 nx = float64(dx) / viewport_width_ * camera_.width() * a;
			float64 ny = -1.0 * float64(dy) / viewport_height_ * camera_.height() * a;
			camera_.frame_.translation().x() += 2 * nx;
			camera_.frame_.translation().y() += 2 * ny;
		}
		need_redraw_ = true;
	}
}

void GLViewer::mouse_wheel_event(float64, float64 dy)
{
	if (dy != 0)
	{
		if (obj_mode())
		{
			auto ntr = inv_camera_ * Eigen::Translation3d(rendering::GLVec3d(0, 0, -inputs_->wheel_sensitivity_ * dy)) *
					   camera_.frame_;
			current_frame_->frame_ = ntr * current_frame_->frame_;
		}
		else
		{
			float64 zcam = 1.0 / std::tan(camera_.field_of_view() / 2.0);
			float64 a = camera_.scene_radius() - camera_.frame_.translation().z() / zcam / camera_.scene_radius();
			camera_.frame_.translation().z() -= inputs_->wheel_sensitivity_ * dy * std::max(0.1, a);
		}
		need_redraw_ = true;
	}
}

void GLViewer::key_press_event(int32 key_code)
{
	unused_parameters(key_code);
}

void GLViewer::key_release_event(int32 key_code)
{
	unused_parameters(key_code);
}

void GLViewer::spin()
{
	if (current_frame_->is_moving_)
	{
		auto tr = current_frame_->frame_.translation().eval();
		current_frame_->frame_.translation().setZero();
		current_frame_->frame_ = current_frame_->spin_ * current_frame_->frame_;
		current_frame_->frame_.translation() = tr;
		need_redraw_ = true;
	}
}

} // namespace ui

} // namespace cgogn
