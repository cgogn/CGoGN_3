
/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/rendering/gl_viewer.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <iostream>

namespace cgogn
{

namespace rendering
{

ViewerInputs::ViewerInputs():
	wheel_sensitivity_(0.0025),
	mouse_sensitivity_(0.005),
	spin_sensitivity_(0.025),
	double_click_timeout_(0.3),
	need_redraw_(true),
	shift_pressed_(false),
	control_pressed_(false),
	alt_pressed_(false),
	meta_pressed_(false)
{}

PureGLViewer::PureGLViewer():
	vp_x_(0),
	vp_y_(0),
	vp_w_(0),
	vp_h_(0),
	inputs_(nullptr),
	need_redraw_(true)

{
	current_frame_ = &cam_;
}

PureGLViewer::~PureGLViewer()
{}

void PureGLViewer::manip(MovingFrame* fr)
{
	if (fr != nullptr)
	{
		current_frame_ = fr;
		inv_cam_ = cam_.frame_.inverse();
	}
	else
	{
		current_frame_ = &(cam_);
	}
}

void PureGLViewer::mouse_press_event(int32 button, float64 x, float64 y)
{
	if (button == 0)
	{
		current_frame_->is_moving_ = false;
		spinning_speed_ = 0;
		current_frame_->spin_ = Transfo3d::Identity();
		need_redraw_ = true;
	}
}

void PureGLViewer::mouse_release_event(int32 button, float64 x, float64 y)
{
	if (button == 0)
	{
		current_frame_->is_moving_ = (spinning_speed_ > inputs_->spin_sensitivity_);
		if (! current_frame_->is_moving_)
		{
			spinning_speed_ = 0;
			current_frame_->spin_ = Transfo3d::Identity();
		}
		need_redraw_ = true;
	}
}

void PureGLViewer::mouse_move_event(float64 x, float64 y)
{
	float64 dx = x - inputs_->last_mouse_x_;
	float64 dy = y - inputs_->last_mouse_y_;

	inputs_->last_mouse_x_ = x;
	inputs_->last_mouse_y_ = y;

	if ((inputs_->mouse_buttons_ & 1) && ((std::abs(dx) + std::abs(dy)) > 0.0))
	{
		GLVec3d axis(dy, dx, 0.0);
		spinning_speed_ = axis.norm();
		axis /= spinning_speed_;
		spinning_speed_ *= inputs_->mouse_sensitivity_;
		if (obj_mode())
		{
			Transfo3d sm(Eigen::AngleAxisd(2.0 * spinning_speed_, axis));
			current_frame_->spin_ = inv_cam_ *  sm * cam_.frame_;
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
	
	if (inputs_->mouse_buttons_ & 2)
	{
		float64 zcam = 1.0 / std::tan(cam_.field_of_view() / 2.0);
		float64 a = cam_.scene_radius() - cam_.frame_.translation().z() / zcam;
		if (obj_mode())
		{
			
			float64 tx = dx / vp_w_ * cam_.width() * a;
			float64 ty = - dy / vp_h_ * cam_.height() * a;
			Transfo3d ntr = inv_cam_ * Eigen::Translation3d(GLVec3d(tx, ty, 0.0)) * cam_.frame_;
			current_frame_->frame_ = ntr * current_frame_->frame_;
		}
		else
		{
			float64 nx = float64(dx) / vp_w_ * cam_.width() * a;
			float64 ny = - 1.0 * float64(dy) / vp_h_ * cam_.height() * a;
			cam_.frame_.translation().x() += 2 * nx;
			cam_.frame_.translation().y() += 2 * ny;
		}
		need_redraw_ = true;
	}
}

void PureGLViewer::spin()
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

void PureGLViewer::mouse_dbl_click_event(int32 buttons, float64 x, float64 y)
{
	if (inputs_->shift_pressed_)
	{
		GLVec3d P;
		if (get_pixel_scene_position(x, y, P))
		{
//			std::cout << P << std::endl;
			set_scene_pivot(P);
			need_redraw_ = true;
		}
	}
}

void PureGLViewer::mouse_wheel_event(float64 dx, float64 dy)
{
	if (dy != 0.0)
	{
		if (obj_mode())
		{
			auto ntr = inv_cam_ * Eigen::Translation3d(GLVec3d(0, 0, -0.0025 * dy)) * cam_.frame_;
			current_frame_->frame_ = ntr * current_frame_->frame_;
		}
		else
		{
			float64 zcam = 1.0 / std::tan(cam_.field_of_view() / 2.0);
			float64 a = cam_.scene_radius() - cam_.frame_.translation().z() / zcam / cam_.scene_radius();
			cam_.frame_.translation().z() -= inputs_->wheel_sensitivity_ * dy * std::max(0.1, a);
		}
		need_redraw_=true;
	}
}

} // namespace cgogn

} // namespace rendering
