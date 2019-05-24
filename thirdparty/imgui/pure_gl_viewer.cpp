
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

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <iostream>
#include <cgogn/rendering_pureGL/pure_gl_viewer.h>

namespace cgogn
{
namespace rendering_pgl
{

PureGLViewer::PureGLViewer():
	vp_x_(0),
	vp_y_(0),
	vp_w_(0),
	vp_h_(0),
	wheel_sensitivity_(0.005),
	mouse_sensitivity_(0.005),
	spin_sensitivity_(0.025)
{
	current_frame_ = &cam_;
}

PureGLViewer::~PureGLViewer()
{}

void PureGLViewer::close_event()
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

void PureGLViewer::key_press_event(int32 key_code)
{
}

void PureGLViewer::key_release_event(int32 key_code)
{
}


void PureGLViewer::mouse_press_event(int32 button, float64 x, float64 y)
{
	std::cout <<" mouse_press_event: " << button <<" : "<< x <<" , "<< y << std::endl;
	if (button == 0)
	{
			current_frame_->is_moving_ = false;
			spinning_speed_ = 0;
			current_frame_->spin_ = Transfo3d::Identity();

	}
}

void PureGLViewer::mouse_release_event(int32 button, float64 x, float64 y)
{
	std::cout <<" mouse_release_event: " << button <<" : "<< x <<" , "<< y << std::endl;
	if (button == 0)
	{
		current_frame_->is_moving_ = (spinning_speed_ > spin_sensitivity_);
		if (! current_frame_->is_moving_)
		{
			spinning_speed_ = 0;
			current_frame_->spin_ = Transfo3d::Identity();
		}

	}
}


void PureGLViewer::mouse_move_event(float64 x, float64 y)
{
	float64 dx = x - last_mouse_x_;
	float64 dy = y - last_mouse_y_;

	last_mouse_x_ = x;
	last_mouse_y_ = y;

	if ((mouse_buttons_ & 1) && ((std::abs(dx)+ std::abs(dy))>0.0))
	{
		Vec3d axis(dy,dx,0.0);
		spinning_speed_ = axis.norm();
		axis /= spinning_speed_;
		spinning_speed_ *= mouse_sensitivity_;
		if (obj_mode())
		{
			Transfo3d sm(Eigen::AngleAxisd(2.0*spinning_speed_,axis));
			current_frame_->spin_ = inv_cam_ *  sm * cam_.frame_;
			auto tr = current_frame_->frame_.translation().eval();
			current_frame_->frame_.translation().setZero();
			current_frame_->frame_ = current_frame_->spin_ * current_frame_->frame_;
			current_frame_->frame_.translation() = tr;
		}
		else
		{
			current_frame_->spin_ = Eigen::AngleAxisd(0.2*spinning_speed_,axis);
			auto tr = current_frame_->frame_.translation().eval();
			current_frame_->frame_.translation().setZero();
			current_frame_->frame_ = Eigen::AngleAxisd(spinning_speed_,axis) * current_frame_->frame_;
			current_frame_->frame_.translation() = tr;
		}	
	}
	
	if (mouse_buttons_ & 2)
	{
		float64 zcam = 1.0/std::tan(cam_.field_of_view()/2.0);
		float64 a = cam_.scene_radius() - cam_.frame_.translation().z()/ zcam;
		if (obj_mode())
		{
			
			float64 tx = dx / vp_w_ * cam_.width() * a;
			float64 ty = - dy / vp_h_ * cam_.height() * a;
			Transfo3d ntr = inv_cam_ * Eigen::Translation3d(Vec3d(tx,ty,0.0)) * cam_.frame_;
			current_frame_->frame_ = ntr * current_frame_->frame_;
		}
		else
		{
			float64 nx = float64(dx) / vp_w_ * cam_.width() * a;
			float64 ny = - 1.0 * float64(dy) / vp_h_ * cam_.height() * a;
			cam_.frame_.translation().x() += 2*nx;
			cam_.frame_.translation().y() += 2*ny;
		}
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
	}
}


void PureGLViewer::mouse_dbl_click_event(int32 buttons, float64 x, float64 y)
{
}


void PureGLViewer::mouse_wheel_event(float64 dx, float64 dy)
{
	if (dy!=0.0)
	{
		if (obj_mode())
		{
			auto ntr = inv_cam_ * Eigen::Translation3d(Vec3d(0,0,0.0025*dy)) * cam_.frame_;
			current_frame_->frame_ = ntr * current_frame_->frame_;
		}
		else
		{
			float64 zcam = 1.0/std::tan(cam_.field_of_view()/2.0);
			float64 a = cam_.scene_radius() - cam_.frame_.translation().z()/zcam/cam_.scene_radius();
			std::cout << a << std::endl;
			cam_.frame_.translation().z() += wheel_sensitivity_*dy*std::max(0.1,a);
		}
	}
}

}
}
