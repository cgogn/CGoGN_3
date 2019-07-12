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

#ifndef CGOGN_UI_CAMERA_H_
#define CGOGN_UI_CAMERA_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/ui/moving_frame.h>

#include <cgogn/rendering/types.h>
#include <cgogn/core/utils/numerics.h>

namespace cgogn
{

namespace ui
{

class CGOGN_UI_EXPORT Camera : public MovingFrame
{
public:

	enum Type { PERSPECTIVE, ORTHOGRAPHIC };

private:

	Type type_;
	float64 field_of_view_;
	float64 asp_ratio_; // width/height
	rendering::GLVec3d pivot_shift_;
	rendering::GLVec3d pivot_point_;
	float64 scene_radius_;
	float64 focal_dist_;
	mutable rendering::GLMat4d proj_;
	mutable rendering::GLMat4d mv_;
	mutable uint32 need_computing_;

	rendering::GLMat4d perspective(float64 znear, float64 zfar) const;
	rendering::GLMat4d ortho(float64 znear, float64 zfar) const;

public:

	inline Camera():
		type_(PERSPECTIVE),
		field_of_view_(0.78),
		asp_ratio_(1.0),
		pivot_shift_(0, 0, 0)
	{}

	inline float64 width() const { return (asp_ratio_ > 1.0) ? asp_ratio_ : 1.0; }
	inline float64 height() const { return (asp_ratio_ > 1.0) ? 1.0 : 1.0 / asp_ratio_; }
	
	inline void set_type(Type type) { type_ = type; need_computing_ = 3; }
	
	inline void set_field_of_view(float64 fov)
	{
		field_of_view_ = fov;
		focal_dist_ = scene_radius_/std::tan(field_of_view_ / 2.0);
		need_computing_ = 3;
	}

	inline float64 field_of_view() { return field_of_view_; }

	inline void set_aspect_ratio(float64 aspect)
	{
		asp_ratio_ = aspect;
		need_computing_ = 3;
	}
	
	inline void set_scene_radius(float64 radius)
	{
		scene_radius_ = radius;
		focal_dist_ = scene_radius_/std::tan(field_of_view_ / 2.0);
		need_computing_ = 3;
	}
	
	inline void change_pivot_point(const rendering::GLVec3d& piv)
	{
		pivot_shift_ = piv-pivot_point_;
		frame_ *= Eigen::Translation3d(pivot_shift_);
		pivot_point_ = piv;
		need_computing_ = 3;
	}

	inline void set_pivot_point(const rendering::GLVec3d& piv)
	{
		pivot_point_ = piv;
	}

	inline void center_scene()
	{
		this->frame_.matrix().block<3, 1>(0, 3).setZero();
		need_computing_ = 3;
	}
	
	inline void show_entire_scene() 
	{
		this->frame_.matrix().block<3, 1>(0, 3).setZero();
		need_computing_ = 3;
	}

	inline void reset()
	{
		this->frame_ = rendering::Transfo3d::Identity();
		this->spin_ = rendering::Transfo3d::Identity();
		need_computing_ = 3;
	}

	inline float64 scene_radius() const { return scene_radius_; }
	
	inline const rendering::GLVec3d& pivot_point() const { return pivot_point_; }

	inline rendering::GLMat4d projection_matrix_d() const
	{
//		Transfo3d tr = this->frame_ * Eigen::Translation3d(-pivot_shift_);
//		float64 d = focal_dist_ - (tr.translation()/*this->frame_.translation()-pivot_shift_*/).z();
		float64 d = focal_dist_ - this->frame_.translation().z();
		float64 znear = std::max(0.001, d - 2.0 * scene_radius_);
		float64 zfar = d + 2.0 * scene_radius_;
		proj_ = ((type_ == PERSPECTIVE) ? perspective(znear, zfar) : ortho(znear, zfar));
		return proj_;
	}

	inline rendering::GLMat4d modelview_matrix_d() const
	{
		rendering::Transfo3d m = Eigen::Translation3d(rendering::GLVec3d(0.0, 0.0, -focal_dist_)) * this->frame_ * Eigen::Translation3d(-pivot_point_);
		mv_ = m.matrix();
		return mv_;
	}

	inline rendering::GLMat4 projection_matrix() const
	{
		return projection_matrix_d().cast<float>();
	}

	inline rendering::GLMat4 modelview_matrix() const
	{
		return modelview_matrix_d().cast<float32>();
	}
};

} // namespace cgogn

} // namespace ui

#endif // CGOGN_UI_CAMERA_H_
