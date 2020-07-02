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

#ifndef CGOGN_UI_CAMERA_H_
#define CGOGN_UI_CAMERA_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/ui/moving_frame.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace ui
{

class CGOGN_UI_EXPORT Camera : public MovingFrame
{
public:
	enum Type
	{
		PERSPECTIVE,
		ORTHOGRAPHIC
	};

private:
	Type type_;
	float64 field_of_view_;
	float64 aspect_ratio_; // width/height
	rendering::GLVec3d pivot_point_;
	float64 scene_radius_;
	float64 focal_distance_;
	mutable rendering::GLMat4d proj_d_;
	mutable rendering::GLMat4d mv_d_;
	mutable rendering::GLMat4 proj_;
	mutable rendering::GLMat4 mv_;

	rendering::GLMat4d perspective(float64 znear, float64 zfar) const;
	rendering::GLMat4d orthographic(float64 znear, float64 zfar) const;

public:
	inline Camera() : type_(PERSPECTIVE), field_of_view_(0.65), aspect_ratio_(1.0)
	{
	}

	inline float64 width() const
	{
		return (aspect_ratio_ > 1.0) ? aspect_ratio_ : 1.0;
	}
	inline float64 height() const
	{
		return (aspect_ratio_ > 1.0) ? 1.0 : 1.0 / aspect_ratio_;
	}

	inline void update_matrices()
	{
		float64 d = focal_distance_ - frame_.translation().z();
		float64 znear = std::max(0.1, d - 2.0 * scene_radius_);
		float64 zfar = d + 2.0 * scene_radius_;
		proj_d_ = ((type_ == PERSPECTIVE) ? perspective(znear, zfar) : orthographic(znear, zfar));
		proj_ = proj_d_.cast<float32>();

		rendering::Transfo3d m = Eigen::Translation3d(rendering::GLVec3d(0.0, 0.0, -focal_distance_)) * frame_ *
								 Eigen::Translation3d(-pivot_point_);
		mv_d_ = m.matrix();
		mv_ = mv_d_.cast<float32>();
	}

	inline void set_type(Type type)
	{
		type_ = type;
	}

	inline void set_field_of_view(float64 fov)
	{
		field_of_view_ = fov;
		focal_distance_ = scene_radius_ / std::tan(field_of_view_ / 2.0);
		update_matrices();
	}

	inline float64 field_of_view()
	{
		return field_of_view_;
	}

	inline void set_aspect_ratio(float64 aspect)
	{
		aspect_ratio_ = aspect;
		update_matrices();
	}

	inline void set_scene_radius(float64 radius)
	{
		scene_radius_ = radius;
		focal_distance_ = scene_radius_ / std::tan(field_of_view_ / 2.0);
		update_matrices();
	}

	inline void change_pivot_point(const rendering::GLVec3d& pivot)
	{
		rendering::GLVec3d pivot_shift = pivot - pivot_point_;
		frame_ *= Eigen::Translation3d(pivot_shift);
		pivot_point_ = pivot;
		update_matrices();
	}

	inline void set_pivot_point(const rendering::GLVec3d& pivot)
	{
		pivot_point_ = pivot;
		update_matrices();
	}

	inline void center_scene()
	{
		frame_.matrix().block<3, 1>(0, 3).setZero();
		update_matrices();
	}

	inline void show_entire_scene()
	{
		frame_.matrix().block<3, 1>(0, 3).setZero();
		update_matrices();
	}

	inline void reset()
	{
		frame_ = rendering::Transfo3d::Identity();
		spin_ = rendering::Transfo3d::Identity();
		update_matrices();
	}

	inline float64 scene_radius() const
	{
		return scene_radius_;
	}

	inline const rendering::GLVec3d& pivot_point() const
	{
		return pivot_point_;
	}

	inline const rendering::GLMat4d& projection_matrix_d() const
	{
		return proj_d_;
	}

	inline const rendering::GLMat4d& modelview_matrix_d() const
	{
		return mv_d_;
	}

	inline const rendering::GLMat4& projection_matrix() const
	{
		return proj_;
	}

	inline const rendering::GLMat4& modelview_matrix() const
	{
		return mv_;
	}
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_UI_CAMERA_H_
