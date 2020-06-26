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

#ifndef CGOGN_UI_GL_VIEWER_H_
#define CGOGN_UI_GL_VIEWER_H_

#include <GL/gl3w.h>

#include <cgogn/ui/camera.h>
#include <cgogn/ui/cgogn_ui_export.h>
#include <cgogn/ui/inputs.h>

namespace cgogn
{

namespace ui
{

class CGOGN_UI_EXPORT GLViewer
{
public:
	GLViewer(Inputs* inputs);
	virtual ~GLViewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(GLViewer);

	inline bool need_redraw() const
	{
		return need_redraw_;
	}
	inline void request_update()
	{
		need_redraw_ = true;
	}

	inline const Camera& camera() const
	{
		return camera_;
	}
	inline void save_camera()
	{
		camera_saved_ = camera_;
	}
	inline void restore_camera()
	{
		camera_ = camera_saved_;
		need_redraw_ = true;
	}

	inline const rendering::GLMat4& projection_matrix() const
	{
		return camera_.projection_matrix();
	}
	inline const rendering::GLMat4d& projection_matrix_d() const
	{
		return camera_.projection_matrix_d();
	}
	inline const rendering::GLMat4& modelview_matrix() const
	{
		return camera_.modelview_matrix();
	}
	inline const rendering::GLMat4d& modelview_matrix_d() const
	{
		return camera_.modelview_matrix_d();
	}

	void set_manipulated_frame(MovingFrame* frame);

	inline void set_scene_radius(float64 radius)
	{
		camera_.set_scene_radius(radius);
	}
	inline void set_scene_center(const rendering::GLVec3d& center)
	{
		scene_center_ = center;
		camera_.set_pivot_point(scene_center_);
	}
	inline void set_scene_center(const rendering::GLVec3& center)
	{
		scene_center_ = center.cast<float64>();
		camera_.set_pivot_point(scene_center_);
	}
	inline void set_scene_pivot(const rendering::GLVec3d& pivot)
	{
		camera_.change_pivot_point(pivot);
	}
	inline void set_scene_pivot(const rendering::GLVec3& pivot)
	{
		camera_.change_pivot_point(pivot.cast<float64>());
	}
	inline void center_scene()
	{
		camera_.center_scene();
	}
	inline void show_entire_scene()
	{
		camera_.show_entire_scene();
	}

	inline int32 viewport_width() const
	{
		return viewport_width_;
	}
	inline int32 viewport_height() const
	{
		return viewport_height_;
	}

	inline bool shift_pressed() const
	{
		return inputs_->shift_pressed_;
	}
	inline bool control_pressed() const
	{
		return inputs_->control_pressed_;
	}
	inline bool alt_pressed() const
	{
		return inputs_->alt_pressed_;
	}
	inline bool meta_pressed() const
	{
		return inputs_->meta_pressed_;
	}

	inline bool mouse_button_pressed(int32 b) const
	{
		return (inputs_->mouse_buttons_ & (1 << b)) > 0;
	}

	inline int32 previous_mouse_x() const
	{
		return inputs_->previous_mouse_x_;
	}
	inline int32 previous_mouse_y() const
	{
		return inputs_->previous_mouse_y_;
	}

	virtual bool pixel_scene_position(int32 x, int32 y, rendering::GLVec3d& P) const = 0;
	virtual std::pair<rendering::GLVec3d, rendering::GLVec3d> pixel_ray(int32 x, int32 y) const = 0;

	inline void set_wheel_sensitivity(float64 s)
	{
		inputs_->wheel_sensitivity_ = s * 0.005;
	}
	inline void set_mouse_sensitivity(float64 s)
	{
		inputs_->mouse_sensitivity_ = s * 0.005;
	}
	inline void set_spin_sensitivity(float64 s)
	{
		inputs_->spin_sensitivity_ = s * 0.025;
	}

protected:
	virtual void resize_event(int32 viewport_width, int32 viewport_height);
	virtual void close_event();

	virtual void mouse_press_event(int32 button, int32 x, int32 y);
	virtual void mouse_release_event(int32 button, int32 x, int32 y);
	virtual void mouse_dbl_click_event(int32 button, int32 x, int32 y);
	virtual void mouse_move_event(int32 x, int32 y);
	virtual void mouse_wheel_event(float64 dx, float64 dy);
	virtual void key_press_event(int32 key_code);
	virtual void key_release_event(int32 key_code);

	inline bool obj_mode() const
	{
		return current_frame_ != &camera_;
	}
	void spin();

	Camera camera_;
	Camera camera_saved_;
	MovingFrame* current_frame_;
	rendering::GLVec3d scene_center_;
	int32 viewport_width_;
	int32 viewport_height_;
	float64 spinning_speed_;

	Inputs* inputs_;

	bool need_redraw_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_UI_GL_VIEWER_H_
