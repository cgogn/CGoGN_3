
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

#ifndef CGOGN_RENDERING_PURE_GL_VIEWER_H_
#define CGOGN_RENDERING_PURE_GL_VIEWER_H_

#include <GL/gl3w.h>

#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>
#include <cgogn/rendering_pureGL/camera.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <Eigen/SVD>

namespace cgogn
{

namespace rendering_pgl
{

struct CGOGN_RENDERING_PUREGL_EXPORT ViewerInputs
{
	ViewerInputs();

	float64 wheel_sensitivity_;
	float64 mouse_sensitivity_;
	float64 spin_sensitivity_;
	float64 double_click_timeout_;

	float64 last_mouse_x_;
	float64 last_mouse_y_;
	uint32 mouse_buttons_;

	bool need_redraw_;
	bool shift_pressed_;
	bool control_pressed_;
	bool alt_pressed_;
	bool meta_pressed_;
};

class CGOGN_RENDERING_PUREGL_EXPORT PureGLViewer
{
protected:

	Camera cam_;
	MovingFrame* current_frame_;
	Transfo3d inv_cam_;
	GLVec3d scene_center_;
	int32 vp_x_;
	int32 vp_y_;
	int32 vp_w_;
	int32 vp_h_;
	float64 spinning_speed_;

	ViewerInputs* inputs_;
	bool need_redraw_;

	void spin();

	inline void set_inputs(ViewerInputs *vin) { inputs_ = vin; }

public:

	PureGLViewer();
	virtual ~PureGLViewer();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(PureGLViewer);

	inline bool obj_mode() const { return  current_frame_ != &cam_; }
	inline void request_update() { need_redraw_ = true; }
	inline Camera camera() { return cam_; }
	inline GLMat4 get_projection_matrix() const { return cam_.get_projection_matrix(); }
	inline GLMat4 get_modelview_matrix() const { return cam_.get_modelview_matrix(); }

	inline void set_scene_radius(float64 radius) { cam_.set_scene_radius(radius); }
	inline void set_scene_center(const GLVec3d& center) { scene_center_ = center; cam_.set_pivot_point(scene_center_); }
	inline void set_scene_pivot(const GLVec3d& piv) { cam_.change_pivot_point(piv); }
	inline void set_scene_center(const GLVec3& center) { scene_center_ = center.cast<float64>(); cam_.set_pivot_point(scene_center_); }
	inline void set_scene_pivot(const GLVec3& piv) { cam_.change_pivot_point(piv.cast<float64>()); }
	inline void center_scene() { cam_.center_scene(); }
	inline void show_entire_scene() { cam_.show_entire_scene(); }

	inline int32 width() const { return vp_w_; }
	inline int32 height() const { return vp_h_; }
	inline void set_vp() { glViewport(vp_x_, vp_y_, vp_w_, vp_h_); }

	void manip(MovingFrame* fr);

	virtual void mouse_press_event(int32 button, float64 x, float64 y);
	virtual void mouse_release_event(int32 button, float64 x, float64 y);
	virtual void mouse_dbl_click_event(int32 button, float64 x, float64 y);
	virtual void mouse_move_event(float64 x, float64 y);
	virtual void mouse_wheel_event(float64 x, float64 y);

	virtual bool get_pixel_scene_position(int32 x, int32 y,GLVec3d& P) = 0;

	inline void set_wheel_sensitivity(float64 s) { inputs_->wheel_sensitivity_ = s * 0.005; }
	inline void set_mouse_sensitivity(float64 s) { inputs_->mouse_sensitivity_ = s * 0.005; }
	inline void set_spin_sensitivity(float64 s) { inputs_->spin_sensitivity_ = s * 0.025; }
};

} // namespace cgogn

} // namespace rendering_pgl

#endif // CGOGN_RENDERING_PURE_GL_VIEWER_H_
