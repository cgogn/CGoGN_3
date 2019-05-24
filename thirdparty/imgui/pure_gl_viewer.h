
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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <cgogn/rendering_pureGL/camera.h>

#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>

namespace cgogn
{
namespace rendering_pgl
{

class CGOGN_RENDERING_PUREGL_EXPORT PureGLViewer
{
//	enum KeyCode : uint32
//	{
//		SHIFT= 0x1000,
//		CONTROL,
//		ALT,
//		META,
//		ARROW_LEFT,
//		ARROW_RIGHT,
//		ARROW_UP,
//		ARROW_DOWN,
//		PAGE_UP,
//		PAGE_DOWN,
//		HOME,
//		END,
//		BACKSPACE,
//		DEL,
//		ENTER,
//		TAB,
//		ESC,
//		F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12
//	};
protected:
	Camera cam_;
	MovingFrame* current_frame_;
	Transfo3d inv_cam_;
	int32 vp_x_;
	int32 vp_y_;
	int32 vp_w_;
	int32 vp_h_;
	float64 spinning_speed_;

	float64 wheel_sensitivity_;
	float64 mouse_sensitivity_;
	float64 spin_sensitivity_;

	float64 last_mouse_x_;
	float64 last_mouse_y_;
	uint32 mouse_buttons_;

protected:
	bool shift_pressed_;
	bool control_pressed_;
	bool alt_pressed_;
	bool meta_pressed_;
	void spin();

public:
	PureGLViewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(PureGLViewer);
	virtual ~PureGLViewer();

	inline bool obj_mode() const { return  current_frame_ != &cam_;}

	inline Camera camera() {return cam_;}

	inline GLMat4 get_projection_matrix() const
	{
		return cam_.get_projection_matrix();
	}

	inline GLMat4 get_modelview_matrix() const
	{
		return cam_.get_modelview_matrix();
	}

	inline void set_scene_radius(float64 radius) { cam_.set_scene_radius(radius); }
	inline void set_scene_center(const GLVec3d& center) {cam_.set_scene_center(center); }

	inline void center_scene() { cam_.center_scene(); }
	inline void show_entire_scene() { cam_.show_entire_scene(); }

	void manip(MovingFrame* fr);

	virtual void mouse_press_event(int32 button, float64 x, float64 y);
	virtual void mouse_release_event(int32 button, float64 x, float64 y);
	virtual void mouse_dbl_click_event(int32 button, float64 x, float64 y);
	virtual void mouse_move_event(float64 x, float64 y);
	virtual void mouse_wheel_event(float64 x, float64 y);
	virtual void key_press_event(int32 key_code);
	virtual void key_release_event(int32 key_code);
	virtual void close_event();

	inline void set_wheel_sensitivity(float64 s) { wheel_sensitivity_ = s*0.005; }
	inline void set_mouse_sensitivity(float64 s) { mouse_sensitivity_ = s*0.005; }
	inline void set_spin_sensitivity(float64 s) { spin_sensitivity_ = s*0.025; }

	virtual bool init() = 0;
	virtual void draw() = 0;
};	

}
}
#endif
