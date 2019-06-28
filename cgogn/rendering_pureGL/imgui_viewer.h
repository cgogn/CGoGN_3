
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

#ifndef CGOGN_RENDERING_IMGUI_VIEWER_H_
#define CGOGN_RENDERING_IMGUI_VIEWER_H_

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <cgogn/rendering_pureGL/pure_gl_viewer.h>
#include <cgogn/rendering_pureGL/fbo.h>
#include <cgogn/rendering_pureGL/types.h>
#include <cgogn/rendering_pureGL/shaders/shader_frame2d.h>
#include <cgogn/rendering_pureGL/shaders/shader_fullscreen_texture.h>

namespace cgogn
{

namespace rendering_pgl
{

class ImGUIApp;

class CGOGN_RENDERING_PUREGL_EXPORT ImGUIViewer : public PureGLViewer
{
	friend class ImGUIApp;

public:

	ImGUIViewer(ImGUIViewer* share = nullptr);
	virtual ~ImGUIViewer() override;

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ImGUIViewer);

	virtual void resize_event(int32 w, int32 h);
	virtual void close_event();
	virtual void key_press_event(int32 key_code);
	virtual void key_release_event(int32 key_code);
	virtual void init() = 0;
	virtual void draw() = 0;

	virtual bool get_pixel_scene_position(int32 x, int32 y, GLVec3d& P) override;

	inline bool over_viewport(int32 x, int32 y)
	{
		y = fr_h_ - y;
		return (x >= vp_x_) && (x < vp_x_+vp_w_) && (y >= vp_y_) && (y < vp_y_+vp_h_);
	}

private:

	bool launch_init();
	void launch_loop();
	void set_viewer_geometry(float64 px, float64 py, float64 pw, float64 ph, float64 frame_width, float64 frame_height);
	void update_viewer_geometry(float64 frame_width, float64 frame_height);

protected:

	float64 fr_w_;
	float64 fr_h_;
	float64 vp_percent_x_;
	float64 vp_percent_y_;
	float64 vp_percent_width_;
	float64 vp_percent_height_;
	float64 last_click_time_;
	std::unique_ptr<ShaderFSTexture::Param> param_fst_;
	std::unique_ptr<FBO> fbo_;
	std::unique_ptr<Texture2D> tex_;

	void internal_init();
};

class CGOGN_RENDERING_PUREGL_EXPORT ImGUIApp
{
public:

	ImGUIApp();
	virtual ~ImGUIApp();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ImGUIApp);
	
//	void resize_event(int32 w, int32 h);
	void close_event();
	virtual bool interface();
	virtual void key_press_event(int32 k);
	virtual void key_release_event(int32 k);

	void set_window_size(int32 w, int32 h);
	void set_window_title(const std::string&  name);

	inline bool shift_pressed() const { return view_inputs_.shift_pressed_; }
	inline bool control_pressed() const { return view_inputs_.control_pressed_; }
	inline bool alt_pressed() const { return view_inputs_.alt_pressed_; }
	inline bool meta_pressed() const { return view_inputs_.meta_pressed_; }
	
	inline float32 device_pixel_ratio() const { return 1.0f; }
//	void set_double_click_timeout(float64 sec);

	void add_view(ImGUIViewer* view);
	void adapt_viewers_geometry();

	inline bool over_frame(int32 x, int32 y)
	{
		// assume positive params
		return (x < win_frame_width_) && (y < win_frame_width_);
	}

	inline void request_interface_update() { interface_need_redraw_ = true; }

	int launch();

protected:

	GLFWwindow* window_;
	ImGuiContext* context_;
	std::string win_name_;
	int32 win_frame_width_;
	int32 win_frame_height_;
	float64 interface_scaling_;
	bool show_imgui_;
	GLColor clear_color_;

	std::unique_ptr<ShaderFrame2d::Param> param_frame_;
	ViewerInputs view_inputs_;
	std::vector<ImGUIViewer*> viewers_;

	ImGUIViewer* focused_;
	bool interface_need_redraw_;

	bool IsItemActiveLastFrame();
};

} // namespace rendering_pgl

} // namespace cgogn

#endif // CGOGN_RENDERING_IMGUI_VIEWER_H_
