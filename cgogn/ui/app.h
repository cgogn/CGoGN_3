/*******************************************************************************
* CGoGN                                                                        *
* Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_UI_WINDOW_H_
#define CGOGN_UI_WINDOW_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/ui/inputs.h>
#include <cgogn/rendering/shaders/shader_frame2d.h>

#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>

namespace cgogn
{

namespace ui
{

class View;
class Module;

class CGOGN_UI_EXPORT App
{
	friend class Module;

public:

	App();
	~App();

	void set_window_size(int32 w, int32 h);
	void set_window_title(const std::string& name);
	
	inline float32 device_pixel_ratio() const { return 1.0f; }

	View* add_view();
	inline View* current_view() const { return current_view_; }

	Module* module(const std::string& name) const;

	void init_modules();
	int launch();
	void stop();

private:

	void close_event();
    void adapt_views_geometry();
	inline bool over_frame(int32 x, int32 y) const
	{
		return (x < window_frame_width_) && (y < window_frame_width_);
	}

	GLFWwindow* window_;
	ImGuiContext* context_;

	std::string window_name_;

	int32 window_frame_width_;
	int32 window_frame_height_;

	float64 interface_scaling_;
	bool show_imgui_;
	
	std::unique_ptr<rendering::ShaderFrame2d::Param> param_frame_;

	Inputs inputs_;

	std::vector<std::unique_ptr<View>> views_;
	View* current_view_;

    mutable std::vector<Module*> modules_;
};

} // namespace cgogn

} // namespace ui

#endif // CGOGN_UI_WINDOW_H_
