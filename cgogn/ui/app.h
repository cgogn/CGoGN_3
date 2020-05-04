/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_UI_WINDOW_H_
#define CGOGN_UI_WINDOW_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/utils/type_traits.h>

#include <cgogn/ui/cgogn_ui_export.h>
#include <cgogn/ui/inputs.h>

#include <cgogn/rendering/shaders/shader_frame2d.h>

#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>
#include <imgui/imgui_internal.h>

#include <boost/synapse/emit.hpp>
#include <boost/synapse/thread_local_queue.hpp>
#include <thread>

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

	static inline float64 fps()
	{
		return fps_;
	}

	static float64 frame_time_;

	View* add_view();
	inline View* current_view() const
	{
		return current_view_;
	}

	template <typename FUNC>
	void foreach_view(const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, View*>::value,
					  "Wrong parameter type: given function should take a View*");
		for (const auto& v : views_)
			f(v.get());
	}

	Module* module(const std::string& name) const;

	void init_modules();
	int launch();
	void stop();

	using timer_tick = struct timer_tick_ (*)();

	template <typename FUNC>
	void start_timer(uint32 interval, FUNC stop_cond) const
	{
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		std::thread timer([this, interval, stop_cond]() {
			while (true)
			{
				if (stop_cond())
					return;
				std::this_thread::sleep_for(std::chrono::milliseconds(interval));
				if (stop_cond())
					return;
				boost::synapse::emit<App::timer_tick>(this);
			}
		});
		timer.detach();
	}

private:
	void close_event();
	void adapt_views_geometry();

	GLFWwindow* window_;
	ImGuiContext* context_;

	std::string window_name_;

	int32 window_width_;
	int32 window_height_;
	int32 framebuffer_width_;
	int32 framebuffer_height_;

	float32 interface_scaling_;

	float64 time_last_50_frames_;
	static float64 fps_;

	bool show_imgui_;
	bool show_demo_;

	std::unique_ptr<rendering::ShaderFrame2d::Param> param_frame_;

	Inputs inputs_;

	std::vector<std::unique_ptr<View>> views_;
	View* current_view_;

	mutable std::vector<Module*> modules_;

	std::shared_ptr<boost::synapse::thread_local_queue> tlq_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_UI_WINDOW_H_
