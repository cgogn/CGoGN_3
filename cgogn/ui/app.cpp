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

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

#include <imgui/imgui_internal.h>

namespace cgogn
{

namespace ui
{

static void glfw_error_callback(int error, const char* description)
{
	std::cerr << "Glfw Error " << error << ": " << description << std::endl;
}

App::App():
	window_(nullptr),
    context_(nullptr),
	window_name_("CGoGN"),
	window_frame_width_(512),
	window_frame_height_(512),
	interface_scaling_(1.0),
	show_imgui_(true),
	focused_(nullptr)
{
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
		std::cerr << "Failed to initialize GFLW!" << std::endl;

	// GL 3.3 + GLSL 150 + Core Profile
	const char* glsl_version = "#version 150";
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	window_ = glfwCreateWindow(window_frame_width_, window_frame_height_, window_name_.c_str(), nullptr, nullptr);
	if (window_ == nullptr)
		std::cerr << "Failed to create Window!" << std::endl;

	glfwMakeContextCurrent(window_);

	bool err = gl3wInit() != 0;
	if (err)
		std::cerr << "Failed to initialize OpenGL loader!" << std::endl;
	
	glfwSwapInterval(1); // Enable vsync

	IMGUI_CHECKVERSION();
	context_ = ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui::StyleColorsDark(); //ImGui::StyleColorsClassic();
	ImGui_ImplGlfw_InitForOpenGL(window_, true);
	ImGui_ImplOpenGL3_Init(glsl_version);
	glfwSetWindowUserPointer(window_, this);

	std::cout << glGetString(GL_VENDOR) << std::endl;
	std::cout << glGetString(GL_RENDERER) << std::endl;
	std::cout << glGetString(GL_VERSION) << std::endl;

	int x, y;
	glfwGetWindowSize(window_, &x, &y);
	std::cout << "Window size: " << x << ", " << y << std::endl;
	glfwGetFramebufferSize(window_, &x, &y);
	std::cout << "Framebuffer size: " << x << ", " << y << std::endl;

	glfwSetWindowSizeCallback(window_, [] (GLFWwindow* wi, int, int)
	{
		App* that = static_cast<App*>(glfwGetWindowUserPointer(wi));

		glfwGetFramebufferSize(wi, &(that->window_frame_width_), &(that->window_frame_height_));

		for (const auto& v : that->views_)
			v->resize_event(that->window_frame_width_, that->window_frame_height_);
	});

	glfwSetMouseButtonCallback(window_, [] (GLFWwindow* wi, int b, int a, int m)
	{
		App* that = static_cast<App*>(glfwGetWindowUserPointer(wi));
		
		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
		{
			that->inputs_.mouse_buttons_ = 0;
			return;
		}

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGui::GetIO().MousePos = ImVec2(cx, cy);

		for (const auto& v : that->views_)
		{
			if (v->over_viewport(cx, cy))
			{
				if (v.get() != that->focused_)
				{
					that->inputs_.mouse_buttons_ = 0;
					that->focused_ = v.get();
				}

				that->inputs_.shift_pressed_ = (m & GLFW_MOD_SHIFT);
				that->inputs_.control_pressed_ = (m & GLFW_MOD_CONTROL);
				that->inputs_.alt_pressed_ = (m & GLFW_MOD_ALT);
				that->inputs_.meta_pressed_ = (m & GLFW_MOD_SUPER);;
				that->inputs_.last_mouse_x_ = cx;
				that->inputs_.last_mouse_y_ = cy;

				double now = glfwGetTime();

				switch(a)
				{
				case GLFW_PRESS:
					that->inputs_.mouse_buttons_ |= 1 << b;
					v->mouse_press_event(b, that->inputs_.last_mouse_x_, that->inputs_.last_mouse_y_);
					if (now - v->last_click_time() < that->inputs_.double_click_timeout_)
						v->mouse_dbl_click_event(b, that->inputs_.last_mouse_x_, that->inputs_.last_mouse_y_);
					v->set_last_click_time(now);
					break;
				case GLFW_RELEASE:
					that->inputs_.mouse_buttons_ &= ~(1 << b);
					v->mouse_release_event(b, that->inputs_.last_mouse_x_, that->inputs_.last_mouse_y_);
					break;
				}
			}
		}

		if (!that->over_frame(cx, cy))
			that->inputs_.mouse_buttons_ = 0;
	});

	glfwSetCursorPosCallback(window_, [] (GLFWwindow* wi, double x, double y)
	{
		App* that = static_cast<App*>(glfwGetWindowUserPointer(wi));
		
		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
		{
			that->inputs_.mouse_buttons_ = 0;
			return;
		}

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGui::GetIO().MousePos = ImVec2(cx, cy);

		for (const auto& v : that->views_)
		{
			if (that->inputs_.mouse_buttons_ && v->over_viewport(cx, cy))
			{
				if (v.get() != that->focused_)
				{
					that->inputs_.mouse_buttons_ = 0;
					that->focused_ = v.get();
				}
				v->mouse_move_event(x, y);
			}
		}

		if (!that->over_frame(cx, cy))
		{
			that->inputs_.mouse_buttons_ = 0;
			that->focused_ = nullptr;
		}

        that->inputs_.last_mouse_x_ = x;
        that->inputs_.last_mouse_y_ = y;
	});

	glfwSetScrollCallback(window_, [] (GLFWwindow* wi, double dx, double dy)
	{
		App* that = static_cast<App*>(glfwGetWindowUserPointer(wi));

		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsAnyWindowFocused())
		{
			that->inputs_.mouse_buttons_ = 0;
			return;
		}

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);

		for (const auto& v : that->views_)
		{
			if (v->over_viewport(cx, cy))
			{
				if (v.get() != that->focused_)
				{
					that->inputs_.mouse_buttons_ = 0;
					that->focused_ = v.get();
				}
				v->mouse_wheel_event(dx, 100 * dy);
			}
		}
	});

	glfwSetCursorEnterCallback(window_, [] (GLFWwindow* wi, int enter)
	{
		App* that = static_cast<App*>(glfwGetWindowUserPointer(wi));
		
		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
		{
			that->inputs_.mouse_buttons_ = 0;
			return;
		}
		
		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGui::GetIO().MousePos = ImVec2(cx, cy);

		if (enter)
		{
			for (const auto& v : that->views_)
			{
				if (v->over_viewport(cx, cy))
				{
					if (v.get() != that->focused_)
					{
						that->inputs_.mouse_buttons_ = 0;
						that->focused_ = v.get();
					}
				}
			}
			that->inputs_.mouse_buttons_ = 0;
		}
		else
		{
			that->inputs_.mouse_buttons_ = 0;
			that->focused_ = nullptr;
		}
	});

	glfwSetKeyCallback(window_, [] (GLFWwindow* wi, int k, int s, int a, int m)
	{
		App* that = static_cast<App*>(glfwGetWindowUserPointer(wi));

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);

        that->inputs_.shift_pressed_ = (m & GLFW_MOD_SHIFT);
        that->inputs_.control_pressed_ = (m & GLFW_MOD_CONTROL);
        that->inputs_.alt_pressed_ = (m & GLFW_MOD_ALT);
        that->inputs_.meta_pressed_ = (m & GLFW_MOD_SUPER);

		if (k == GLFW_KEY_ESCAPE)
        {
            that->close_event();
			glfwSetWindowShouldClose(that->window_, GLFW_TRUE);
			return;
        }

        switch(a)
        {
            case GLFW_PRESS:
                if (k == GLFW_KEY_SPACE)
                    that->show_imgui_ = !that->show_imgui_;
                else if (k == GLFW_KEY_KP_ADD && that->inputs_.shift_pressed_)
                {
                    that->interface_scaling_ += 0.1f;
                    ImGui::GetIO().FontGlobalScale = that->interface_scaling_;
                }
                else if (k == GLFW_KEY_KP_SUBTRACT && that->inputs_.shift_pressed_)
                {
                    that->interface_scaling_ -= 0.1f;
                    ImGui::GetIO().FontGlobalScale = that->interface_scaling_;
                }
                break;
            case GLFW_RELEASE:
                break;
        }

		for (const auto& v : that->views_)
		{
			if (v->over_viewport(cx, cy))
			{
				that->focused_ = v.get();

				switch(a)
				{
				case GLFW_PRESS:
					if ((k == GLFW_KEY_F) && that->inputs_.control_pressed_  && !that->inputs_.shift_pressed_)
					{
						GLFWmonitor* monitor = glfwGetPrimaryMonitor();
						const GLFWvidmode* mode = glfwGetVideoMode(monitor);
						glfwSetWindowMonitor(wi, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
						glfwSetInputMode(wi, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
						return;
					}
					if ((k == GLFW_KEY_F) && that->inputs_.control_pressed_ && that->inputs_.shift_pressed_)
					{
						int count;
						GLFWmonitor** monitors = glfwGetMonitors(&count);
						if (count > 1)
						{
							const GLFWvidmode* mode = glfwGetVideoMode(monitors[1]);
							glfwSetWindowMonitor(wi, monitors[1], 0, 0, mode->width, mode->height, mode->refreshRate);
							glfwSetInputMode(wi, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
						}
						else
							std::cerr << "Only one monitor" << std::endl;
						return;
					}
					if ((k == GLFW_KEY_W) && that->inputs_.control_pressed_)
					{
						glfwSetWindowMonitor(wi, nullptr, 100, 100, 1024, 1024, 0);
						glfwSetInputMode(wi, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
						return;
					}
					v->key_press_event(k);
					break;
				case GLFW_RELEASE:
					v->key_release_event(k);
					break;
				}
			}
		}
	});

    focused_ = add_view();
}

App::~App()
{}

void App::set_window_size(int32 w, int32 h)
{
	glfwSetWindowSize(window_, w, h);
}

void App::set_window_title(const std::string& name)
{
	window_name_ = name;
	if (window_)
		glfwSetWindowTitle(window_, window_name_.c_str());
}

View* App::add_view()
{
    if (views_.size() < 4)
    {
        glfwMakeContextCurrent(window_);
        views_.push_back(std::make_unique<View>(&inputs_));
        adapt_views_geometry();
        return views_.back().get();
    }
    return nullptr;
}

Module* App::module(const std::string& name) const
{
	auto it = std::find_if(
		modules_.begin(),
		modules_.end(),
		[&] (Module* m) { return m->name().compare(name) == 0; }
	);
	if (it != modules_.end())
		return *it;
	return nullptr;
}

void App::close_event()
{
	for (const auto& v : views_)
		v->close_event();
	
	cgogn::rendering::ShaderProgram::clean_all();
}

bool App::interface()
{
	return false;
}

void App::adapt_views_geometry()
{
	switch(views_.size())
	{
	case 1:
		views_[0]->set_view_ratio(0, 0, 1, 1);
		break;
	case 2:
		views_[0]->set_view_ratio(0, 0, 0.5, 1);
		views_[1]->set_view_ratio(0.5, 0, 0.5, 1);
		break;
	case 3:
		views_[0]->set_view_ratio(0, 0, 0.5, 0.5);
		views_[1]->set_view_ratio(0.5, 0, 0.5, 0.5);
		views_[2]->set_view_ratio(0, 0, 1, 0.5);
		break;
	case 4:
		views_[0]->set_view_ratio(0, 0, 0.5, 0.5);
		views_[1]->set_view_ratio(0.5, 0, 0.5, 0.5);
		views_[2]->set_view_ratio(0, 0.5, 0.5, 0.5);
		views_[3]->set_view_ratio(0.5, 0.5, 0.5, 0.5);
		break;
	}
}

int App::launch()
{
	for (const auto& v : views_)
		v->resize_event(window_frame_width_, window_frame_height_);

	param_frame_ = rendering::ShaderFrame2d::generate_param();
	param_frame_->sz_ = 5.0f;
	param_frame_->color_ = rendering::GLColor(0.25f, 0.25f, 0.25f, 1);

	while (!glfwWindowShouldClose(window_))
	{
		glfwPollEvents();
		glfwMakeContextCurrent(window_);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        for (const auto& v : views_)
        {
            v->draw();
            param_frame_->draw(v->width(), v->height());
        }

        if (show_imgui_)
        {
            ImGui::SetCurrentContext(context_);
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            interface();
            for (Module* m : modules_)
                m->interface();
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        glfwSwapBuffers(window_);
	}
    
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
	
	glfwDestroyWindow(window_);
	return EXIT_SUCCESS;
}

} // namespace cgogn

} // namespace ui
