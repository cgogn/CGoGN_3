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

#include <cgogn/rendering/imgui_viewer.h>
#include <cgogn/rendering/fbo.h>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

#include <imgui/imgui_internal.h>

#include <iostream>
#include <chrono>
#include <thread>

namespace cgogn
{

namespace rendering
{

static void glfw_error_callback(int error, const char* description)
{
	std::cerr << "Glfw Error " << error << ": " << description << std::endl;
}

/*****************************************************************************/
/*                                 ImGUIViewer                               */
/*****************************************************************************/

ImGUIViewer::ImGUIViewer(ImGUIViewer* share):
	viewport_percent_x_(0),
	viewport_percent_y_(0),
	viewport_percent_width_(1),
	viewport_percent_height_(1),
	last_click_time_(0),
	param_fst_(nullptr),
	fbo_(nullptr),
	tex_(nullptr)
{
	bool err = gl3wInit();
}

ImGUIViewer::~ImGUIViewer()
{}

void ImGUIViewer::set_view_geometry(float64 px, float64 py, float64 pw, float64 ph, float64 frame_width, float64 frame_height)
{
	frame_w_ = frame_width;
	frame_h_ = frame_height;

	viewport_percent_x_ = px;
	viewport_percent_y_ = py;
	viewport_percent_width_ = pw;
	viewport_percent_height_ = ph;

	viewport_x_ = int32(viewport_percent_x_ * frame_width);
	viewport_y_ = int32(viewport_percent_y_ * frame_height);
	viewport_w_ = int32(viewport_percent_width_ * frame_width);
	viewport_h_ = int32(viewport_percent_height_ * frame_height);

	resize_event(viewport_w_, viewport_h_);
}

void ImGUIViewer::update_view_geometry(float64 frame_width, float64 frame_height)
{
	frame_w_ = frame_width;
	frame_h_ = frame_height;

	viewport_x_ = int32(viewport_percent_x_ * frame_width);
	viewport_y_ = int32(viewport_percent_y_ * frame_height);
	viewport_w_ = int32(viewport_percent_width_ * frame_width);
	viewport_h_ = int32(viewport_percent_height_ * frame_height);

	resize_event(viewport_w_, viewport_h_);
}

void ImGUIViewer::resize_event(int32 , int32)
{}

void ImGUIViewer::key_press_event(int32)
{}

void ImGUIViewer::key_release_event(int32)
{}

void ImGUIViewer::close_event()
{}

bool ImGUIViewer::pixel_scene_position(int32 x, int32 y, GLVec3d& P)
{
	float z[4];
	GLint xs, ys;
	float64 xogl;
	float64 yogl;
	float64 zogl;

	if (fbo_ != nullptr)
	{
		xs = GLint(double(x - viewport_x_) / double(viewport_w_) * fbo_->width());
		ys = GLint(double((frame_h_ - y) - viewport_y_) / double(viewport_h_) * fbo_->height());
		fbo_->bind();
		glReadBuffer(GL_DEPTH_ATTACHMENT);
		glReadPixels(xs, ys, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, z);
		fbo_->release();
		if (*z >= 1.0f)
			return false;
		xogl = (float64(xs) / viewport_w_) * 2.0 - 1.0;
		yogl = (float64(ys) / viewport_h_) * 2.0 - 1.0;
		zogl = float64(*z) * 2.0 - 1.0;
	}
	else
	{
		xs = x;
		ys = frame_h_ - y;
		glEnable(GL_DEPTH_TEST);
		glReadBuffer(GL_FRONT);
		glReadPixels(xs, ys, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, z);
		if (*z >= 1.0f)
			return false;
		xogl = (float64(xs - viewport_x_) / viewport_w_) * 2.0 - 1.0;
		yogl = (float64(ys - viewport_y_) / viewport_h_) * 2.0 - 1.0;
		zogl = float64(*z) * 2.0 - 1.0;
	}

	GLVec4d Q(xogl, yogl, zogl, 1.0);
	GLMat4d m = camera().projection_matrix_d() * camera().modelview_matrix_d();
	GLMat4d im = m.inverse();

	GLVec4d P4 = im * Q;
	if (Q.w() != 0.0)
	{
		P.x() = P4.x() / P4.w();
		P.y() = P4.y() / P4.w();
		P.z() = P4.z() / P4.w();
		return true;
	}

	return false;
}

void ImGUIViewer::internal_init()
{
	tex_ = std::make_unique<Texture2D>();
	tex_->alloc(1, 1, GL_RGBA8, GL_RGBA);
	std::vector<Texture2D*> vt{tex_.get()};
	fbo_ = std::make_unique<FBO>(vt, true, nullptr);
	fbo_->resize(width(), height());
	param_fst_ = ShaderFSTexture::generate_param();
	param_fst_->texture_ = fbo_->texture(0);
}

/*****************************************************************************/
/*                                  ImGUIApp                                 */
/*****************************************************************************/

ImGUIApp::ImGUIApp():
	window_(nullptr),
	win_name_("CGoGN"),
	win_frame_width_(512),
	win_frame_height_(512),
	interface_scaling_(1.0f),
	show_imgui_(true),
	focused_(nullptr),
	interface_need_redraw_(true)
{
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
	{
		std::cerr << "Failed to initialize GFLW!" << std::endl;
	}

	// GL 3.3 + GLSL 150 + Core Profile
	const char* glsl_version = "#version 150";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

	window_ = glfwCreateWindow(win_frame_width_, win_frame_height_, win_name_.c_str(), nullptr, nullptr);
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

	int x,y;
	glfwGetWindowSize(window_, &x, &y);
	std::cout << x << " , " << y << std::endl;
	glfwGetFramebufferSize(window_, &x, &y);
	std::cout << x << " , " << y << std::endl;

	glfwSetWindowSizeCallback(window_, [] (GLFWwindow* wi, int, int)
	{
		ImGUIApp* that = static_cast<ImGUIApp*>(glfwGetWindowUserPointer(wi));
		glfwGetFramebufferSize(wi, &(that->win_frame_width_), &(that->win_frame_height_));

		for (ImGUIViewer* v: that->viewers_)
		{
			v->update_view_geometry(that->win_frame_width_, that->win_frame_height_);
			v->camera_.set_aspect_ratio(double(v->viewport_w_) / v->viewport_h_);
			v->need_redraw_ = true;
			v->fbo_->resize(v->viewport_w_, v->viewport_h_);
			v->resize_event(v->viewport_w_, v->viewport_h_);
		}

		that->interface_need_redraw_ = true;
	});

	glfwSetMouseButtonCallback(window_, [] (GLFWwindow* wi, int b, int a, int m)
	{
		ImGUIApp* that = static_cast<ImGUIApp*>(glfwGetWindowUserPointer(wi));
		
		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
		{
			that->interface_need_redraw_ = true;
			that->inputs_.mouse_buttons_ = 0;
			return;
		}

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGui::GetIO().MousePos = ImVec2(cx, cy);

		for (ImGUIViewer* v: that->viewers_)
		{
			if (v->over_viewport(cx, cy))
			{
				if (v != that->focused_)
				{
					that->inputs_.mouse_buttons_ = 0;
					that->focused_ = v;
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
					if (now - v->last_click_time_ < that->inputs_.double_click_timeout_)
						v->mouse_dbl_click_event(b, that->inputs_.last_mouse_x_, that->inputs_.last_mouse_y_);
					v->last_click_time_ = now;
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

	glfwSetScrollCallback(window_, [] (GLFWwindow* wi, double dx, double dy)
	{
		ImGUIApp* that = static_cast<ImGUIApp*>(glfwGetWindowUserPointer(wi));

		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsAnyWindowFocused())
		{
			that->interface_need_redraw_ = true;
			that->inputs_.mouse_buttons_ = 0;
			return;
		}

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);

		for (ImGUIViewer* v: that->viewers_)
		{
			if (v->over_viewport(cx, cy))
			{
				if (v != that->focused_)
				{
					that->inputs_.mouse_buttons_ = 0;
					that->focused_ = v;
				}
				v->mouse_wheel_event(dx, 100 * dy);
			}
		}
	});

	glfwSetCursorPosCallback(window_, [] (GLFWwindow* wi, double x, double y)
	{
		ImGUIApp* that = static_cast<ImGUIApp*>(glfwGetWindowUserPointer(wi));
		
		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
		{
			that->interface_need_redraw_ = true;
			that->inputs_.mouse_buttons_ = 0;
			return;
		}

		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGui::GetIO().MousePos = ImVec2(cx, cy);

		for (ImGUIViewer* v: that->viewers_)
		{
			if (that->inputs_.mouse_buttons_ && v->over_viewport(cx,cy))
			{
				if (v != that->focused_)
				{
					that->inputs_.mouse_buttons_ = 0;
					that->focused_ = v;
				}
				v->mouse_move_event(x, y);
			}
		}

		if (!that->over_frame(cx, cy))
		{
			that->inputs_.mouse_buttons_ = 0;
			that->focused_ = nullptr;
		}
	});

	glfwSetCursorEnterCallback(window_, [] (GLFWwindow* wi, int enter)
	{
		ImGUIApp* that = static_cast<ImGUIApp*>(glfwGetWindowUserPointer(wi));
		
		if (ImGui::GetIO().WantCaptureMouse || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
		{
			that->interface_need_redraw_ = true;
			that->inputs_.mouse_buttons_ = 0;
			return;
		}
		
		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGui::GetIO().MousePos = ImVec2(cx, cy);

		if (enter)
		{
			for (ImGUIViewer* v: that->viewers_)
			{
				if (v->over_viewport(cx, cy))
				{
					if (v != that->focused_)
					{
						that->inputs_.mouse_buttons_ = 0;
						that->focused_ = v;
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
		double cx, cy;
		glfwGetCursorPos(wi, &cx, &cy);
		ImGUIApp* that = static_cast<ImGUIApp*>(glfwGetWindowUserPointer(wi));

		if (k == GLFW_KEY_ESCAPE)
			exit(0);

		for (ImGUIViewer* v: that->viewers_)
		{
			if (v->over_viewport(cx, cy))
			{
				that->focused_ = v;
				that->inputs_.shift_pressed_ = (m & GLFW_MOD_SHIFT);
				that->inputs_.control_pressed_ = (m & GLFW_MOD_CONTROL);
				that->inputs_.alt_pressed_ = (m & GLFW_MOD_ALT);
				that->inputs_.meta_pressed_ = (m & GLFW_MOD_SUPER);

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

		switch(a)
		{
			case GLFW_PRESS:
				that->key_press_event(k);
			break;
			case GLFW_RELEASE:
				that->key_release_event(k);
			break;
		}
	});
}

ImGUIApp::~ImGUIApp()
{}

void ImGUIApp::close_event()
{
	for (ImGUIViewer* v: viewers_)
		v->close_event();
}

bool ImGUIApp::interface()
{
	return true;
}

void ImGUIApp::key_press_event(int32)
{}

void ImGUIApp::key_release_event(int32)
{}

void ImGUIApp::set_window_size(int32 w, int32 h)
{
	glfwSetWindowSize(window_, w, h);
}

void ImGUIApp::set_window_title(const std::string& name)
{
	win_name_ = name;
	if (window_)
		glfwSetWindowTitle(window_, win_name_.c_str());
}

void ImGUIApp::add_view(ImGUIViewer* view)
{
	glfwMakeContextCurrent(window_);
	view->internal_init();
	view->init();
	view->set_inputs(&inputs_);
	viewers_.push_back(view);
	focused_ = view;
}

void ImGUIApp::adapt_views_geometry()
{
	switch(viewers_.size())
	{
	case 1:
		viewers_[0]->set_view_geometry(0, 0, 1, 1, win_frame_width_, win_frame_height_);
		break;
	case 2:
		viewers_[0]->set_view_geometry(0, 0, 0.5, 1, win_frame_width_, win_frame_height_);
		viewers_[1]->set_view_geometry(0.5, 0, 0.5, 1, win_frame_width_, win_frame_height_);
		break;
	case 3:
		viewers_[0]->set_view_geometry(0, 0, 0.5, 0.5, win_frame_width_, win_frame_height_);
		viewers_[1]->set_view_geometry(0.5, 0, 0.5, 0.5, win_frame_width_, win_frame_height_);
		viewers_[2]->set_view_geometry(0, 0, 1, 0.5, win_frame_width_, win_frame_height_);
		break;
	case 4:
		viewers_[0]->set_view_geometry(0, 0, 0.5, 0.5, win_frame_width_, win_frame_height_);
		viewers_[1]->set_view_geometry(0.5, 0, 0.5, 0.5, win_frame_width_, win_frame_height_);
		viewers_[2]->set_view_geometry(0, 0.5, 0.5, 0.5, win_frame_width_, win_frame_height_);
		viewers_[3]->set_view_geometry(0.5, 0.5, 0.5, 0.5, win_frame_width_, win_frame_height_);
		break;
	}
}

int ImGUIApp::launch()
{
	interface_need_redraw_ = true;
	adapt_views_geometry();

	for (ImGUIViewer* v : viewers_)
	{
		v->update_view_geometry(win_frame_width_, win_frame_height_);
		v->fbo_->resize(v->viewport_w_, v->viewport_h_);
		v->need_redraw_ = true;
		v->camera_.set_aspect_ratio(double(v->viewport_w_) / v->viewport_h_);
		v->resize_event(v->viewport_w_, v->viewport_h_);
	}

	param_frame_ = ShaderFrame2d::generate_param();
	param_frame_->sz_ = 5.0f;
	param_frame_->color_ = GLColor(0.5f, 0.5f, 0.5f, 1);

	while (!glfwWindowShouldClose(window_))
	{
		glfwPollEvents();
		glfwMakeContextCurrent(window_);

		for (ImGUIViewer* v : viewers_)
			interface_need_redraw_ |= v->need_redraw_;

		if (interface_need_redraw_ )
		{
			bool interface_need_redraw_now = false;
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			for (ImGUIViewer* v : viewers_)
			{
				v->spin();
				v->set_viewport();
				if (v->need_redraw_)
				{
					v->fbo_->bind();
					glEnable(GL_DEPTH_TEST);
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
					GLenum idbuf = GL_COLOR_ATTACHMENT0;
					glDrawBuffers(1, &idbuf);
					v->draw();
				}

				v->fbo_->release();
				glDisable(GL_DEPTH_TEST);
				v->param_fst_->draw();
				param_frame_->draw(v->width(), v->height());

				v->need_redraw_ = v->camera().is_moving_;
			}

			if (show_imgui_ && interface_need_redraw_)
			{
				ImGui_ImplOpenGL3_NewFrame();
				ImGui_ImplGlfw_NewFrame();
				ImGui::NewFrame();
				interface_need_redraw_now = interface();
				ImGui::Render();
				ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
			}

			glfwSwapBuffers(window_);
			interface_need_redraw_ = interface_need_redraw_now;
		}
		else
			std::this_thread::sleep_for(std::chrono::milliseconds(16));
	}
	
	glfwDestroyWindow(window_);
	return EXIT_SUCCESS;
}

bool ImGUIApp::IsItemActiveLastFrame()
{
	if (context_->ActiveIdPreviousFrame)
		return context_->ActiveIdPreviousFrame == context_->CurrentWindow->DC.LastItemId;
	return false;
}

} // namespace cgogn

} // namespace rendering
