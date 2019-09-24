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

#include <cgogn/ui/view.h>

namespace cgogn
{

namespace ui
{

View::View(Inputs* inputs, View* share) :
	GLViewer(inputs),
    frame_h_(512),
    frame_w_(512),
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

	tex_ = std::make_unique<rendering::Texture2D>();
	tex_->alloc(1, 1, GL_RGBA8, GL_RGBA);
	std::vector<rendering::Texture2D*> vt{tex_.get()};
	fbo_ = std::make_unique<rendering::FBO>(vt, true, nullptr);
	fbo_->resize(width(), height());
	param_fst_ = rendering::ShaderFSTexture::generate_param();
	param_fst_->texture_ = fbo_->texture(0);
}

View::~View()
{}

void View::set_view_ratio(float64 px, float64 py, float64 pw, float64 ph)
{
	viewport_percent_x_ = px;
	viewport_percent_y_ = py;
	viewport_percent_width_ = pw;
	viewport_percent_height_ = ph;
}

void View::resize_event(int32 frame_width, int32 frame_height)
{
	frame_w_ = frame_width;
	frame_h_ = frame_height;

	viewport_x_ = int32(viewport_percent_x_ * frame_width);
	viewport_y_ = int32(viewport_percent_y_ * frame_height);
	viewport_w_ = int32(viewport_percent_width_ * frame_width);
	viewport_h_ = int32(viewport_percent_height_ * frame_height);

	fbo_->resize(viewport_w_, viewport_h_);

	for (Module* m : linked_modules_)
		m->resize_event(this, viewport_w_, viewport_h_);

	GLViewer::resize_event(viewport_w_, viewport_h_);
}

void View::close_event()
{
	for (Module* m : linked_modules_)
		m->close_event();
}

void View::mouse_press_event(int32 button, float64 x, float64 y)
{
	for (Module* m : linked_modules_)
		m->mouse_press_event(this, button, x, y);

	GLViewer::mouse_press_event(button, x, y);
}

void View::mouse_release_event(int32 button, float64 x, float64 y)
{
	for (Module* m : linked_modules_)
		m->mouse_release_event(this, button, x, y);

	GLViewer::mouse_release_event(button, x, y);
}

void View::mouse_dbl_click_event(int32 button, float64 x, float64 y)
{
	for (Module* m : linked_modules_)
		m->mouse_dbl_click_event(this, button, x, y);
	
	GLViewer::mouse_dbl_click_event(button, x, y);
}

void View::mouse_move_event(float64 x, float64 y)
{
	for (Module* m : linked_modules_)
		m->mouse_move_event(this, x, y);
	
	GLViewer::mouse_move_event(x, y);
}

void View::mouse_wheel_event(float64 x, float64 y)
{
	for (Module* m : linked_modules_)
		m->mouse_wheel_event(this, x, y);
	
	GLViewer::mouse_wheel_event(x, y);
}

void View::key_press_event(int32 key_code)
{
	for (Module* m : linked_modules_)
		m->key_press_event(this, key_code);
	
	GLViewer::key_press_event(key_code);
}

void View::key_release_event(int32 key_code)
{
	for (Module* m : linked_modules_)
		m->key_release_event(this, key_code);
	
	GLViewer::key_release_event(key_code);
}

void View::draw()
{
	spin();
	glViewport(viewport_x_, viewport_y_, viewport_w_, viewport_h_);
	
	if (need_redraw_)
	{
		fbo_->bind();
		glEnable(GL_DEPTH_TEST);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		GLenum idbuf = GL_COLOR_ATTACHMENT0;
		glDrawBuffers(1, &idbuf);
		for (Module* m : linked_modules_)
			m->draw(this);
		fbo_->release();
		glDisable(GL_DEPTH_TEST);
	}

	param_fst_->draw();
}

void View::link_module(Module* m)
{
	linked_modules_.push_back(m);
	m->linked_views_.push_back(this);
}

bool View::pixel_scene_position(int32 x, int32 y, rendering::GLVec3d& P) const
{
	float z[4];
	GLint xs, ys;
	float64 xogl;
	float64 yogl;
	float64 zogl;

	if (fbo_)
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

	rendering::GLVec4d Q(xogl, yogl, zogl, 1.0);
	rendering::GLMat4d im = (camera().projection_matrix_d() * camera().modelview_matrix_d()).inverse();
	rendering::GLVec4d P4 = im * Q;
	if (P4.w() != 0.0)
	{
		P.x() = P4.x() / P4.w();
		P.y() = P4.y() / P4.w();
		P.z() = P4.z() / P4.w();
		return true;
	}

	return false;
}

rendering::GLVec3d View::unproject(const rendering::GLVec3d& P) const
{
	float64 xogl = ((P.x() - viewport_x_) / viewport_w_) * 2.0 - 1.0;
	float64 yogl = ((P.y() - viewport_y_) / viewport_h_) * 2.0 - 1.0;
	float64 zogl = P.z() * 2.0 - 1.0;
	rendering::GLVec4d Q(xogl, yogl, zogl, 1.0);
	rendering::GLMat4d im = (camera().projection_matrix_d() * camera().modelview_matrix_d()).inverse();
	rendering::GLVec4d res = im * Q;
	res /= res.w();
	return res.head(3);
}

} // namespace cgogn

} // namespace ui
