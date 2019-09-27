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

#ifndef CGOGN_UI_VIEW_H_
#define CGOGN_UI_VIEW_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/ui/gl_viewer.h>
#include <cgogn/ui/module.h>

#include <cgogn/rendering/types.h>
#include <cgogn/rendering/fbo.h>
#include <cgogn/rendering/shaders/shader_fullscreen_texture.h>

namespace cgogn
{

namespace ui
{

class App;

class CGOGN_UI_EXPORT View : public GLViewer
{
	friend class App;

public:

	View(Inputs* inputs, View* share = nullptr);
	virtual ~View() override;

protected:

	virtual void resize_event(int32 frame_width, int32 frame_height) override;
	virtual void close_event() override;

	virtual void mouse_press_event(int32 button, float64 x, float64 y) override;
	virtual void mouse_release_event(int32 button, float64 x, float64 y) override;
	virtual void mouse_dbl_click_event(int32 button, float64 x, float64 y) override;
	virtual void mouse_move_event(float64 x, float64 y) override;
	virtual void mouse_wheel_event(float64 x, float64 y) override;
	virtual void key_press_event(int32 key_code) override;
	virtual void key_release_event(int32 key_code) override;

	void draw();

public:

    void link_module(ViewModule* m);
    void link_module(ProviderModule* m);

	inline bool over_viewport(int32 x, int32 y) const
	{
		y = frame_h_ - y;
		return (x >= viewport_x_) && (x < viewport_x_ + viewport_w_) && (y >= viewport_y_) && (y < viewport_y_ + viewport_h_);
	}

	void set_view_ratio(float64 px, float64 py, float64 pw, float64 ph);

    inline float64 last_click_time() const { return last_click_time_; }
    inline void set_last_click_time(float64 t) { last_click_time_ = t; }

	virtual bool pixel_scene_position(int32 x, int32 y, rendering::GLVec3d& P) const override;
	rendering::GLVec3d unproject(const rendering::GLVec3d& P) const;

protected:

	float64 frame_w_;
	float64 frame_h_;

	float64 viewport_percent_x_;
	float64 viewport_percent_y_;
	float64 viewport_percent_width_;
	float64 viewport_percent_height_;

	float64 last_click_time_;

	std::unique_ptr<rendering::ShaderFSTexture::Param> param_fst_;
	std::unique_ptr<rendering::FBO> fbo_;
	std::unique_ptr<rendering::Texture2D> tex_;

    std::vector<ViewModule*> linked_view_modules_;
    std::vector<ProviderModule*> linked_provider_modules_;
};

} // namespace cgogn

} // namespace ui

#endif // CGOGN_UI_VIEW_H_
