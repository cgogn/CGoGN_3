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

#ifndef CGOGN_UI_VIEW_H_
#define CGOGN_UI_VIEW_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/ui/gl_viewer.h>
#include <cgogn/ui/module.h>

#include <cgogn/rendering/fbo.h>
#include <cgogn/rendering/shaders/shader_fullscreen_texture.h>
#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace ui
{

class App;

class CGOGN_UI_EXPORT View : public GLViewer
{
	friend class App;

public:
	View(Inputs* inputs, const std::string& name);
	virtual ~View() override;

protected:
	void resize_event(int32 window_width, int32 window_height, int32 frame_buffer_width, int32 frame_buffer_height);
	virtual void close_event() override;

	virtual void mouse_press_event(int32 button, int32 x, int32 y) override;
	virtual void mouse_release_event(int32 button, int32 x, int32 y) override;
	virtual void mouse_dbl_click_event(int32 button, int32 x, int32 y) override;
	virtual void mouse_move_event(int32 x, int32 y) override;
	virtual void mouse_wheel_event(float64 dx, float64 dy) override;
	virtual void key_press_event(int32 key_code) override;
	virtual void key_release_event(int32 key_code) override;

	void draw();

public:
	const std::string& name() const
	{
		return name_;
	}

	void link_module(ViewModule* m);
	void link_module(ProviderModule* m);

	inline bool contains(int32 x, int32 y) const
	{
		return x >= x_offset_ && x < x_offset_ + width_ && y >= y_offset_ && y < y_offset_ + height_;
	}

	void set_view_ratio(float64 px, float64 py, float64 pw, float64 ph);

	void update_scene_bb();
	void lock_scene_bb();
	void unlock_scene_bb();

	virtual bool pixel_scene_position(int32 x, int32 y, rendering::GLVec3d& P) const override;
	virtual std::pair<rendering::GLVec3d, rendering::GLVec3d> pixel_ray(int32 x, int32 y) const override;
	rendering::GLVec3d unproject(int32 x, int32 y, float64 z) const;

	inline void stop_event()
	{
		event_stopped_ = true;
	}

	bool scene_bb_locked_;

protected:
	std::string name_;

	float64 ratio_x_offset_;
	float64 ratio_y_offset_;
	float64 ratio_width_;
	float64 ratio_height_;

	int32 x_offset_;
	int32 y_offset_;
	int32 width_;
	int32 height_;

	// viewport_width_ & viewport_heght_ are in GLViewer
	int32 viewport_x_offset_;
	int32 viewport_y_offset_;

	std::unique_ptr<rendering::ShaderFSTexture::Param> param_fst_;
	std::unique_ptr<rendering::FBO> fbo_;
	std::unique_ptr<rendering::Texture2D> tex_;

	std::vector<ViewModule*> linked_view_modules_;
	std::vector<ProviderModule*> linked_provider_modules_;

	bool event_stopped_;
	bool closing_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_UI_VIEW_H_
