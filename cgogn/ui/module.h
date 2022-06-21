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

#ifndef CGOGN_UI_MODULE_H_
#define CGOGN_UI_MODULE_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

namespace cgogn
{

namespace ui
{

class App;
class View;

class CGOGN_UI_EXPORT Module
{
	friend class App;

public:
	Module(const App& app, const std::string& name);
	virtual ~Module();

	const std::string& name() const
	{
		return name_;
	}

protected:
	virtual void init();
	virtual void main_menu();
	virtual void ui_interface();
	virtual void popups();

	virtual void close_event();

	const App& app_;
	std::string name_;
};

/*****************************************************************************/
// ViewModule
/*****************************************************************************/

class CGOGN_UI_EXPORT ViewModule : public Module
{
	friend class View;

public:
	ViewModule(const App& app, const std::string& name);
	virtual ~ViewModule();

	const std::vector<View*>& linked_views() const
	{
		return linked_views_;
	}

protected:
	virtual void mouse_press_event(View* view, int32 button, int32 x, int32 y);
	virtual void mouse_release_event(View* view, int32 button, int32 x, int32 y);
	virtual void mouse_dbl_click_event(View* view, int32 button, int32 x, int32 y);
	virtual void mouse_move_event(View* view, int32 x, int32 y);
	virtual void mouse_wheel_event(View* view, int32 dx, int32 dy);
	virtual void key_press_event(View* view, int32 key_code);
	virtual void key_release_event(View* view, int32 key_code);

	virtual void draw(View* view);

	std::vector<View*> linked_views_;
};

/*****************************************************************************/
// ProviderModule
/*****************************************************************************/

class CGOGN_UI_EXPORT ProviderModule : public Module
{
	friend class View;

public:
	ProviderModule(const App& app, const std::string& name);
	virtual ~ProviderModule();

	virtual std::pair<geometry::Vec3, geometry::Vec3> meshes_bb() const = 0;

protected:
	std::vector<View*> linked_views_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_UI_MODULE_H_
