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
#include <cgogn/ui/module.h>

namespace cgogn
{

namespace ui
{

/*****************************************************************************/
// Module
/*****************************************************************************/

Module::Module(const App& app, const std::string& name) : app_(app), name_(name)
{
	app.modules_.push_back(this);
}

Module::~Module()
{
}

void Module::init()
{
}

void Module::main_menu()
{
}

void Module::interface()
{
}

void Module::close_event()
{
}

/*****************************************************************************/
// ViewModule
/*****************************************************************************/

ViewModule::ViewModule(const App& app, const std::string& name) : Module(app, name)
{
}

ViewModule::~ViewModule()
{
}

void ViewModule::mouse_press_event(View*, int32, int32, int32)
{
}
void ViewModule::mouse_release_event(View*, int32, int32, int32)
{
}
void ViewModule::mouse_dbl_click_event(View*, int32, int32, int32)
{
}
void ViewModule::mouse_move_event(View*, int32, int32)
{
}
void ViewModule::mouse_wheel_event(View*, int32, int32)
{
}
void ViewModule::key_press_event(View*, int32)
{
}
void ViewModule::key_release_event(View*, int32)
{
}

void ViewModule::draw(View*)
{
}

/*****************************************************************************/
// ProviderModule
/*****************************************************************************/

ProviderModule::ProviderModule(const App& app, const std::string& name) : Module(app, name)
{
}

ProviderModule::~ProviderModule()
{
}

} // namespace ui

} // namespace cgogn
