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

#ifndef CGOGN_UI_INPUTS_H_
#define CGOGN_UI_INPUTS_H_

#include <cgogn/ui/cgogn_ui_export.h>
#include <cgogn/core/utils/numerics.h>

namespace cgogn
{

namespace ui
{

struct CGOGN_UI_EXPORT Inputs
{
    Inputs() :
        wheel_sensitivity_(0.0025),
        mouse_sensitivity_(0.005),
        spin_sensitivity_(0.025),
        double_click_timeout_(0.3),
        need_redraw_(true),
        shift_pressed_(false),
        control_pressed_(false),
        alt_pressed_(false),
        meta_pressed_(false)
    {}

	float64 wheel_sensitivity_;
	float64 mouse_sensitivity_;
	float64 spin_sensitivity_;
	float64 double_click_timeout_;

	int32 last_mouse_x_;
	int32 last_mouse_y_;
	uint32 mouse_buttons_;

	bool need_redraw_;
	bool shift_pressed_;
	bool control_pressed_;
	bool alt_pressed_;
	bool meta_pressed_;
};

} // namespace cgogn

} // namespace ui

#endif // CGOGN_UI_INPUTS_H_
