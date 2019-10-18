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

#ifndef CGOGN_UI_MOVING_FRAME_H_
#define CGOGN_UI_MOVING_FRAME_H_

#include <cgogn/ui/cgogn_ui_export.h>

#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace ui
{

struct CGOGN_UI_EXPORT MovingFrame
{
	rendering::Transfo3d frame_;
	rendering::Transfo3d spin_;
	bool is_moving_;

	MovingFrame():
		frame_(rendering::Transfo3d::Identity()),
		spin_(rendering::Transfo3d::Identity()),
		is_moving_(false)
	{}

	// GLVec3d local_coordinates(GLVec3d glob);
};

} // namespace cgogn

} // namespace ui

#endif // CGOGN_UI_MOVING_FRAME_H_
