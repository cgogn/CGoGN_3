
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

#ifndef CGOGN_RENDERING_MFRAME_H_
#define CGOGN_RENDERING_MFRAME_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace rendering
{

struct CGOGN_RENDERING_EXPORT MovingFrame
{
	Transfo3d frame_;
	Transfo3d spin_;
	bool is_moving_;

	MovingFrame():
		frame_(Transfo3d::Identity()),
		spin_(Transfo3d::Identity()),
		is_moving_(false)
	{}

	//GLVec3d local_coordinates(GLVec3d glob);
};

} // namespace cgogn

} // namespace rendering

#endif // CGOGN_RENDERING_MFRAME_H_
