
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

#include <cgogn/rendering_pureGL/camera.h>

namespace cgogn
{
namespace rendering_pgl
{
GLMat4d Camera::perspective(float64 znear, float64 zfar) const
{
	float64 range_inv = 1.0 / (znear - zfar);
	float64 f = 1.0/std::tan(field_of_view_/2.0);
	auto m05 = (asp_ratio_>1) ? std::make_pair(f/asp_ratio_,f) : std::make_pair(f,f*asp_ratio_);
	GLMat4d m;
	m << m05.first,  0,  0,  0,
		  0, m05.second,  0,  0,
		  0,  0, (znear+zfar)*range_inv, 2*znear*zfar*range_inv,
		  0,  0, -1 ,0;
	return m;
}

GLMat4d Camera::ortho(float64 znear, float64 zfar) const
{
	float64 range_inv = 1.0 / (znear - zfar);
	auto m05 = (asp_ratio_<1) ? std::make_pair(1.0/asp_ratio_,1.0) : std::make_pair(1.0,1.0/asp_ratio_);
	GLMat4d m;
	m << m05.first,  0,  0,  0,
		  0, m05.second,  0,  0,
		  0,  0, 2*range_inv, 0,
		  0,  0, (znear+zfar)*range_inv,0;
	return m;
}


}
}
