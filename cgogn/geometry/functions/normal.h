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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_NORMAL_H_
#define CGOGN_GEOMETRY_FUNCTIONS_NORMAL_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/vector_ops.h>

namespace cgogn
{

namespace geometry
{

/**
 * normal of the plane spanned by 3 points in 3D
 */
template <typename VEC3>
VEC3 normal(const VEC3& p1, const VEC3& p2, const VEC3& p3)
{
	static_assert(vector_traits<VEC3>::SIZE == 3, "The dimension of the vector must be equal to 3.");
	VEC3 v1 = p2 - p1;
	VEC3 v2 = p3 - p1;
	return cross(v1, v2);
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_NORMAL_H_
