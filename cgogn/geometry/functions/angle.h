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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_ANGLE_H_
#define CGOGN_GEOMETRY_FUNCTIONS_ANGLE_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>

namespace cgogn
{

namespace geometry
{

inline Scalar cos_angle(const Vec3& a, const Vec3& b)
{
    Scalar ab = std::sqrt(a.squaredNorm() * b.squaredNorm());
	Scalar res = ab > 1e-20 ? a.dot(b) / ab : Scalar(1);
	return std::clamp(res, Scalar(-1), Scalar(1));
}

inline Scalar angle(const Vec3& a, const Vec3& b)
{
	return std::acos(cos_angle(a, b));
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_ANGLE_H_
