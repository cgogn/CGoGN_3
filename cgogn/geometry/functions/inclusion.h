/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_INCLUSION_H_
#define CGOGN_GEOMETRY_FUNCTIONS_INCLUSION_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/normal.h>

namespace cgogn
{

namespace geometry
{

inline bool in_sphere(const Vec3& point, const Vec3& center, Scalar radius)
{
	return (point - center).norm() < radius;
}

inline double triple_product(const Vec3& U, const Vec3& V, const Vec3& W)
{
	return U.dot(V.cross(W));
}

inline bool in_triangle(const Vec3& P, const Vec3& normal, const Vec3& Ta, const Vec3& Tb, const Vec3& Tc)
{
	if (triple_product(P - Ta, Tb - Ta, normal) >= 0 || triple_product(P - Tb, Tc - Tb, normal) >= 0 ||
		triple_product(P - Tc, Ta - Tc, normal) >= 0)
		return false;

	return true;
}

inline bool in_triangle(const Vec3& P, const Vec3& Ta, const Vec3& Tb, const Vec3& Tc)
{
	return in_triangle(P, normal(Ta, Tb, Tc), Ta, Tb, Tc);
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_INCLUSION_H_
