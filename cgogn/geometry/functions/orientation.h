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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_ORIENTATION_H_
#define CGOGN_GEOMETRY_FUNCTIONS_ORIENTATION_H_

#include <cgogn/geometry/functions/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

enum Orientation2D
{
	ALIGNED = 0,
	RIGHT,
	LEFT
};

enum Orientation3D
{
	ON = 0,
	OVER,
	UNDER
};

/**
 * get the orientation of point P w.r.t. the plane spanned by the 3 given points
 * @param P the point
 * @param A plane point 1
 * @param B plane point 2
 * @param C plane point 3
 * @return the orientation
 */
inline Orientation3D test_orientation_3D(const Vec3& P, const Vec3& A, const Vec3& B, const Vec3& C)
{
	// plane
	Vec3 n = normal(A, B, C);
	n.normalize();
	Scalar d = -(A.dot(n));

	const Scalar dist = n.dot(P) + d;

	if (cgogn::almost_equal_relative(dist, Scalar(0)))
		return Orientation3D::ON;

	if (dist < -Scalar(0))
		return Orientation3D::UNDER;

	return Orientation3D::OVER;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_ORIENTATION_H_
