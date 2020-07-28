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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_DISTANCE_H_
#define CGOGN_GEOMETRY_FUNCTIONS_DISTANCE_H_

#include <cgogn/core/utils/assert.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

/**
 * @brief distance plane point
 * @param normal normal of the plane
 * @param d distance to the origin of the plane
 * @param p point to compute distance to plane
 * @return distance
 */
inline Scalar distance_plane_point(const Vec3& normal, Scalar d, const Vec3& p)
{
	return abs(normal.dot(p) + d);
}

/**
 * @brief squared distance line point (optimized version for testing many points with the same line)
 * @param A one point of line
 * @param AB normalized vector or line
 * @param P point to compute distance to line
 * @return distance
 */
inline Scalar squared_distance_normalized_line_point(const Vec3& A, const Vec3& AB_norm, const Vec3& P)
{
	return ((A - P).cross(AB_norm)).squaredNorm();
}

/**
 * @brief squared distance line point
 * @param A one point of line
 * @param B second point of line
 * @param P point to compute distance to line
 * @return distance
 */
inline Scalar squared_distance_line_point(const Vec3& A, const Vec3& B, const Vec3& P)
{
	Vec3 AB = B - A;
	cgogn_message_assert(AB.squaredNorm() > 0.0, "line must be defined by 2 different points");
	AB.normalize();
	return squared_distance_normalized_line_point(A, AB, P);
}

/**
 * compute squared distance from line to segment
 * @param A point of line
 * @param AB vector of line
 * @param AB2 AB*AB (for optimization if call several times with AB
 * @param P first point of segment
 * @param Q second point of segment
 * @return the squared distance
 */
inline Scalar squared_distance_line_seg(const Vec3& A, const Vec3& AB, Scalar AB2, const Vec3& P, const Vec3& Q)
{
	Vec3 PQ = Q - P;
	Scalar PQ_n2 = PQ.squaredNorm();

	// if P == Q compute distance to P
	if (PQ_n2 == Scalar(0)) // P == Q
	{
		Vec3 V = AB / AB.norm();
		return squared_distance_normalized_line_point(A, V, P);
	}

	Scalar X = AB.dot(PQ);
	Vec3 AP = P - A;

	Scalar beta = (AB2 * (AP.dot(PQ)) - X * (AP.dot(AB))) / (X * X - AB2 * PQ_n2);

	if (beta < Scalar(0))
	{
		Vec3 W = AB.cross(AP);
		return W.squaredNorm() / AB2;
	}

	if (beta > Scalar(1))
	{
		Vec3 AQ = Q - A;
		Vec3 W = AB.cross(AQ);
		return W.squaredNorm() / AB2;
	}

	Vec3 temp = AB.cross(PQ);
	Scalar num = AP.dot(temp);
	Scalar den = temp.squaredNorm();

	return (num * num) / den;
}

/**
 * compute squared distance from line to segment
 * @warning if used many times with same line prefer version, with A, AB and AB2 parameter
 * @param A point of line
 * @param B point of line
 * @param P first point of segment
 * @param Q second point of segment
 * @return the squared distance
 */
inline Scalar squared_distance_line_seg(const Vec3& A, const Vec3& B, const Vec3& P, const Vec3& Q)
{
	Vec3 AB = B - A;
	return squared_distance_line_seg(A, AB, AB.dot(AB), P, Q);
}

/**
 * compute squared distance from segment to point
 * @param A point of segment
 * @param AB vector of segment
 * @param P the point
 * @return the squared distance
 */
inline Scalar squared_distance_seg_point(const Vec3& A, const Vec3& AB, const Vec3& P)
{
	Vec3 AP = P - A;

	// squared vector length
	Scalar AB2 = AB.dot(AB);

	// position of projection of P on [A,B]
	Scalar t = AP.dot(AB) / AB2;

	// before A, distance is PA
	if (t <= Scalar(0.))
		return AP.squaredNorm();

	// after B, distance is PB
	if (t >= Scalar(1.))
	{
		Vec3 BP = P - (AB + A);
		return BP.squaredNorm();
	}

	// between A & B, distance is projection on (AB)
	Vec3 X = AB.cross(AP);
	return X.squaredNorm() / AB2;
}

/**
 * compute squared distance from point to triangle
 * @param P the point
 * @param A triangle point 1
 * @param B triangle point 2
 * @param C triangle point 3
 * @return the squared distance
 */
Scalar squared_distance_point_triangle(const Vec3& P, const Vec3& A, const Vec3& B, const Vec3& C);

/**
 * compute the barycentric coordinates of the point in the triangle ABC that is closest to point P
 * @param P the point
 * @param A triangle point 1
 * @param B triangle point 2
 * @param C triangle point 3
 * @param u barycentric coordinate 1 of closest point
 * @param v barycentric coordinate 2 of closest point
 * @param w barycentric coordinate 3 of closest point
 */
void closest_point_in_triangle(const Vec3& P, const Vec3& A, const Vec3& B, const Vec3& C, Scalar& u, Scalar& v,
							   Scalar& w);

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_DISTANCE_H_
