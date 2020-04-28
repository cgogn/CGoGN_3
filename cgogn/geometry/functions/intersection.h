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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_INTERSECTION_H_
#define CGOGN_GEOMETRY_FUNCTIONS_INTERSECTION_H_

#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/functions/inclusion.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cmath>

namespace cgogn
{

namespace geometry
{

enum Intersection
{
	NO_INTERSECTION = 0,
	VERTEX_INTERSECTION,
	EDGE_INTERSECTION,
	FACE_INTERSECTION
};

inline bool intersection_ray_triangle(const Vec3& P, const Vec3& Dir, const Vec3& Ta, const Vec3& Tb, const Vec3& Tc,
									  Vec3* inter = nullptr)
{
	Vec3 u = Ta - P;
	Vec3 v = Tb - P;
	Vec3 w = Tc - P;

	Scalar x = Dir.dot(u.cross(v));
	Scalar y = Dir.dot(v.cross(w));
	Scalar z = Dir.dot(w.cross(u));

	uint32 np = 0;
	uint32 nn = 0;
	uint32 nz = 0;

	if (x > Scalar(0))
		++np;
	else if (x < Scalar(0))
		++nn;
	else
		++nz;

	if (y > Scalar(0))
		++np;
	else if (y < Scalar(0))
		++nn;
	else
		++nz;

	if (z > Scalar(0))
		++np;
	else if (z < Scalar(0))
		++nn;
	else
		++nz;

	// line intersect the triangle
	if (((np != 0) && (nn != 0)) || (nz == 3))
		return false;

	Scalar sum = x + y + z;
	Scalar alpha = y / sum;
	Scalar beta = z / sum;
	Scalar gamma = Scalar(1) - alpha - beta;
	Vec3 I = Ta * alpha + Tb * beta + Tc * gamma;

	// it's a ray not a line !
	if (Dir.dot(I - P) < 0.0)
		return false;

	if (inter)
		*inter = I;

	return true;
}

/**
 * \param[in] center the position of the center of the sphere.
 * \param[in] radius the radius of the sphere
 * \param[in] p1 first point of the segment
 * \param[in] p2 second point of the segment
 * \param[out] alpha ratio of the segment inside the sphere
 */
inline bool intersection_sphere_segment(const Vec3& center, Scalar radius, const Vec3& p1, const Vec3& p2,
										Scalar& alpha)
{
	if (in_sphere(p1, center, radius) && !in_sphere(p2, center, radius))
	{
		Vec3 p = p1 - center;
		Vec3 qminusp = p2 - center - p;
		Scalar s = p.dot(qminusp);
		Scalar n2 = qminusp.squaredNorm();
		alpha = (-s + std::sqrt(s * s + n2 * (radius * radius - p.squaredNorm()))) / n2;
		return true;
	}

	return false;
}

inline Intersection intersection_segment_segment(const Vec3& PA, const Vec3& PB, const Vec3& QA, const Vec3& QB,
												 Vec3& Inter)
{
	Vec3 vp1p2 = PB - PA;
	Vec3 vq1q2 = QB - QA;
	Vec3 vp1q1 = QA - PA;

	Scalar delta = vp1p2[0] * vq1q2[1] - vp1p2[1] * vq1q2[0];
	Scalar coeff = vp1q1[0] * vq1q2[1] - vp1q1[1] * vq1q2[0];

	if (delta == 0) // parallel
	{
		// test if colinear
		if (coeff == 0)
		{
			// colinear
			// TODO : check if there is a common point between the two edges
			Inter = QA;
			return EDGE_INTERSECTION;
		}
		else
			return NO_INTERSECTION;
	}
	else
		Inter = Vec3((PA[0] * delta + vp1p2[0] * coeff) / delta, (PA[1] * delta + vp1p2[1] * coeff) / delta,
					 (PA[2] * delta + vp1p2[2] * coeff) / delta);

	// test if inter point is outside the edges
	if ((Inter[0] < PA[0] && Inter[0] < PB[0]) || (Inter[0] > PA[0] && Inter[0] > PB[0]) ||
		(Inter[0] < QA[0] && Inter[0] < QB[0]) || (Inter[0] > QA[0] && Inter[0] > QB[0]) ||
		(Inter[1] < PA[1] && Inter[1] < PB[1]) || (Inter[1] > PA[1] && Inter[1] > PB[1]) ||
		(Inter[1] < QA[1] && Inter[1] < QB[1]) || (Inter[1] > QA[1] && Inter[1] > QB[1]))
		return NO_INTERSECTION;

	if (PA.isApprox(Inter) || PB.isApprox(Inter) || QA.isApprox(Inter) || QB.isApprox(Inter))
		return VERTEX_INTERSECTION;

	return EDGE_INTERSECTION;
}

inline bool intersection_line_plane(const Vec3& point_line, const Vec3& dir_line, const Vec3& point_plane,
							 const Vec3& normal_plane, Vec3* inter = nullptr)
{
	const Scalar PRECISION = std::numeric_limits<Scalar>::epsilon();

	Scalar b = normal_plane.dot(dir_line);

	if (std::abs(b) < PRECISION)
		return false;

	Scalar a = normal_plane.dot(point_plane - point_line);
	if (inter)
		*inter = point_line + (a / b) * dir_line;

	return true;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_INTERSECTION_H_
