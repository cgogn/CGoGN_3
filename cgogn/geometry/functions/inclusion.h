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

namespace cgogn
{

namespace geometry
{

inline bool in_sphere(const Vec3& point, const Vec3& center, Scalar radius)
{
	return (point - center).norm() < radius;
}






template <typename VEC3a, typename VEC3b, typename VEC3c>
double triple_product(const Eigen::MatrixBase<VEC3a>& U, const Eigen::MatrixBase<VEC3b>& V, const Eigen::MatrixBase<VEC3c>& W)
{
//	static_assert(is_same_vector<VEC3a,VEC3b,VEC3c>::value, "parameters must have same type");
//	static_assert(is_dim_of<VEC3a, 3>::value, "The size of the vector must be equal to 3.");

	return U.dot(V.cross(W));
}


template <typename VEC3a, typename VEC3b, typename VEC3c, typename VEC3d, typename VEC3e>
bool in_triangle(const Eigen::MatrixBase<VEC3a>& P, const Eigen::MatrixBase<VEC3b>& normal,
				 const Eigen::MatrixBase<VEC3c>& Ta,  const Eigen::MatrixBase<VEC3d>& Tb, const VEC3e& Tc)
{
//	static_assert(is_same_vector<VEC3a,VEC3b,VEC3c,VEC3d,VEC3e>::value, "parameters must have same type");
//	static_assert(is_dim_of<VEC3a, 3>::value, "The size of the vector must be equal to 3.");

	if (triple_product(P-Ta, Tb-Ta, normal) >= 0 ||
		triple_product(P-Tb, Tc-Tb, normal) >= 0  ||
		triple_product(P-Tc, Ta-Tc, normal) >= 0  )
		return false;

	return true;
}



template <typename VEC3a, typename VEC3b, typename VEC3c, typename VEC3d>
bool in_triangle(const Eigen::MatrixBase<VEC3a>& P, const Eigen::MatrixBase<VEC3b>& Ta,  const Eigen::MatrixBase<VEC3c>& Tb, const Eigen::MatrixBase<VEC3d>& Tc)
{
//	static_assert(is_same_vector<VEC3a,VEC3b,VEC3c,VEC3d>::value, "parameters must have same type");
//	static_assert(is_dim_of<VEC3a, 3>::value, "The size of the vector must be equal to 3.");

	return in_triangle(P, normal(Ta, Tb, Tc), Ta, Tb,Tc );
}










} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_INCLUSION_H_
