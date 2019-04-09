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

#include <cgogn/geometry/types/geometry_traits.h>

namespace cgogn
{

namespace geometry
{

/**
 * normal of the plane spanned by 3 points in 3D
 */
template <typename VEC3a, typename VEC3b, typename VEC3c>
inline typename vector_traits<VEC3a>::Type normal(const Eigen::MatrixBase<VEC3a>& p1, const Eigen::MatrixBase<VEC3b>& p2, const Eigen::MatrixBase<VEC3c>& p3)
{
	static_assert(is_same_vector<VEC3a,VEC3b,VEC3c>::value, "parameters must have same type");
	static_assert(is_dim_of<VEC3a, 3>::value, "The size of the vector must be equal to 3.");
	return (p2-p1).cross(p3-p1);
}




template <typename VEC3>
inline auto normal(const VEC3& p1, const VEC3& p2, const VEC3& p3)
 -> typename std::enable_if <is_vec_non_eigen<VEC3>::value, VEC3>::type
{
	static_assert(is_dim_of<VEC3, 3>::value, "The size of the vector must be equal to 3.");
	return copy_to_vec<VEC3>(normal(eigenize(p1),eigenize(p2),eigenize(p3)));
}


} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_NORMAL_H_
