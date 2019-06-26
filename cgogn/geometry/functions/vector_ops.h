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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_VECTOR_OPS_H_
#define CGOGN_GEOMETRY_FUNCTIONS_VECTOR_OPS_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

// template <typename VEC,
// 		  typename = typename std::enable_if<is_eigen<VEC>::value>::type>
// void
// normalize(VEC& v)
// {
// 	using Scalar = typename vector_traits<VEC>::Scalar;
// 	const Scalar norm2 = v.squaredNorm();
// 	if (norm2 > Scalar(0))
// 		v /= std::sqrt(norm2);
// }

// template <typename VEC,
// 		  typename = typename std::enable_if<is_eigen<VEC>::value>::type>
// typename vector_traits<VEC>::Scalar
// norm(VEC& v)
// {
// 	return v.norm();
// }

// template <typename VEC,
// 		  typename = typename std::enable_if<is_eigen<VEC>::value>::type>
// typename vector_traits<VEC>::Scalar
// squared_norm(VEC& v)
// {
// 	v.squaredNorm();
// }

template <typename VEC,
		  typename std::enable_if<(vector_traits<VEC>::SIZE > 1) && is_eigen<VEC>::value>::type* = nullptr>
void
set_zero(VEC& v)
{
	v.setZero();
}

template <typename VEC,
		  typename std::enable_if<(vector_traits<VEC>::SIZE == 1)>::type* = nullptr>
void
set_zero(VEC& v)
{
	v = 0;
}

// template <typename VEC3,
// 		  typename = typename std::enable_if<is_eigen<VEC3>::value>::type>
// VEC3
// cross(const VEC3& v1, const VEC3& v2)
// {
// 	static_assert (vector_traits<VEC3>::SIZE == 3, "vec_ops: cross product is only defined for vectors of dimension 3");
// 	return v1.cross(v2);
// }

// void transpose(const VEC3& v)
// {

// }

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_VECTOR_OPS_H_
