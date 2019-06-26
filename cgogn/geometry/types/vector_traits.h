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

#ifndef CGOGN_GEOMETRY_TYPES_VECTOR_TRAITS_H_
#define CGOGN_GEOMETRY_TYPES_VECTOR_TRAITS_H_

#include <Eigen/Dense>

namespace cgogn
{

namespace geometry
{

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Vec4 = Eigen::Vector4d;

using Vec2f = Eigen::Vector2f;
using Vec3f = Eigen::Vector3f;
using Vec4f = Eigen::Vector4f;

using Vec2i = Eigen::Vector2i;
using Vec3i = Eigen::Vector3i;
using Vec4i = Eigen::Vector4i;

using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix3d;


template <typename VEC, typename Enable = void>
struct vector_traits;

// specialization for uniform manip of vec & scalar (vbo)

template <typename T>
struct vector_traits<T, typename std::enable_if<std::is_integral<T>::value || std::is_floating_point<T>::value>::type>
{
	static const std::size_t SIZE = 1;
	using Scalar = T;
};


// specialization : Eigen::Vector

template <typename T>
std::true_type cgogn_check_eigen_type(const Eigen::MatrixBase<T>*);
std::false_type cgogn_check_eigen_type(...);
template <typename T>
struct is_eigen : public decltype(cgogn_check_eigen_type(std::declval<T*>()))
{};

template <typename V>
struct vector_traits<V, typename std::enable_if<is_eigen<V>::value>::type>
{
	static const std::size_t SIZE = Eigen::internal::traits<V>::RowsAtCompileTime;
	using Scalar = typename Eigen::internal::traits<V>::Scalar;
};

template <typename V>
struct vector_traits<Eigen::MatrixBase<V>, typename std::enable_if<is_eigen<Eigen::MatrixBase<V>>::value>::type>
{
	static const std::size_t SIZE = Eigen::internal::traits<V>::RowsAtCompileTime;
	using Scalar = typename Eigen::internal::traits<V>::Scalar;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_VECTOR_TRAITS_H_
