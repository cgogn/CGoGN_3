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

#ifndef CGOGN_CORE_UTILS_TUPLES_H_
#define CGOGN_CORE_UTILS_TUPLES_H_

#include <tuple>

namespace cgogn
{

template <typename V, typename T>
struct is_in_tuple;

template <typename V, typename T0, typename... T>
struct is_in_tuple<V, std::tuple<T0, T...>>
{
	static const bool value = is_in_tuple<V, std::tuple<T...>>::value;
};

template <typename V, typename... T>
struct is_in_tuple<V, std::tuple<V, T...>>
{
	static const bool value = true;
};

template <typename V>
struct is_in_tuple<V, std::tuple<>>
{
	static const bool value = false;
};

template <typename V, typename T>
inline constexpr bool is_in_tuple_v = is_in_tuple<V, T>::value;

template <class T, class Tuple>
struct tuple_type_index;

template <class T, class... Types>
struct tuple_type_index<T, std::tuple<T, Types...>>
{
	static const std::size_t value = 0;
};

template <class T, class U, class... Types>
struct tuple_type_index<T, std::tuple<U, Types...>>
{
	static const std::size_t value = 1 + tuple_type_index<T, std::tuple<Types...>>::value;
};

} // namespace cgogn

#endif // CGOGN_CORE_UTILS_TUPLES_H_
