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

#ifndef CGOGN_CORE_UTILS_TYPE_TRAITS_H_
#define CGOGN_CORE_UTILS_TYPE_TRAITS_H_

#include <functional>
#include <tuple>

namespace cgogn
{

// template <bool>
// struct void_
// {
// 	typedef void type;
// };

// template <typename T, typename = void>
// struct is_mesh_view
// {
// 	static const bool value = false;
// };

// template <typename T>
// struct is_mesh_view<T, typename void_<T::is_mesh_view>::type>
// {
// 	static bool const value = true;
// };

namespace internal
{

namespace type_traits
{

/**
 * function_traits
 * Traits class to inspect function characteristics (return type, arity, parameters types)
 * Warning : when dealing with a member function, the pointer to the current object is ignored.
 */

// specialization for lambda functions
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{
};

// General case
template <typename ReturnType, typename... Args>
struct function_traits<ReturnType(Args...)>
{
	static const size_t arity = sizeof...(Args);

	using result_type = ReturnType;

	template <size_t i>
	struct arg
	{
		static_assert(i < sizeof...(Args),
					  "Trying to access to an argument whose index is higher than the function arity.");
		using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
		// the i-th argument is equivalent to the i-th tuple element of a tuple composed of those arguments.
	};
};

// specialization for function pointers
template <typename ReturnType, typename... Args>
struct function_traits<ReturnType (*)(Args...)> : public function_traits<ReturnType(Args...)>
{
};

// specialization for function references
template <typename ReturnType, typename... Args>
struct function_traits<ReturnType (&)(Args...)> : public function_traits<ReturnType(Args...)>
{
};

// specialization for member function pointers
template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType (ClassType::*)(Args...)> : public function_traits<ReturnType(Args...)>
{
};

// specialization for const member function pointers
template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType (ClassType::*)(Args...) const> : public function_traits<ReturnType(Args...)>
{
};

} // namespace type_traits

} // namespace internal

template <typename F>
using func_arity = std::integral_constant<std::size_t, internal::type_traits::function_traits<F>::arity>;

template <typename F>
using func_parameter_type = typename internal::type_traits::function_traits<F>::template arg<0>::type;

template <typename F, std::size_t i>
using func_ith_parameter_type = typename internal::type_traits::function_traits<F>::template arg<i>::type;

template <typename F, typename T>
using is_func_parameter_same = std::is_same<func_parameter_type<F>, T>;

template <typename F, std::size_t i, typename T>
using is_ith_func_parameter_same = std::is_same<func_ith_parameter_type<F, i>, T>;

template <typename F>
using func_return_type = typename internal::type_traits::function_traits<F>::result_type;

template <typename F, typename T>
using is_func_return_same = std::is_same<func_return_type<F>, T>;

} // namespace cgogn

#endif // CGOGN_CORE_UTILS_TYPE_TRAITS_H_
