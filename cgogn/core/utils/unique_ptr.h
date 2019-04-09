/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *                                                                  *                                                                              *
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

#ifndef CGOGN_CORE_UTILS_UNIQUE_PTR_H_
#define CGOGN_CORE_UTILS_UNIQUE_PTR_H_

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

/*
 *  std::make_unique is c++14 but std::make_shared is c++11
 *  below you can find the iso implementation of make_unique
 *  see https://isocpp.org/files/papers/N3656.txt
 */

namespace cgogn
{

template <class T>
struct _Unique_if
{
	using _Single_object = std::unique_ptr<T>;
};

template <class T>
struct _Unique_if<T[]>
{
	using _Unknown_bound = std::unique_ptr<T[]>;
};

template <class T, size_t N>
struct _Unique_if<T[N]>
{
	using _Known_bound = void;
};

template <class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args&&... args)
{
	return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n)
{
	using U = typename std::remove_extent<T>::type;
	return std::unique_ptr<T>(new U[n]());
}

template <class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;

/**
 * @brief dynamic_cast_unique_ptr
 * @param ptr, a unique_ptr that might loose the ownership.
 * @return nullptr if the cast wasn't successfull, a valid unique_ptr<TO> otherwise
 * Warning : if the cast is successfull, the unique_ptr "ptr" loses the ownership (ptr.release() is called), whereas if the cast fails ptr still has the ownership.
 */
template <typename TO, typename FROM>
inline std::unique_ptr<TO> dynamic_cast_unique_ptr(std::unique_ptr<FROM>&& ptr)
{
	TO* const res = dynamic_cast<TO*>(ptr.get());
	if (res != nullptr)
		ptr.release();
	return std::unique_ptr<TO>(res);
}

} // namespace cgogn

#endif // CGOGN_CORE_UTILS_UNIQUE_PTR_H_
