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

#ifndef CGOGN_CORE_UTILS_NUMERICS_H_
#define CGOGN_CORE_UTILS_NUMERICS_H_

#include <cstdint>
#include <type_traits>
#include <cmath>
#include <limits>
#include <algorithm>

#include <cgogn/core/utils/assert.h>

namespace cgogn
{

namespace numerics
{

using int8 = std::int8_t;
using int16 = std::int16_t;
using int32 = std::int32_t;
using int64 = std::int64_t;

using uint8 = std::uint8_t;
using uint16 = std::uint16_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

using float32 = float;
using float64 = double;

// FYI MSVC doesn't support std::numeric_limits<uint32>::max() when declaring static const variables
static const uint32 INVALID_INDEX = UINT32_MAX;

template <class Scalar>
inline auto almost_equal_relative(Scalar x, Scalar y, const Scalar max_rel_diff = std::numeric_limits<Scalar>::epsilon() ) -> typename std::enable_if<std::is_floating_point<Scalar>::value, bool>::type
{
	const Scalar diff = std::fabs(x - y);
	x = std::fabs(x);
	y = std::fabs(y);

	return diff <= std::max(x, y) * max_rel_diff;
}

template <class Scalar>
inline auto almost_equal_absolute(Scalar x, Scalar y, const Scalar epsilon = std::numeric_limits<Scalar>::epsilon() ) -> typename std::enable_if<std::is_floating_point<Scalar>::value, bool>::type
{
	cgogn_assert(epsilon > 0);
	return std::fabs(y - x) < epsilon;
}

template <class Scalar, class Integer>
inline Scalar scale_expand_within_0_1(Scalar x, const Integer n)
{
	static_assert(std::is_floating_point<Scalar>::value, "Floating point number required.");
	static_assert(std::is_integral<Integer>::value, "Integer number required.");

	for (Integer i = 1; i <= n; i++)
		x = Real((Scalar(1) - std::cos(Scalar(M_PI) * x)) / Scalar(2));
	for (Integer i = -1; i >= n; i--)
		x = Real(std::acos(Scalar(1) - Scalar(2) * x) / M_PI);
	return x;
}

template  <class Scalar, class Integer>
inline Scalar scale_expand_towards_1(Scalar x, const Integer n)
{
	static_assert(std::is_floating_point<Scalar>::value, "Floating point number required.");
	static_assert(std::is_integral<Integer>::value, "Integer number required.");

	for (Integer i = 1; i <= n; i++)
		x = Real(std::sin(x * Scalar(M_PI_2)));
	for (Integer i = -1; i >= n; i--)
		x = Real(std::asin(x) * Scalar(M_2_PI));
	return x;
}

template <class Scalar>
inline Scalar scale_to_0_1(const Scalar x, const Scalar min, const Scalar max)
{
	static_assert(std::is_floating_point<Scalar>::value, "Floating point number required.");

	return (x - min) / (max - min);
}

template <class Scalar>
inline Scalar scale_and_clamp_to_0_1(const Scalar x, const Scalar min, const Scalar max)
{
	static_assert(std::is_floating_point<Scalar>::value, "Floating point number required.");

	const Scalar v = (x - min) / (max - min);
	return v < Scalar(0) ? Scalar(0) : (v > Scalar(1) ? Scalar(1) : v);
}

template <class Scalar>
inline void scale_centering_around_0(Scalar& min, Scalar& max)
{
	static_assert(std::is_floating_point<Scalar>::value, "Floating point number required.");

	min = std::min(min, -max);
	max = std::max(max, -min);
}

template <class Scalar>
inline Scalar scale_to_0_1_around_one_half(const Scalar x, const Scalar min, const Scalar max)
{
	static_assert(std::is_floating_point<Scalar>::value, "Floating point number required.");

	const Scalar ma = std::max(max, -min);
	const Scalar mi = std::min(min, -max);
	return (x - mi) / (ma - mi);
}

template <typename Scalar>
inline Scalar clamp(const Scalar x, const Scalar min, const Scalar max)
{
	static_assert(std::is_floating_point<Scalar>::value || std::is_integral<Scalar>::value,
		"'clamp' only accept floating-point or integer inputs");

	return std::min(max, std::max(min, x));
}

template<typename T, std::size_t bytes, typename enable = void>
struct fixed_precision {};

template<typename T>
struct fixed_precision<T,1ul,typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value>::type>
{
	using type = int8;
};

template<typename T>
struct fixed_precision<T,2ul,typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value>::type>
{
	using type = int16;
};

template<typename T>
struct fixed_precision<T,4ul,typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value>::type>
{
	using type = int32;
};

template<typename T>
struct fixed_precision<T,8ul,typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value>::type>
{
	using type = int64;
};

template<typename T>
struct fixed_precision<T,1ul,typename std::enable_if<std::is_unsigned<T>::value && std::is_integral<T>::value>::type>
{
	using type = uint8;
};

template<typename T>
struct fixed_precision<T,2ul,typename std::enable_if<std::is_unsigned<T>::value && std::is_integral<T>::value>::type>
{
	using type = uint16;
};

template<typename T>
struct fixed_precision<T,4ul,typename std::enable_if<std::is_unsigned<T>::value && std::is_integral<T>::value>::type>
{
	using type = uint32;
};

template<typename T>
struct fixed_precision<T,8ul,typename std::enable_if<std::is_unsigned<T>::value && std::is_integral<T>::value>::type>
{
	using type = uint64;
};

template<typename T>
struct fixed_precision<T,1ul,typename std::enable_if<std::is_floating_point<T>::value>::type>
{
	using type = float32; // float is at least 32 bits, but we need this specialization for the compilation
};

template<typename T>
struct fixed_precision<T,2ul,typename std::enable_if<std::is_floating_point<T>::value>::type>
{
	using type = float32; // float is at least 32 bits, but we need this specialization for the compilation
};

template<typename T>
struct fixed_precision<T,4ul,typename std::enable_if<std::is_floating_point<T>::value>::type>
{
	using type = float32;
};

template<typename T>
struct fixed_precision<T,8ul,typename std::enable_if<std::is_floating_point<T>::value>::type>
{
	using type = float64;
};

template<typename T, std::size_t sz>
struct fixed_precision<T, sz,typename std::enable_if<!std::is_arithmetic<T>::value>::type>
{
	using type = T;
};

} // namespace numerics

using namespace numerics;

} // namespace cgogn

#endif // CGOGN_CORE_UTILS_NUMERICS_H_
