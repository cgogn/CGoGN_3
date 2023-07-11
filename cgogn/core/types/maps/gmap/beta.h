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

#ifndef CGOGN_CORE_TYPES_GMAP_BETA_H_
#define CGOGN_CORE_TYPES_GMAP_BETA_H_

#include <cgogn/core/types/maps/dart.h>
#include <cgogn/core/types/mesh_traits.h>


namespace cgogn
{

template <int8 Arg, int8... Args, typename MESH>
auto beta(const MESH& m, Dart d) -> std::enable_if_t < std::is_convertible_v<MESH&, GMapBase&>, Dart >
{
	static_assert((Arg >= 0 && Arg <= mesh_traits<MESH>::dimension), "Bad phi value");

	Dart res;
	if constexpr (Arg == 0)
		res = beta0(m, d);
	if constexpr (Arg == 1)
		res = beta1(m, d);
	if constexpr (Arg == 2)
		res = beta2(m, d);
	if constexpr (Arg == 3)
		res = beta3(m, d);

	if constexpr (sizeof...(Args) > 0)
		return beta<Args...>(m, res);
	else
		return res;
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_GMAP_BETA_H_
