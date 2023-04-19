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

#ifndef CGOGN_CORE_TYPES_CMAP_PHI_H_
#define CGOGN_CORE_TYPES_CMAP_PHI_H_

namespace cgogn
{

template <int8 Arg, int8... Args, typename MESH>
inline Dart phi(const MESH& m, Dart d)
{
	static_assert((Arg >= -1 && Arg <= mesh_traits<MESH>::dimension), "Bad phi value");

	Dart res;
	if constexpr (Arg == -1)
		res = phi_1(m, d);
	if constexpr (Arg == 1)
		res = phi1(m, d);
	if constexpr (Arg == 2)
		res = phi2(m, d);
	if constexpr (Arg == 3)
		res = phi3(m, d);

	if constexpr (sizeof...(Args) > 0)
		return phi<Args...>(m, res);
	else
		return res;

	// static_assert(((Args >= -1 && Args <= mesh_traits<MESH>::dimension) && ...), "Bad phi value");

	// Dart res = d;
	// for (int8 i : {Args...})
	// {
	// 	switch (i)
	// 	{
	// 	case -1:
	// 		if constexpr (mesh_traits<MESH>::dimension >= 1)
	// 			res = phi_1(m, res);
	// 		break;
	// 	case 1:
	// 		if constexpr (mesh_traits<MESH>::dimension >= 1)
	// 			res = phi1(m, res);
	// 		break;
	// 	case 2:
	// 		if constexpr (mesh_traits<MESH>::dimension >= 2)
	// 			res = phi2(m, res);
	// 		break;
	// 	case 3:
	// 		if constexpr (mesh_traits<MESH>::dimension >= 3)
	// 			res = phi3(m, res);
	// 		break;
	// 	}
	// }
	// return res;
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_PHI_H_
