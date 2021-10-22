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

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/cmap/cph3.h>
#include <cgogn/core/types/cmap/graph.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// Dart phiX(const MESH& m, Dart d);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

inline Dart phi1(const CMap1& m, Dart d)
{
	return (*(m.phi1_))[d.index];
}

inline Dart phi_1(const CMap1& m, Dart d)
{
	return (*(m.phi_1_))[d.index];
}

inline Dart phi2(const CMap2& m, Dart d)
{
	return (*(m.phi2_))[d.index];
}

inline Dart phi3(const CMap3& m, Dart d)
{
	return (*(m.phi3_))[d.index];
}

inline Dart alpha0(const Graph& m, Dart d)
{
	return (*m.alpha0_)[d.index];
}

inline Dart alpha1(const Graph& m, Dart d)
{
	return (*m.alpha1_)[d.index];
}

inline Dart alpha_1(const Graph& m, Dart d)
{
	return (*m.alpha_1_)[d.index];
}

//////////
// CPH3 //
//////////

Dart phi1(const CPH3& m, Dart d);
Dart phi_1(const CPH3& m, Dart d);
Dart phi2(const CPH3& m, Dart d);
Dart phi3(const CPH3& m, Dart d);

/////////////
// GENERIC //
/////////////

// template <uint64 N, typename MESH>
// inline Dart phi(const MESH& m, Dart d)
// {
// 	unused_parameters(m);
// 	static_assert(N % 10 <= mesh_traits<MESH>::dimension, "Composition of PHI: invalid index (phi1/phi2/phi3 only)");

// 	if constexpr (N % 10 == 1 && mesh_traits<MESH>::dimension >= 1)
// 		return phi1(m, phi<N / 10>(m, d));
// 	else
// 	{
// 		if constexpr (N % 10 == 2 && mesh_traits<MESH>::dimension >= 2)
// 			return phi2(m, phi<N / 10>(m, d));
// 		else
// 		{
// 			if constexpr (N % 10 == 3 && mesh_traits<MESH>::dimension >= 3)
// 				return phi3(m, phi<N / 10>(m, d));
// 			else
// 				return d;
// 		}
// 	}
// }

template <int8... Args, typename MESH>
inline Dart phi(const MESH& m, Dart d)
{
	static_assert(((Args >= -1 && Args <= mesh_traits<MESH>::dimension) && ...), "Bad phi value");

	Dart res = d;
	for (int8 i : {Args...})
	{
		switch (i)
		{
		case -1:
			if constexpr (mesh_traits<MESH>::dimension >= 1)
				res = phi_1(m, res);
			break;
		case 1:
			if constexpr (mesh_traits<MESH>::dimension >= 1)
				res = phi1(m, res);
			break;
		case 2:
			if constexpr (mesh_traits<MESH>::dimension >= 2)
				res = phi2(m, res);
			break;
		case 3:
			if constexpr (mesh_traits<MESH>::dimension >= 3)
				res = phi3(m, res);
			break;
		}
	}
	return res;
}

/*****************************************************************************/

// template <typename MESH>
// Dart phiX_sew(const MESH& m, Dart d, Dart e);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

inline void phi1_sew(CMap1& m, Dart d, Dart e)
{
	Dart f = phi1(m, d);
	Dart g = phi1(m, e);
	(*(m.phi1_))[d.index] = g;
	(*(m.phi1_))[e.index] = f;
	(*(m.phi_1_))[g.index] = d;
	(*(m.phi_1_))[f.index] = e;
}

inline void phi1_unsew(CMap1& m, Dart d)
{
	Dart e = phi1(m, d);
	Dart f = phi1(m, e);
	(*(m.phi1_))[d.index] = f;
	(*(m.phi1_))[e.index] = e;
	(*(m.phi_1_))[f.index] = d;
	(*(m.phi_1_))[e.index] = e;
}

inline void phi2_sew(CMap2& m, Dart d, Dart e)
{
	cgogn_assert(phi2(m, d) == d);
	cgogn_assert(phi2(m, e) == e);
	(*(m.phi2_))[d.index] = e;
	(*(m.phi2_))[e.index] = d;
}

inline void phi2_unsew(CMap2& m, Dart d)
{
	Dart e = phi2(m, d);
	(*(m.phi2_))[d.index] = d;
	(*(m.phi2_))[e.index] = e;
}

inline void phi3_sew(CMap3& m, Dart d, Dart e)
{
	cgogn_assert(phi3(m, d) == d);
	cgogn_assert(phi3(m, e) == e);
	(*(m.phi3_))[d.index] = e;
	(*(m.phi3_))[e.index] = d;
}

inline void phi3_unsew(CMap3& m, Dart d)
{
	Dart e = phi3(m, d);
	(*(m.phi3_))[d.index] = d;
	(*(m.phi3_))[e.index] = e;
}

inline void alpha0_sew(Graph& m, Dart d, Dart e)
{
	(*m.alpha0_)[d.index] = e;
	(*m.alpha0_)[e.index] = d;
}

inline void alpha0_unsew(Graph& m, Dart d)
{
	Dart e = alpha0(m, d);
	(*m.alpha0_)[d.index] = d;
	(*m.alpha0_)[e.index] = e;
}

inline void alpha1_sew(Graph& m, Dart d, Dart e)
{
	Dart f = alpha1(m, d);
	Dart g = alpha1(m, e);
	(*m.alpha1_)[d.index] = g;
	(*m.alpha1_)[e.index] = f;
	(*m.alpha_1_)[g.index] = d;
	(*m.alpha_1_)[f.index] = e;
}

inline void alpha1_unsew(Graph& m, Dart d)
{
	Dart e = alpha1(m, d);
	Dart f = alpha_1(m, d);
	(*m.alpha1_)[f.index] = e;
	(*m.alpha1_)[d.index] = d;
	(*m.alpha_1_)[e.index] = f;
	(*m.alpha_1_)[d.index] = d;
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_PHI_H_
