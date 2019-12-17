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

struct CMapBase;
struct CPH3;

/*****************************************************************************/

// template <typename CMAP>
// Dart phiX(const CMAP& m, Dart d);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CMAP>
auto phi1(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
{
	return (*(m.phi1_))[d.index];
}

template <typename CMAP>
auto phi_1(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
{
	return (*(m.phi_1_))[d.index];
}

template <typename CMAP>
auto phi2(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
{
	return (*(m.phi2_))[d.index];
}

template <typename CMAP>
auto phi3(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
{
	return (*(m.phi3_))[d.index];
}

template <typename CMAP>
auto alpha0(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
{
	return (*m.alpha0_)[d.index];
}

template <typename CMAP>
auto alpha1(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
{
	return (*m.alpha1_)[d.index];
}

template <typename CMAP>
auto alpha_1(const CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>, Dart>
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

template <uint64 N, typename CMAP>
inline Dart phi(const CMAP& m, Dart d)
{
	static_assert(N % 10 <= CMAP::dimension, "Composition of PHI: invalid index (phi1/phi2/phi3 only)");

	if constexpr (N % 10 == 1 && CMAP::dimension >= 1)
		return phi1(m, phi<N / 10>(m, d));
	else
	{
		if constexpr (N % 10 == 2 && CMAP::dimension >= 2)
			return phi2(m, phi<N / 10>(m, d));
		else
		{
			if constexpr (N % 10 == 3 && CMAP::dimension >= 3)
				return phi3(m, phi<N / 10>(m, d));
			else
				return d;
		}
	}
}

/*****************************************************************************/

// template <typename CMAP>
// Dart phiX_sew(const CMAP& m, Dart d, Dart e);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CMAP>
auto phi1_sew(CMAP& m, Dart d, Dart e) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	Dart f = phi1(m, d);
	Dart g = phi1(m, e);
	(*(m.phi1_))[d.index] = g;
	(*(m.phi1_))[e.index] = f;
	(*(m.phi_1_))[g.index] = d;
	(*(m.phi_1_))[f.index] = e;
}

template <typename CMAP>
auto phi1_unsew(CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	Dart e = phi1(m, d);
	Dart f = phi1(m, e);
	(*(m.phi1_))[d.index] = f;
	(*(m.phi1_))[e.index] = e;
	(*(m.phi_1_))[f.index] = d;
	(*(m.phi_1_))[e.index] = e;
}

template <typename CMAP>
auto phi2_sew(CMAP& m, Dart d, Dart e) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	cgogn_assert(phi2(m, d) == d);
	cgogn_assert(phi2(m, e) == e);
	(*(m.phi2_))[d.index] = e;
	(*(m.phi2_))[e.index] = d;
}

template <typename CMAP>
auto phi2_unsew(CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	Dart e = phi2(m, d);
	(*(m.phi2_))[d.index] = d;
	(*(m.phi2_))[e.index] = e;
}

template <typename CMAP>
auto phi3_sew(CMAP& m, Dart d, Dart e) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	cgogn_assert(phi3(m, d) == d);
	cgogn_assert(phi3(m, e) == e);
	(*(m.phi3_))[d.index] = e;
	(*(m.phi3_))[e.index] = d;
}

template <typename CMAP>
auto phi3_unsew(CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	Dart e = phi3(m, d);
	(*(m.phi3_))[d.index] = d;
	(*(m.phi3_))[e.index] = e;
}

template <typename CMAP>
auto alpha0_sew(CMAP& m, Dart d, Dart e) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	(*m.alpha0_)[d.index] = e;
	(*m.alpha0_)[e.index] = d;
}

template <typename CMAP>
auto alpha0_unsew(CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	Dart e = alpha0(m, d);
	(*m.alpha0_)[d.index] = d;
	(*m.alpha0_)[e.index] = e;
}

template <typename CMAP>
auto alpha1_sew(CMAP& m, Dart d, Dart e) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
{
	Dart f = alpha1(m, d);
	Dart g = alpha1(m, e);
	(*m.alpha1_)[d.index] = g;
	(*m.alpha1_)[e.index] = f;
	(*m.alpha_1_)[g.index] = d;
	(*m.alpha_1_)[f.index] = e;
}

template <typename CMAP>
auto alpha1_unsew(CMAP& m, Dart d) -> typename std::enable_if_t<std::is_base_of_v<CMapBase, CMAP>>
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
