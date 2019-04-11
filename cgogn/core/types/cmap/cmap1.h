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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP1_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP1_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap0.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap1 : public CMap0
{
	std::shared_ptr<Attribute<Dart>> phi1_;
	std::shared_ptr<Attribute<Dart>> phi_1_;

	using Vertex = Cell<DART>;
	using Edge = Cell<DART>;
	using Face = Cell<PHI1>;
	using CC = Face;

	using Cells = std::tuple<Vertex, Edge, Face>;

	CMap1() : CMap0()
	{
		phi1_ = add_relation("phi1");
		phi_1_ = add_relation("phi_1");
	}

	inline Dart phi1(Dart d) const
	{
		return (*phi1_)[d.index];
	}

	inline Dart phi_1(Dart d) const
	{
		return (*phi_1_)[d.index];
	}

	template <uint64 N>
	inline Dart phi(Dart d) const
	{
		static_assert(N % 10 <= 1, "Composition of PHI: invalid index");
		if (N >= 10)
			return phi1(phi<N / 10>(d));
		if (N == 1)
			return phi1(d);
		return d;
	}

	inline void phi1_sew(Dart d, Dart e)
	{
		Dart f = phi1(d);
		Dart g = phi1(e);
		(*phi1_)[d.index] = g;
		(*phi1_)[e.index] = f;
		(*phi_1_)[g.index] = d;
		(*phi_1_)[f.index] = e;
	}

	inline void phi1_unsew(Dart d)
	{
		Dart e = phi1(d);
		Dart f = phi1(e);
		(*phi1_)[d.index] = f;
		(*phi1_)[e.index] = e;
		(*phi_1_)[f.index] = d;
		(*phi_1_)[e.index] = e;
	}

	template <typename CELL, typename FUNC>
	inline void foreach_dart_of_orbit(CELL c, const FUNC& f) const
	{
		static_assert(is_in_tuple<CELL, Cells>::value, "Cell not supported in a CMap1");
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		static const Orbit orbit = CELL::ORBIT;
		switch (orbit)
		{
			case DART: f(c.dart); break;
			case PHI1: foreach_dart_of_PHI1(c.dart, f); break;
		}
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI1(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		Dart it = d;
		do
		{
			if (!f(it))
				break;
			it = phi1(it);
		} while (it != d);
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP1_H_
