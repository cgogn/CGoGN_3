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

#ifndef CGOGN_CORE_TYPES_CMAP_GRAPH_H_
#define CGOGN_CORE_TYPES_CMAP_GRAPH_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap_base.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT Graph : public CMapBase
{
	std::shared_ptr<Attribute<Dart>> alpha0_;
	std::shared_ptr<Attribute<Dart>> alpha1_;
	std::shared_ptr<Attribute<Dart>> alpha_1_;

	using Vertex = Cell<PHI21>;
	using Edge = Cell<PHI2>;

	using Cells = std::tuple<Vertex, Edge>;

	Graph() : CMapBase()
	{
		alpha0_ = add_relation("alpha0");
		alpha1_ = add_relation("alpha1");
		alpha_1_ = add_relation("alpha_1");
	}

	inline Dart alpha0(Dart d) const
	{
		return (*alpha0_)[d.index];
	}

	inline Dart alpha1(Dart d) const
	{
		return (*alpha1_)[d.index];
	}

	inline Dart alpha_1(Dart d) const
	{
		return (*alpha_1_)[d.index];
	}

	/* alpha0 is an involution */
	inline void alpha0_sew(Dart d, Dart e)
	{
		(*alpha0_)[d.index] = e;
		(*alpha0_)[e.index] = d;
	}

	inline void alpha0_unsew(Dart d)
	{
		Dart e = alpha0(d);
		(*alpha0_)[d.index] = d;
		(*alpha0_)[e.index] = e;
	}

	/* alpha1 is a permutation */
	inline void alpha1_sew(Dart d, Dart e)
	{
		Dart f = alpha1(d);
		Dart g = alpha1(e);
		(*alpha1_)[d.index] = g;
		(*alpha1_)[e.index] = f;
		(*alpha_1_)[g.index] = d;
		(*alpha_1_)[f.index] = e;
	}

	inline void alpha1_unsew(Dart d)
	{
		Dart e = alpha1(d);
		Dart f = alpha_1(d);
		(*alpha1_)[f.index] = e;
		(*alpha1_)[d.index] = d;
		(*alpha_1_)[e.index] = f;
		(*alpha_1_)[d.index] = d;
	}

	template <typename CELL, typename FUNC>
	inline void foreach_dart_of_orbit(CELL c, const FUNC& f) const
	{
		static_assert(is_in_tuple<CELL, Cells>::value, "Cell not supported in a Graph");
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		static const Orbit orbit = CELL::ORBIT;
		switch (orbit)
		{
			case DART: f(c.dart); break;
			case PHI2: foreach_dart_of_ALPHA0(c.dart, f); break;
			case PHI21: foreach_dart_of_ALPHA1(c.dart, f); break;
			default: break;
		}
	}

	template <typename FUNC>
	inline void foreach_dart_of_ALPHA0(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		if (f(d))
			f(alpha0(d));
	}

	template <typename FUNC>
	inline void foreach_dart_of_ALPHA1(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		Dart it = d;
		do
		{
			if (!f(it))
				break;
			it = alpha1(it);
		} while (it != d);
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_GRAPH_H_
