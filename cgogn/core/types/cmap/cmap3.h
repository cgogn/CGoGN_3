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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP3_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP3_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap2.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap3 : public CMap2
{
	std::shared_ptr<Attribute<Dart>> phi3_;

	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using Edge = Cell<PHI2_PHI3>;
	using Face = Cell<PHI1_PHI3>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Cell<PHI1_PHI2_PHI3>;

	using Cells = std::tuple<Vertex, Vertex2, Edge, Face, Volume>;

	CMap3() : CMap2()
	{
		phi3_ = add_relation("phi3");
	}

	inline Dart phi3(Dart d) const
	{
		return (*phi3_)[d.index];
	}

	template <uint64 N>
	inline Dart phi(Dart d) const
	{
		static_assert(N % 10 <= 3, "Composition of PHI: invalid index (phi1/phi2/phi3 only)");
		switch (N % 10)
		{
			case 1: return phi1(phi<N / 10>(d));
			case 2: return phi2(phi<N / 10>(d));
			case 3: return phi3(phi<N / 10>(d));
			default: return d;
		}
	}

	inline void phi3_sew(Dart d, Dart e)
	{
		cgogn_assert(phi3(d) == d);
		cgogn_assert(phi3(e) == e);
		(*phi3_)[d.index] = e;
		(*phi3_)[e.index] = d;
	}

	inline void phi3_unsew(Dart d)
	{
		Dart e = phi3(d);
		(*phi3_)[d.index] = d;
		(*phi3_)[e.index] = e;
	}

	template <typename CELL, typename FUNC>
	inline void foreach_dart_of_orbit(CELL c, const FUNC& f) const
	{
		static_assert(is_in_tuple<CELL, Cells>::value, "Cell not supported in a CMap3");
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		static const Orbit orbit = CELL::ORBIT;
		switch (orbit)
		{
			case DART: f(c.dart); break;
			case PHI21: foreach_dart_of_PHI21(c.dart, f); break;
			case PHI1_PHI2: foreach_dart_of_PHI1_PHI2(c.dart, f); break;
			case PHI1_PHI3: foreach_dart_of_PHI1_PHI3(c.dart, f); break;
			case PHI2_PHI3: foreach_dart_of_PHI2_PHI3(c.dart, f); break;
			case PHI21_PHI31: foreach_dart_of_PHI21_PHI31(c.dart, f); break;
			case PHI1_PHI2_PHI3: foreach_dart_of_PHI1_PHI2_PHI3(c.dart, f); break;
			default: break;
		}
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI1_PHI3(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		foreach_dart_of_PHI1(d, [&] (Dart fd) -> bool
		{
			if (f(fd)) return f(phi3(fd));
			return false;
		});
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI2_PHI3(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		Dart it = d;
		do
		{
			if (!f(it)) break;
			it = phi2(it);
			if (!f(it)) break;
			it = phi3(it);
		} while (it != d);
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI21_PHI31(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		DartMarkerStore marker(*this);
		const std::vector<Dart>& marked_darts = marker.marked_darts();

		marker.mark(d);
		for (uint32 i = 0; i < marked_darts.size(); ++i)
		{
			const Dart curr_dart = marked_darts[i];
//			if ( !(is_boundary(curr_dart) && is_boundary(phi3(curr_dart))) )
				if (!f(curr_dart))
					break;

			const Dart d_1 = phi_1(curr_dart);
			const Dart d2_1 = phi2(d_1); // turn in volume
			const Dart d3_1 = phi3(d_1); // change volume

			if (!marker.is_marked(d2_1))
				marker.mark(d2_1);
			if (!marker.is_marked(d3_1))
				marker.mark(d3_1);
		}
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI1_PHI2_PHI3(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		DartMarkerStore marker(*this);

		std::vector<Dart> visited_face2;
		visited_face2.push_back(d); // Start with the face of d

		// For every face added to the list
		for (uint32 i = 0; i < visited_face2.size(); ++i)
		{
			const Dart e = visited_face2[i];
			if (!marker.is_marked(e))	// Face2 has not been visited yet
			{
				// mark visited darts (current face2)
				// and add non visited phi2-adjacent face2 to the list of face2
				Dart it = e;
				do
				{
					if (!f(it)) // apply the function to the darts of the face2
						return;
					marker.mark(it);				// Mark
					const Dart adj2 = this->phi2(it);	// Get phi2-adjacent face2
					if (!marker.is_marked(adj2))
						visited_face2.push_back(adj2);	// Add it
					it = phi1(it);
				} while (it != e);
				// add phi3-adjacent face2 to the list
				visited_face2.push_back(phi3(it));
			}
		}
	}

	Volume close_hole(Dart d, bool set_indices = true);

	inline void sew_volumes(Dart d0, Dart d1)
	{
		// cgogn_message_assert(codegree(*this, CMap3::Face(d0)) == codegree(*this, CMap3::Face(d1)), "the two faces to sow do not have the same codegree");

		Dart it0 = d0;
		Dart it1 = d1;
		do
		{
			phi3_sew(it0, it1);
			it0 = phi1(it0);
			it1 = phi_1(it1);
		} while (it0 != d0);
	}

	uint32 close(bool set_indices = true);
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP3_H_
