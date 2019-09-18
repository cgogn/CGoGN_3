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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP2_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP2_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap1.h>
#include <cgogn/core/types/cmap/dart_marker.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap2 : public CMap1
{
	std::shared_ptr<Attribute<Dart>> phi2_;

	using Vertex = Cell<PHI21>;
	using Edge = Cell<PHI2>;
	using Face = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Volume;

	using Cells = std::tuple<Vertex, Edge, Face, Volume>;

	CMap2() : CMap1()
	{
		phi2_ = add_relation("phi2");
	}

	inline Dart phi2(Dart d) const
	{
		return (*phi2_)[d.index];
	}

	template <uint64 N>
	inline Dart phi(Dart d) const
	{
		static_assert(N % 10 <= 2, "Composition of PHI: invalid index (phi1/phi2 only)");
		switch (N % 10)
		{
			case 1: return phi1(phi<N / 10>(d));
			case 2: return phi2(phi<N / 10>(d));
			default: return d;
		}
	}

	inline void phi2_sew(Dart d, Dart e)
	{
		cgogn_assert(phi2(d) == d);
		cgogn_assert(phi2(e) == e);
		(*phi2_)[d.index] = e;
		(*phi2_)[e.index] = d;
	}

	inline void phi2_unsew(Dart d)
	{
		Dart e = phi2(d);
		(*phi2_)[d.index] = d;
		(*phi2_)[e.index] = e;
	}

	template <typename CELL, typename FUNC>
	inline void foreach_dart_of_orbit(CELL c, const FUNC& f) const
	{
		static_assert(is_in_tuple<CELL, Cells>::value, "Cell not supported in a CMap2");
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		static const Orbit orbit = CELL::ORBIT;
		switch (orbit)
		{
			case DART: f(c.dart); break;
			case PHI1: foreach_dart_of_PHI1(c.dart, f); break;
			case PHI2: foreach_dart_of_PHI2(c.dart, f); break;
			case PHI21: foreach_dart_of_PHI21(c.dart, f); break;
			case PHI1_PHI2: foreach_dart_of_PHI1_PHI2(c.dart, f); break;
			default: break;
		}
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI2(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		if (f(d))
			f(phi2(d));
	}

	template <typename FUNC>
	inline void foreach_dart_of_PHI21(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		Dart it = d;
		do
		{
			if (!f(it))
				break;
			it = phi2(phi_1(it));
		} while (it != d);
	}

	template <typename FUNC>
	void foreach_dart_of_PHI1_PHI2(Dart d, const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		DartMarkerStore marker(*this);

		std::vector<Dart> visited_faces;
		visited_faces.push_back(d); // Start with the face of d

		// For every face added to the list
		for (uint32 i = 0; i < visited_faces.size(); ++i)
		{
			const Dart e = visited_faces[i];
			if (!marker.is_marked(e))	// Face has not been visited yet
			{
				// mark visited darts (current face)
				// and add non visited adjacent faces to the list of face
				Dart it = e;
				do
				{
					if (!f(it)) // apply the function to the darts of the face
						return;
					marker.mark(it);				// Mark
					const Dart adj = phi2(it);		// Get adjacent face
					if (!marker.is_marked(adj))
						visited_faces.push_back(adj);	// Add it
					it = phi1(it);
				} while (it != e);
			}
		}
	}

	Dart close_hole(Dart d, bool set_indices = true);

	uint32 close(bool set_indices = true);
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP2_H_
