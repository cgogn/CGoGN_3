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

#ifndef CGOGN_CORE_CMAP_ORBIT_TRAVERSAL_H_
#define CGOGN_CORE_CMAP_ORBIT_TRAVERSAL_H_

#include <cgogn/core/types/cmap/cell.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/phi.h>

#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename CELL, typename FUNC>
auto foreach_dart_of_orbit(const MESH& m, CELL c, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	static const Orbit orbit = CELL::ORBIT;

	if constexpr (orbit == DART)
	{
		unused_parameters(m);
		f(c.dart);
		return;
	}
	if constexpr (orbit == PHI1)
	{
		foreach_dart_of_PHI1(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI2)
	{
		foreach_dart_of_PHI2(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI21)
	{
		foreach_dart_of_PHI21(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI1_PHI2)
	{
		foreach_dart_of_PHI1_PHI2(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI1_PHI3)
	{
		foreach_dart_of_PHI1_PHI3(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI2_PHI3)
	{
		foreach_dart_of_PHI2_PHI3(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI21_PHI31)
	{
		foreach_dart_of_PHI21_PHI31(m, c.dart, f);
		return;
	}
	if constexpr (orbit == PHI1_PHI2_PHI3)
	{
		foreach_dart_of_PHI1_PHI2_PHI3(m, c.dart, f);
		return;
	}
}

///////////
// Graph //
///////////

template <typename CELL, typename FUNC>
void foreach_dart_of_orbit(const Graph& m, CELL c, const FUNC& f)
{
	static_assert(is_in_tuple<CELL, typename Graph::Cells>::value, "Cell not supported in a Graph");
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	static const Orbit orbit = CELL::ORBIT;
	switch (orbit)
	{
	case DART:
		f(c.dart);
		break;
	case PHI2:
		foreach_dart_of_ALPHA0(m, c.dart, f);
		break;
	case PHI21:
		foreach_dart_of_ALPHA1(m, c.dart, f);
		break;
	default:
		break;
	}
}

/*****************************************************************************/
// orbits traversals
/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI1(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = d;
	do
	{
		if (!f(it))
			break;
		it = phi1(m, it);
	} while (it != d);
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (f(d))
		f(phi2(m, d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI21(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = d;
	do
	{
		if (!f(it))
			break;
		it = phi<-1, 2>(m, it);
	} while (it != d);
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI1_PHI2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	DartMarkerStore<MESH> marker(m);

	std::vector<Dart> visited_faces;
	visited_faces.push_back(d); // Start with the face of d

	// For every face added to the list
	for (uint32 i = 0; i < uint32(visited_faces.size()); ++i)
	{
		const Dart e = visited_faces[i];
		if (!marker.is_marked(e)) // Face has not been visited yet
		{
			// mark visited darts (current face)
			// and add non visited adjacent faces to the list of face
			Dart it = e;
			do
			{
				if (!f(it)) // apply the function to the darts of the face
					return;
				marker.mark(it);			  // Mark
				const Dart adj = phi2(m, it); // Get adjacent face
				if (!marker.is_marked(adj))
					visited_faces.push_back(adj); // Add it
				it = phi1(m, it);
			} while (it != e);
		}
	}
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI1_PHI3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_PHI1(m, d, [&](Dart fd) -> bool {
		if (f(fd))
			return f(phi3(m, fd));
		return false;
	});
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI2_PHI3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = d;
	do
	{
		if (!f(it))
			break;
		it = phi2(m, it);
		if (!f(it))
			break;
		it = phi3(m, it);
	} while (it != d);
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI21_PHI31(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		//			if ( !(is_boundary(curr_dart) && is_boundary(phi3(curr_dart))) )
		if (!f(curr_dart))
			break;

		const Dart d_1 = phi_1(m, curr_dart);
		const Dart d2_1 = phi2(m, d_1); // turn in volume
		const Dart d3_1 = phi3(m, d_1); // change volume

		if (!marker.is_marked(d2_1))
			marker.mark(d2_1);
		if (!marker.is_marked(d3_1))
			marker.mark(d3_1);
	}
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_PHI1_PHI2_PHI3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	DartMarkerStore<MESH> marker(m);

	std::vector<Dart> visited_face2;
	visited_face2.push_back(d); // Start with the face of d

	// For every face added to the list
	for (uint32 i = 0; i < visited_face2.size(); ++i)
	{
		const Dart e = visited_face2[i];
		if (!marker.is_marked(e)) // Face2 has not been visited yet
		{
			// mark visited darts (current face2)
			// and add non visited phi2-adjacent face2 to the list of face2
			Dart it = e;
			do
			{
				if (!f(it)) // apply the function to the darts of the face2
					return;
				marker.mark(it);			   // Mark
				const Dart adj2 = phi2(m, it); // Get phi2-adjacent face2
				if (!marker.is_marked(adj2))
					visited_face2.push_back(adj2); // Add it
				it = phi1(m, it);
			} while (it != e);
			// add phi3-adjacent face2 to the list
			visited_face2.push_back(phi3(m, it));
		}
	}
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_ALPHA0(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (f(d))
		f(alpha0(m, d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_ALPHA1(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = d;
	do
	{
		if (!f(it))
			break;
		it = alpha1(m, it);
	} while (it != d);
}

} // namespace cgogn

#endif // CGOGN_CORE_CMAP_ORBIT_TRAVERSAL_H_
