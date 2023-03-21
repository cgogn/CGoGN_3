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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_edge(MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_edge(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_incident_edge(m, c, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_edge(const MESH& m, CELL c, const FUNC& func, CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Edge = typename mesh_traits<MESH>::Edge;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, Graph&> && mesh_traits<MESH>::dimension == 1 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex>)
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap1&> && mesh_traits<MESH>::dimension == 1 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2 &&
					   (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::Face>))
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Face2(c.dart),
							  [&](Dart d) -> bool { return func(Edge(d)); });
	}
	else
	{
		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Edge>(m))
		{
			CellMarkerStore<MESH, Edge> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Edge e(d);
				if (!marker.is_marked(e))
				{
					marker.mark(e);
					return func(e);
				}
				return true;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				if (!marker.is_marked(d))
				{
					Edge e(d);
					foreach_dart_of_orbit(m, e, [&](Dart d) -> bool {
						marker.mark(d);
						return true;
					});
					return func(e);
				}
				return true;
			});
		}
	}
}

//////////////////////
/// IncidenceGraph ///
//////////////////////

template <typename CELL, typename FUNC>
auto foreach_incident_edge(const IncidenceGraph& ig, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this IncidenceGraph");
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, IncidenceGraph::Vertex>)
	{
		for (auto& ep : (*ig.vertex_incident_edges_)[c.index_])
		{
			if (!func(ep))
				break;
		}
	}
	else if constexpr (std::is_same_v<CELL, IncidenceGraph::Face>)
	{
		for (auto& ep : (*ig.face_incident_edges_)[c.index_])
		{
			if (!func(ep))
				break;
		}
	}
}

/*****************************************************************************/

// template <typename MESH, typename FUNC>
// void foreach_adjacent_edge_through_face(MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_adjacent_edge_through_face(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_adjacent_edge_through_face(m, e, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_edge_through_face(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func,
										CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Edge = typename mesh_traits<MESH>::Edge;

	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		using Face = typename mesh_traits<MESH>::Face;

		Dart d1 = e.dart;
		Dart d2 = phi2(m, d1);
		if (!is_boundary(m, d1))
			foreach_dart_of_orbit(m, Face(d1), [&](Dart d) -> bool { return d != d1 ? func(Edge(d)) : true; });
		if (!is_boundary(m, d2))
			foreach_dart_of_orbit(m, Face(d2), [&](Dart d) -> bool { return d != d2 ? func(Edge(d)) : true; });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		using Face2 = typename mesh_traits<MESH>::Face2;

		foreach_dart_of_orbit(m, e, [&](Dart ed) -> bool {
			foreach_dart_of_orbit(m, Face2(ed), [&](Dart d) -> bool { return d != ed ? func(Edge(d)) : true; });
			return true;
		});
	}
}

//////////////////////
/// IncidenceGraph ///
//////////////////////

template <typename FUNC>
void foreach_adjacent_edge_through_face(const IncidenceGraph& ig, IncidenceGraph::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	bool stop = false;
	CellMarkerStore<IncidenceGraph, IncidenceGraph::Edge> marker(ig);
	marker.mark(e);
	foreach_incident_face(ig, e, [&](IncidenceGraph::Face f0) -> bool {
		foreach_incident_edge(ig, f0, [&](IncidenceGraph::Edge e1) -> bool {
			if (!marker.is_marked(e1))
			{
				marker.mark(e1);
				stop = !func(e1);
			}
			return !stop;
		});
		return !stop;
	});
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Edge> incident_edges(MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::Edge> incident_edges(const MESH& m, CELL c)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	std::vector<Edge> edges;
	edges.reserve(32u);
	foreach_incident_edge(m, c, [&](Edge e) -> bool {
		edges.push_back(e);
		return true;
	});
	return edges;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_
