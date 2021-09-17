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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_incident_vertex(const MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_incident_vertex(m, c, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& m, CELL c, const FUNC& func, CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, Graph&> && mesh_traits<MESH>::dimension == 1 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Vertex(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap1&> && mesh_traits<MESH>::dimension == 1 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Vertex(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2 &&
					   (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::Face>))
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Vertex(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge>))
	{
		foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Edge2(c.dart),
							  [&](Dart d) -> bool { return func(Vertex(d)); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Face2(c.dart),
							  [&](Dart d) -> bool { return func(Vertex(d)); });
	}
	else
	{
		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Vertex>(m))
		{
			CellMarkerStore<MESH, Vertex> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Vertex v(d);
				if (!marker.is_marked(v))
				{
					marker.mark(v);
					return func(v);
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
					Vertex v(d);
					foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
						marker.mark(d);
						return true;
					});
					return func(v);
				}
				return true;
			});
		}
	}
}

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename CELL, typename FUNC>
auto foreach_incident_vertex(const IncidenceGraph& ig, CELL c, const FUNC& func)
{
	using Vertex = IncidenceGraph::Vertex;

	static_assert(is_in_tuple<CELL, mesh_traits<IncidenceGraph>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, mesh_traits<IncidenceGraph>::Edge>)
	{
		std::pair<Vertex, Vertex> evs = (*ig.edge_incident_vertices_)[c.index_];
		if (func(evs.first))
			func(evs.second);
	}
	else if constexpr (std::is_same_v<CELL, mesh_traits<IncidenceGraph>::Face>)
	{
		CellMarkerStore<IncidenceGraph, Vertex> marker(ig);
		for (auto& ep : (*ig.face_incident_edges_)[c.index_])
		{
			std::pair<Vertex, Vertex> evs = (*ig.edge_incident_vertices_)[ep.index_];
			bool stop = false;
			if (!marker.is_marked(evs.first))
			{
				marker.mark(evs.first);
				stop = !func(evs.first);
			}
			if (!marker.is_marked(evs.second) && !stop)
			{
				marker.mark(evs.second);
				stop = !func(evs.second);
			}
			if (stop)
				break;
		}
	}
}

/*****************************************************************************/

// template <typename MESH, typename FUNC>
// void foreach_adjacent_vertex_through_edge(MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_adjacent_vertex_through_edge(m, v, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func,
										  CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, Graph&> && mesh_traits<MESH>::dimension == 1)
	{
		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool { return func(Vertex(alpha0(m, d))); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool { return func(Vertex(phi2(m, d))); });
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Vertex>(m))
		{
			CellMarkerStore<MESH, Vertex> marker(m);
			foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
				Vertex av(phi2(m, d));
				if (!marker.is_marked(av))
				{
					marker.mark(av);
					return func(av);
				}
				return true;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
				Vertex av(phi2(m, d));
				if (!marker.is_marked(av.dart))
				{
					foreach_dart_of_orbit(m, av, [&](Dart d) -> bool {
						marker.mark(d);
						return true;
					});
					return func(v);
				}
				return true;
			});
		}
	}
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::vector<typename mesh_traits<MESH>::Vertex> incident_vertices(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::Vertex> incident_vertices(const MESH& m, CELL c)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices;
	vertices.reserve(32u);
	foreach_incident_vertex(m, c, [&](Vertex v) -> bool {
		vertices.push_back(v);
		return true;
	});
	return vertices;
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::vector<typename mesh_traits<MESH>::Vertex> append_incident_vertices(const MESH& m, CELL c, std::vector<typename
// mesh_traits<MESH>::Vertex>& vertices);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
void append_incident_vertices(const MESH& m, CELL c, std::vector<typename mesh_traits<MESH>::Vertex>& vertices)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	foreach_incident_vertex(m, c, [&vertices](Vertex v) -> bool {
		vertices.push_back(v);
		return true;
	});
}

/*****************************************************************************/

// template <typename MESH>
// std::vector<typename mesh_traits<MESH>::Vertex>
// adjacent_vertices_through_edge(MESH& m, typename mesh_traits<MESH>::Vertex v);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH>
std::vector<typename mesh_traits<MESH>::Vertex> adjacent_vertices_through_edge(const MESH& m,
																			   typename mesh_traits<MESH>::Vertex v)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices;
	vertices.reserve(32u);
	foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
		vertices.push_back(av);
		return true;
	});
	return vertices;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_
