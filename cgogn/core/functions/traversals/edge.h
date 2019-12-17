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
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

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

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_edge(MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

// this version works for any CMap and CELL

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_edge(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMapBase>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_indexed<Edge>(m))
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
			Edge e(d);
			if (!marker.is_marked(d))
			{
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

// below are specialized versions for cases that do not need marking (ordered neighborhoods)

////////////////////////////
// Graph (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, Graph>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, v, [&](Dart d) -> bool { return func(Edge(d)); });
}

////////////////////////////
// CMap1 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap1>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, f, [&](Dart d) -> bool { return func(Edge(d)); });
}

////////////////////////////
// CMap2 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap2>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, v, [&](Dart d) -> bool { return func(Edge(d)); });
}

template <typename MESH, typename FUNC>
auto foreach_incident_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap2>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, f, [&](Dart d) -> bool { return func(Edge(d)); });
}

////////////////////////////
// CMap3 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap3>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Face2(f.dart), [&](Dart d) -> bool { return func(Edge(d)); });
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_
