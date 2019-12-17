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

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_incident_vertex(const MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

// this version works for any CMap and CELL

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMapBase>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_indexed<Vertex>(m))
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
			Vertex v(d);
			if (!marker.is_marked(d))
			{
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

// below are specialized versions for cases that do not need marking (ordered neighborhoods)

////////////////////////////
// Graph (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_vertex(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, Graph>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, e, [&](Dart d) -> bool { return func(Vertex(d)); });
}

////////////////////////////
// CMap1 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_vertex(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap1>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, f, [&](Dart d) -> bool { return func(Vertex(d)); });
}

////////////////////////////
// CMap2 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_vertex(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap2>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, e, [&](Dart d) -> bool { return func(Vertex(d)); });
}

template <typename MESH, typename FUNC>
auto foreach_incident_vertex(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap2>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, f, [&](Dart d) -> bool { return func(Vertex(d)); });
}

////////////////////////////
// CMap3 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_vertex(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap3>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Edge2(e.dart),
						  [&](Dart d) -> bool { return func(Vertex(d)); });
}

template <typename MESH, typename FUNC>
auto foreach_incident_vertex(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap3>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Face2(f.dart),
						  [&](Dart d) -> bool { return func(Vertex(d)); });
}

/*****************************************************************************/

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_adjacent_vertex_through_edge(MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& f);

/*****************************************************************************/

////////////////////////////
// Graph (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, Graph>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(v, [&](Dart d) -> bool { return func(Vertex(alpha0(m, d))); });
}

////////////////////////////
// CMap2 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH, CMap2>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m, v, [&](Dart d) -> bool { return func(Vertex(phi2(m, d))); });
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_
