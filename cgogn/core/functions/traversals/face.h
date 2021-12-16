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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_

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
// void foreach_incident_face(MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_face(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_incident_face(m, c, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_face(const MESH& m, CELL c, const FUNC& func, CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2 &&
				  (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex> ||
				   std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge> ||
				   std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>))
	{
		foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
			if (!is_boundary(m, d))
				return func(Face(d));
			return true;
		});
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
	{
		Dart d = c.dart;
		do
		{
			if (!func(Face(d)))
				break;
			d = phi3(m, phi2(m, d));
		} while (d != c.dart);
	}
	else
	{
		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Face>(m))
		{
			CellMarkerStore<MESH, Face> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Face f(d);
				if constexpr (mesh_traits<MESH>::dimension == 2) // faces can be boundary cells
				{
					if (!marker.is_marked(f) && !is_boundary(m, d))
					{
						marker.mark(f);
						return func(f);
					}
				}
				else
				{
					if (!marker.is_marked(f))
					{
						marker.mark(f);
						return func(f);
					}
				}
				return true;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				if constexpr (mesh_traits<MESH>::dimension == 2) // faces can be boundary cells
				{
					if (!is_boundary(m, d) && !marker.is_marked(d))
					{
						Face f(d);
						foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(f);
					}
				}
				else
				{
					if (!marker.is_marked(d))
					{
						Face f(d);
						foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(f);
					}
				}
				return true;
			});
		}
	}
}

/*****************************************************************************/

// template <typename MESH, typename FUNC>
// void foreach_adjacent_face_through_edge(MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_adjacent_face_through_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_adjacent_face_through_edge(m, f, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_face_through_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func,
										CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
			if (!is_boundary(m, d))
				return func(Face(phi2(m, d)));
			return true;
		});
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		using Face2 = typename mesh_traits<MESH>::Face2;
		using Edge = typename mesh_traits<MESH>::Edge;

		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Face>(m))
		{
			CellMarkerStore<MESH, Face> marker(m);
			marker.mark(f);
			foreach_dart_of_orbit(m, Face2(f.dart), [&](Dart d) -> bool {
				bool cont = true;
				foreach_incident_face(m, Edge(d), [&](Face iface) -> bool {
					if (!marker.is_marked(iface))
					{
						cont = func(iface);
						marker.mark(iface);
						return cont;
					}
					return true;
				});
				return cont;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
				marker.mark(d);
				return true;
			});
			foreach_dart_of_orbit(m, Face2(f.dart), [&](Dart d) -> bool {
				bool cont = true;
				foreach_incident_face(m, Edge(d), [&](Face iface) -> bool {
					if (!marker.is_marked(iface.dart))
					{
						cont = func(iface);
						foreach_dart_of_orbit(m, iface, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return cont;
					}
					return true;
				});
				return cont;
			});
		}
	}
}

//////////////////////
/// IncidenceGraph ///
//////////////////////

template <typename CELL, typename FUNC>
auto foreach_incident_face(const IncidenceGraph& ig, CELL c, const FUNC& func)
{
	using Face = mesh_traits<IncidenceGraph>::Face;

	static_assert(is_in_tuple<CELL, mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this IncidenceGraph");
	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, mesh_traits<IncidenceGraph>::Vertex>)
	{
		CellMarkerStore<IncidenceGraph, Face> marker(ig);
		for (auto& ep : (*ig.vertex_incident_edges_)[c.index_])
		{
			bool stop = false;
			for (auto& fp : (*ig.edge_incident_faces_)[ep.index_])
			{
				stop = !func(fp);
				if (stop)
					break;
			}
			if (stop)
				break;
		}
	}
	else if constexpr (std::is_same_v<CELL, mesh_traits<IncidenceGraph>::Edge>)
	{
		for (auto& fp : (*ig.edge_incident_faces_)[c.index_])
		{
			if (!func(fp))
				break;
		}
	}
}

template <typename FUNC>
auto foreach_adjacent_face_through_edge(const IncidenceGraph& ig, IncidenceGraph::Face f, const FUNC& func)
{
	using Face = typename mesh_traits<IncidenceGraph>::Face;

	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	bool stop = false;
	CellMarkerStore<IncidenceGraph, IncidenceGraph::Face> marker(ig);
	marker.mark(f);
	foreach_incident_edge(ig, f, [&](IncidenceGraph::Edge e) -> bool {
		foreach_incident_face(ig, e, [&](IncidenceGraph::Face f1) -> bool {
			if(!marker.is_marked(f1))
			{
				marker.mark(f1);
				stop = func(f1);
			}

			return !stop;
		});
		return !stop;
	});
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Face> incident_faces(MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::Face> incident_faces(const MESH& m, CELL c)
{
	using Face = typename mesh_traits<MESH>::Face;
	std::vector<Face> faces;
	faces.reserve(32u);
	foreach_incident_face(m, c, [&](Face f) -> bool {
		faces.push_back(f);
		return true;
	});
	return faces;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_
