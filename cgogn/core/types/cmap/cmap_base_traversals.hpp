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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_TRAVERSALS_HPP_
#define CGOGN_CORE_CMAP_CMAP_BASE_TRAVERSALS_HPP_


#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/utils/thread.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/core/types/cmap/phi.h>


namespace cgogn
{

struct CMap1;
struct CMap2;
struct CMap3;
struct Graph;

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

	if constexpr (std::is_convertible_v<MESH&, Graph&> && mesh_traits<MESH>::dimension == 1)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex>)
		{
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
		}
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap1&> && mesh_traits<MESH>::dimension == 1)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
		}
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
		}
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Face2(c.dart),
								[&](Dart d) -> bool { return func(Edge(d)); });
		}
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
			Dart d2 = phi2(m, d);
			if (!is_boundary(m, d2))
				return func(Face(d2));
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


template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_halfedge(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, HalfEdge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(HalfEdge(d)); });
}




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

	if constexpr (std::is_convertible_v<MESH&, Graph&> && mesh_traits<MESH>::dimension == 1)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
		{
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Vertex(d)); });
		}
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap1&> && mesh_traits<MESH>::dimension == 1)
	{

		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Vertex(d)); });
		}
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2 )
		{
			if constexpr  (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
			{
				foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Vertex(d)); });
			}
		}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
		{
			if constexpr  (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge> ||
						std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge>)
			{
				foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Edge2(c.dart),
									[&](Dart d) -> bool { return func(Vertex(d)); });
			}
		}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Face2(c.dart),
								[&](Dart d) -> bool { return func(Vertex(d)); });
		}
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


template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_incident_volume(m, c, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func, CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Volume = typename mesh_traits<MESH>::Volume;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		func(Volume(c.dart));
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
		{
			Dart d = c.dart;
			do
			{
				if (!is_boundary(m, d))
				{
					if (!func(Volume(d)))
						break;
				}
				d = phi3(m, phi2(m, d));
			} while (d != c.dart);
		}
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			Dart d = c.dart;
			if (!is_boundary(m, d))
				if (!func(Volume(d)))
					return;
			d = phi3(m, d);
			if (!is_boundary(m, d))
				func(Volume(d));
		}
	}
	else
	{
		if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<Volume>(m))
		{
			CellMarkerStore<MESH, Volume> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Volume v(d);
				if constexpr (mesh_traits<MESH>::dimension == 3) // volumes can be boundary cells
				{
					if (!is_boundary(m, d) && !marker.is_marked(v))
					{
						marker.mark(v);
						return func(v);
					}
				}
				else
				{
					if (!marker.is_marked(v))
					{
						marker.mark(v);
						return func(v);
					}
				}
				return true;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				if constexpr (mesh_traits<MESH>::dimension == 3) // volumes can be boundary cells
				{
					if (!is_boundary(m, d) && !marker.is_marked(d))
					{
						Volume v(d);
						foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(v);
					}
				}
				else
				{
					if (!marker.is_marked(d))
					{
						Volume v(d);
						foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(v);
					}
				}
				return true;
			});
		}
	}
}

template <typename MESH, typename FUNC>
auto foreach_cell(const MESH& m, const FUNC& func) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	foreach_cell(m, func, CMapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_cell(const MESH& m, const FUNC& f, CMapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, CELL>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (traversal_policy == CMapBase::TraversalPolicy::AUTO && is_indexed<CELL>(m))
	{
		CellMarker<MESH, CELL> cm(m);
		for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		{
			const CELL c(d);
			if (!is_boundary(m, d) && !cm.is_marked(c))
			{
				cm.mark(c);
				if (!f(c))
					break;
			}
		}
	}
	else
	{
		DartMarker dm(m);
		for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		{
			if (!is_boundary(m, d) && !dm.is_marked(d))
			{
				const CELL c(d);
				foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
					dm.mark(d);
					return true;
				});
				if (!f(c))
					break;
			}
		}
	}
}

template <typename MESH, typename FUNC>
auto parallel_foreach_cell(const MESH& m, const FUNC& f) 
-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, CELL>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	ThreadPool* pool = thread_pool();
	uint32 nb_workers = pool->nb_workers();
	if (nb_workers == 0)
		return foreach_cell(m, f);

	using VecCell = std::vector<uint32>;
	using Future = std::future<void>;

	std::array<std::vector<VecCell*>, 2> cells_buffers;
	std::array<std::vector<Future>, 2> futures;
	cells_buffers[0].reserve(nb_workers);
	cells_buffers[1].reserve(nb_workers);
	futures[0].reserve(nb_workers);
	futures[1].reserve(nb_workers);

	Buffers<uint32>* buffers = uint32_buffers();

	Dart it = m.begin();
	Dart last = m.end();

	uint32 i = 0u; // buffer id (0/1)
	uint32 j = 0u; // thread id (0..nb_workers)

	if (is_indexed<CELL>(m))
	{
		CellMarker<MESH, CELL> cm(m);
		while (it.index < last.index)
		{
			// fill buffer
			cells_buffers[i].push_back(buffers->buffer());
			VecCell& cells = *cells_buffers[i].back();
			cells.reserve(PARALLEL_BUFFER_SIZE);
			for (uint32 k = 0u; k < PARALLEL_BUFFER_SIZE && it.index < last.index;)
			{
				CELL c(it);
				if (!is_boundary(m, it) && !cm.is_marked(c))
				{
					cm.mark(c);
					cells.push_back(c.dart.index);
					++k;
				}
				it = m.next(it);
			}
			// launch thread
			futures[i].push_back(pool->enqueue([&cells, &f]() {
				for (uint32 index : cells)
					f(CELL(Dart(index)));
			}));
			// next thread
			if (++j == nb_workers)
			{ // again from 0 & change buffer
				j = 0u;
				i = (i + 1u) % 2u;
				for (auto& fu : futures[i])
					fu.wait();
				for (auto& b : cells_buffers[i])
					buffers->release_buffer(b);
				futures[i].clear();
				cells_buffers[i].clear();
			}
		}
	}
	else
	{
		DartMarker dm(m);
		while (it.index < last.index)
		{
			// fill buffer
			cells_buffers[i].push_back(buffers->buffer());
			VecCell& cells = *cells_buffers[i].back();
			cells.reserve(PARALLEL_BUFFER_SIZE);
			for (uint32 k = 0u; k < PARALLEL_BUFFER_SIZE && it.index < last.index;)
			{
				if (!is_boundary(m, it) && !dm.is_marked(it))
				{
					CELL c(it);
					foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
						dm.mark(d);
						return true;
					});
					cells.push_back(c.dart.index);
					++k;
				}
				it = m.next(it);
			}
			// launch thread
			futures[i].push_back(pool->enqueue([&cells, &f]() {
				for (uint32 index : cells)
					f(CELL(Dart(index)));
			}));
			// next thread
			if (++j == nb_workers)
			{ // again from 0 & change buffer
				j = 0u;
				i = (i + 1u) % 2u;
				for (auto& fu : futures[i])
					fu.wait();
				for (auto& b : cells_buffers[i])
					buffers->release_buffer(b);
				futures[i].clear();
				cells_buffers[i].clear();
			}
		}
	}

	// clean all at the end
	for (auto& fu : futures[0u])
		fu.wait();
	for (auto& b : cells_buffers[0u])
		buffers->release_buffer(b);
	for (auto& fu : futures[1u])
		fu.wait();
	for (auto& b : cells_buffers[1u])
		buffers->release_buffer(b);
}


} // namespace cgogn


#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
