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

#ifndef CGOGN_CORE_TYPES_MAPS_GMAP_GMAP_BASE_H_
#define CGOGN_CORE_TYPES_MAPS_GMAP_GMAP_BASE_H_

#include <cgogn/core/types/maps/map_base.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT GMapBase : public MapBase
{
};

struct GMap0;
struct GMap1;
struct GMap2;
struct GMap3;

/*************************************************************************/
// Variable arguments beta function
/*************************************************************************/

template <int8 Arg, int8... Args, typename MESH>
auto beta(const MESH& m, Dart d) -> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>, Dart>
{
	static_assert((Arg >= 0 && Arg <= mesh_traits<MESH>::dimension), "Bad beta value");

	Dart res;
	if constexpr (Arg == 0)
		res = beta0(m, d);
	if constexpr (Arg == 1)
		res = beta1(m, d);
	if constexpr (Arg == 2)
		res = beta2(m, d);
	if constexpr (Arg == 3)
		res = beta3(m, d);

	if constexpr (sizeof...(Args) > 0)
		return beta<Args...>(m, res);
	else
		return res;
}

template <int8 Arg, int8... Args, typename MESH>
auto phi(const MESH& m, Dart d) -> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>, Dart>
{
	static_assert((Arg >= -1 && Arg <= mesh_traits<MESH>::dimension), "Bad phi value");

	Dart res;
	if constexpr (Arg == -1)
		res = phi_1(m, d);
	if constexpr (Arg == 1)
		res = phi1(m, d);
	if constexpr (Arg == 2)
		res = phi2(m, d);
	if constexpr (Arg == 3)
		res = phi3(m, d);

	if constexpr (sizeof...(Args) > 0)
		return phi<Args...>(m, res);
	else
		return res;
}

/*************************************************************************/
// Orbit traversals
/*************************************************************************/

// do not downgrade the concrete MESH type to CMapBase because some phi functions can be specialized for derived MESH
// types

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 0)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (f(d))
		f(beta0(m, d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA1(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (f(d))
		f(beta1(m, d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		it = beta0(m, it);
		if (!f(it))
			return;
		it = beta1(m, it);
	} while ((it != d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA01(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = d;
	do
	{
		if (!f(it))
			return;
		it = beta<0, 1>(m, it);
	} while (it != d);
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA21(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		it = beta<2, 1>(m, it);
	} while ((it != d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA1_BETA2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 2)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		it = beta1(m, it);
		if (!f(it))
			return;
		it = beta2(m, it);
	} while ((it != d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 2)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	if (!f(it))
		return;
	it = beta0(m, it);

	if (!f(it))
		return;
	it = beta2(m, it);

	if (!f(it))
		return;
	it = beta0(m, it);

	f(it);
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA2_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension == 3)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		if (!f(beta0(m, it)))
			return;
		it = beta2(m, it);
		if (!f(it))
			return;
		if (!f(beta0(m, it)))
			return;
		it = beta3(m, it);
	} while ((it != d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1_BETA2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 2)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	// NON OPTIMIZED (BUT SO SIMPLE TO CODE!)
	// TODO: optimized version

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		if (!f(curr_dart))
			break;

		const Dart d0 = beta0(m, curr_dart);
		const Dart d1 = beta1(m, curr_dart);
		const Dart d2 = beta2(m, curr_dart);

		if (!marker.is_marked(d0))
			marker.mark(d0);
		if (!marker.is_marked(d1))
			marker.mark(d1);
		if (!marker.is_marked(d2))
			marker.mark(d2);
	}
	// DartMarkerStore<MESH> marker(m);

	// std::vector<Dart> visited_faces;
	// visited_faces.reserve(32);
	// visited_faces.push_back(d); // Start with the face of d

	//// For every face added to the list
	// for (uint32 i = 0; i < uint32(visited_faces.size()); ++i)
	//{
	//	const Dart e = visited_faces[i];
	//	if (!marker.is_marked(e)) // Face has not been visited yet
	//	{
	//		foreach_dart_of_orbit(m, GMap1::Face(d), [&](Dart fd) -> bool
	//		{
	//			if (!f(fd))
	//				return false;
	//			marker.mark(fd);				  // Mark
	//			const Dart adj = beta2(m, fd);	  // Get adjacent face
	//			if (!marker.is_marked(adj))		  // no need to test FP because of marker
	//				visited_faces.push_back(adj); // Add it
	//			return true;
	//		});
	//	}
	// }
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension == 3)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	foreach_dart_of_BETA0_BETA1(m, d, [&](Dart fd) -> bool {
		if (f(fd))
		{
			Dart e = beta3(m, fd);
			if (e != fd)
				return f(e);
		}
		return false;
	});
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA1_BETA2_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		if (!f(curr_dart))
			break;

		const Dart d1 = beta1(m, curr_dart);
		const Dart d2 = beta2(m, curr_dart);
		const Dart d3 = beta3(m, curr_dart);

		if (!marker.is_marked(d1))
			marker.mark(d1);
		if (!marker.is_marked(d2))
			marker.mark(d2);
		if (!marker.is_marked(d3))
			marker.mark(d3);
	}
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1_BETA2_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	// NON OPTIMIZED (BUT SO SIMPLE TO CODE!)
	// TODO: optimized version

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		if (!f(curr_dart))
			break;

		const Dart d0 = beta0(m, curr_dart);
		const Dart d1 = beta1(m, curr_dart);
		const Dart d2 = beta2(m, curr_dart);
		const Dart d3 = beta3(m, curr_dart);

		if (!marker.is_marked(d0))
			marker.mark(d0);
		if (!marker.is_marked(d1))
			marker.mark(d1);
		if (!marker.is_marked(d2))
			marker.mark(d2);
		if (!marker.is_marked(d3))
			marker.mark(d3);
	}
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_dart_of_orbit(const MESH& m, CELL c, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	static const Orbit orbit = CELL::ORBIT;

	if constexpr (orbit == DART)
	{
		unused_parameters(m);
		f(c.dart_);
		return;
	}
	if constexpr (orbit == BETA0)
	{
		foreach_dart_of_BETA0(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA1)
	{
		foreach_dart_of_BETA1(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1)
	{
		foreach_dart_of_BETA0_BETA1(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA2)
	{
		foreach_dart_of_BETA0_BETA2(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA1_BETA2)
	{
		foreach_dart_of_BETA1_BETA2(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1_BETA2)
	{
		foreach_dart_of_BETA0_BETA1_BETA2(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1_BETA3)
	{
		foreach_dart_of_BETA0_BETA1_BETA3(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA2_BETA3)
	{
		foreach_dart_of_BETA0_BETA2_BETA3(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA1_BETA2_BETA3)
	{
		foreach_dart_of_BETA1_BETA2_BETA3(m, c.dart_, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1_BETA2_BETA3)
	{
		foreach_dart_of_BETA0_BETA1_BETA2_BETA3(m, c.dart_, f);
		return;
	}
}

/*************************************************************************/
// Local vertex traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_incident_vertex(m, c, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& m, CELL c, const FUNC& func, MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap1&> && mesh_traits<MESH>::dimension == 1 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_BETA01(m, c.dart_, [&](Dart d) -> bool { return func(Vertex(d)); });
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap2&> && (mesh_traits<MESH>::dimension == 2))
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
		{
			if (func(Vertex(c.dart_)))
				func(Vertex(phi2(m, c.dart_)));
			return;
		}
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge>)
		{
			// TODO
			func(Vertex(c.dart_));
			func(Vertex(beta0(m, c.dart_)));
			// CHECK
			return;
		}
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_BETA01(m, c.dart_, [&](Dart d) -> bool { return func(Vertex(d)); });
			return;
		}
	}

	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3 &&
				  (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge> ||
				   std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge>))
	{
		if (func(Vertex(c.dart_)))
			func(Vertex(beta0(m, c.dart_)));
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_BETA01(m, c.dart_, [&](Dart d) -> bool { return func(Vertex(d)); });
		return;
	}
	if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<Vertex>(m))
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

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_adjacent_vertex_through_edge(m, v, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func,
										  MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		foreach_dart_of_BETA21(m, v.dart_, [&](Dart d) -> bool { return func(Vertex(phi2(m, d))); });
	}
	else if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<Vertex>(m))
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
				if (!marker.is_marked(av.dart_))
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

/*************************************************************************/
// Local half-edge traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_halfedge(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, HalfEdge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(HalfEdge(d)); });
}

/*************************************************************************/
// Local edge traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_edge(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_incident_edge(m, c, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_edge(const MESH& m, CELL c, const FUNC& func, MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Edge = typename mesh_traits<MESH>::Edge;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap1&> && mesh_traits<MESH>::dimension == 1 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		foreach_dart_of_BETA01(m, c.dart_, [&](Dart d) -> bool { return func(Edge(d)); });
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex>)
		{
			foreach_dart_of_BETA21(m, c.dart_, [&](Dart d) -> bool { return func(Edge(d)); });
			return;
		}
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge>)
		{
			func(Edge(c.dart_));
			func(Edge(phi2(m, c.dart_)));
			// foreach_dart_of_orbit(m, c, [&](Dart d) -> bool { return func(Edge(d)); });
			return;
		}
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
		{
			foreach_dart_of_BETA01(m, c.dart_, [&](Dart d) -> bool { return func(Edge(d)); });
			return;
		}
	}

	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3 &&
				  (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>))
	{
		foreach_dart_of_BETA01(m, c.dart_, [&](Dart d) -> bool { return func(Edge(d)); });
		return;
	}

	if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<Edge>(m))
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

template <typename MESH, typename FUNC>
auto foreach_adjacent_edge_through_face(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_adjacent_edge_through_face(m, e, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_edge_through_face(const MESH& m, typename mesh_traits<MESH>::Edge e, const FUNC& func,
										MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Edge = typename mesh_traits<MESH>::Edge;

	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		using Face = typename mesh_traits<MESH>::Face;

		Dart d1 = e.dart_;
		Dart d2 = phi2(m, d1);
		if (!is_boundary(m, d1))
			foreach_dart_of_BETA01(m, d1, [&](Dart d) -> bool { return d != d1 ? func(Edge(d)) : true; });
		if (!is_boundary(m, d2))
			foreach_dart_of_BETA01(m, d2, [&](Dart d) -> bool { return d != d2 ? func(Edge(d)) : true; });
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		using Face2 = typename mesh_traits<MESH>::Face2;

		foreach_dart_of_BETAs(m, e, beta<2, 3>, [&](Dart ed) -> bool {
			foreach_dart_of_BETA01(m, ed, [&](Dart d) -> bool { return d != ed ? func(Edge(d)) : true; });
			return true;
		});
	}
}

/*************************************************************************/
// Local face traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_face(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_incident_face(m, c, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_face(const MESH& m, CELL c, const FUNC& func, MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex>)
		{
			foreach_dart_of_BETA21(m, c.dart_, [&](Dart d) -> bool {
				if (!is_boundary(m, d))
					return func(Face(d));
				return true;
			});
			return;
		}

		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::HalfEdge>)
		{
			if (!is_boundary(m, c.dart_))
				func(Face(c.dart_));
			return;
		}

		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
		{
			if (!is_boundary(m, c.dart_))
				if (!func(Face(c.dart_)))
					return;
			Dart d2 = beta<2, 0>(m, c.dart_);
			if (!is_boundary(m, d2))
				func(Face(d2));
			return;
		}
	}
	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
	{
		Dart d = c.dart_;
		do
		{
			if (!func(Face(d)))
				break;
			d = phi3(m, phi2(m, d));
		} while (d != c.dart_);
		return;
	}

	if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<Face>(m))
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

template <typename MESH, typename FUNC>
auto foreach_adjacent_face_through_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_adjacent_face_through_edge(m, f, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_face_through_edge(const MESH& m, typename mesh_traits<MESH>::Face f, const FUNC& func,
										MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		foreach_dart_of_BETA01(m, f.dart_, [&](Dart d) -> bool {
			Dart d2 = phi2(m, d);
			if (!is_boundary(m, d2))
				return func(Face(d2));
			return true;
		});
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3)
	{
		using Face2 = typename mesh_traits<MESH>::Face2;
		using Edge = typename mesh_traits<MESH>::Edge;

		if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<Face>(m))
		{
			CellMarkerStore<MESH, Face> marker(m);
			marker.mark(f);
			foreach_dart_of_orbit(m, Face2(f.dart_), [&](Dart d) -> bool {
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
			foreach_dart_of_orbit(m, Face2(f.dart_), [&](Dart d) -> bool {
				bool cont = true;
				foreach_incident_face(m, Edge(d), [&](Face iface) -> bool {
					if (!marker.is_marked(iface.dart_))
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

/*************************************************************************/
// Local volume traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	foreach_incident_volume(m, c, func, MapBase::TraversalPolicy::AUTO);
}

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func, MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	using Volume = typename mesh_traits<MESH>::Volume;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, GMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		func(Volume(c.dart_));
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
	{
		Dart d = c.dart_;
		do
		{
			if (!is_boundary(m, d))
			{
				if (!func(Volume(d)))
					break;
			}
			d = phi3(m, phi2(m, d));
		} while (d != c.dart_);
		return;
	}
	if constexpr (std::is_convertible_v<MESH&, GMap3&> && mesh_traits<MESH>::dimension == 3 &&
				  std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		Dart d = c.dart_;
		if (!is_boundary(m, d))
			if (!func(Volume(d)))
				return;
		d = phi3(m, d);
		if (!is_boundary(m, d))
			func(Volume(d));
	}

	if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<Volume>(m))
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

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MAPS_GMAP_GMAP_BASE_H_
