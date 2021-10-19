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

#ifndef CGOGN_CORE_FUNCTIONS_CONVERT_H_
#define CGOGN_CORE_FUNCTIONS_CONVERT_H_

#include <cgogn/core/functions/attributes.h>

#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/volume.h>

namespace cgogn
{

template <typename SURFACE, typename NONMANIFOLD, typename FUNC1, typename FUNC2, typename FUNC3>
void non_manifold_from_surface(SURFACE& s, NONMANIFOLD& nm, const FUNC1& on_vertex_added, const FUNC2& on_edge_added,
							   const FUNC3& on_face_added)
{
	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NonManifoldFace = typename mesh_traits<NONMANIFOLD>::Face;

	static_assert(is_ith_func_parameter_same<FUNC1, 0, NonManifoldVertex>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC1, 1, SurfaceVertex>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC2, 0, NonManifoldEdge>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC2, 1, SurfaceEdge>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC3, 0, NonManifoldFace>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC3, 1, SurfaceFace>::value, "Wrong function parameter type");

	auto non_manifold_vertex = add_attribute<NonManifoldVertex, SurfaceVertex>(s, "__non_manifold_vertex");
	foreach_cell(s, [&](SurfaceVertex sv) -> bool {
		NonManifoldVertex nmv = add_vertex(nm);
		value<NonManifoldVertex>(s, non_manifold_vertex, sv) = nmv;
		on_vertex_added(nmv, sv);
		return true;
	});
	auto non_manifold_edge = add_attribute<NonManifoldEdge, SurfaceEdge>(s, "__non_manifold_edge");
	foreach_cell(s, [&](SurfaceEdge se) -> bool {
		std::vector<SurfaceVertex> iv = incident_vertices(s, se);
		NonManifoldEdge nme = add_edge(nm, value<NonManifoldVertex>(s, non_manifold_vertex, iv[0]),
									   value<NonManifoldVertex>(s, non_manifold_vertex, iv[1]));
		value<NonManifoldEdge>(s, non_manifold_edge, se) = nme;
		on_edge_added(nme, se);
		return true;
	});
	foreach_cell(s, [&](SurfaceFace sf) -> bool {
		std::vector<SurfaceEdge> ie = incident_edges(s, sf);
		std::vector<NonManifoldEdge> edges;
		std::transform(ie.begin(), ie.end(), std::back_inserter(edges),
					   [&](SurfaceEdge e) { return value<NonManifoldEdge>(s, non_manifold_edge, e); });
		NonManifoldFace nmf = add_face(nm, edges);
		on_face_added(nmf, sf);
		return true;
	});
	remove_attribute<SurfaceVertex>(s, non_manifold_vertex);
	remove_attribute<SurfaceEdge>(s, non_manifold_edge);
}

template <typename SURFACE, typename GRAPH, typename FUNC1, typename FUNC2>
void graph_from_surface(SURFACE& s, GRAPH& g, const FUNC1& on_vertex_added, const FUNC2& on_edge_added)
{
	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

	static_assert(is_ith_func_parameter_same<FUNC1, 0, GraphVertex>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC1, 1, SurfaceVertex>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC2, 0, GraphEdge>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC2, 1, SurfaceEdge>::value, "Wrong function parameter type");

	auto graph_vertex = add_attribute<GraphVertex, SurfaceVertex>(s, "__graph_vertex");
	foreach_cell(s, [&](SurfaceVertex sv) -> bool {
		GraphVertex gv = add_vertex(g);
		value<GraphVertex>(s, graph_vertex, sv) = gv;
		on_vertex_added(gv, sv);
		return true;
	});
	foreach_cell(s, [&](SurfaceEdge se) -> bool {
		std::vector<SurfaceVertex> iv = incident_vertices(s, se);
		GraphEdge ge =
			connect_vertices(g, value<GraphVertex>(s, graph_vertex, iv[0]), value<GraphVertex>(s, graph_vertex, iv[1]));
		on_edge_added(ge, se);
		return true;
	});
	remove_attribute<SurfaceVertex>(s, graph_vertex);
}

template <typename NONMANIFOLD, typename GRAPH, typename FUNC1, typename FUNC2>
void graph_from_non_manifold(NONMANIFOLD& nm, GRAPH& g, const FUNC1& on_vertex_added, const FUNC2& on_edge_added)
{
	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

	static_assert(is_ith_func_parameter_same<FUNC1, 0, GraphVertex>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC1, 1, NonManifoldVertex>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC2, 0, GraphEdge>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC2, 1, NonManifoldEdge>::value, "Wrong function parameter type");

	auto graph_vertex = add_attribute<GraphVertex, NonManifoldVertex>(nm, "__graph_vertex");
	foreach_cell(nm, [&](NonManifoldVertex nmv) -> bool {
		GraphVertex gv = add_vertex(g);
		value<GraphVertex>(nm, graph_vertex, nmv) = gv;
		on_vertex_added(gv, nmv);
		return true;
	});
	foreach_cell(nm, [&](NonManifoldEdge nme) -> bool {
		std::vector<NonManifoldVertex> iv = incident_vertices(nm, nme);
		GraphEdge ge = connect_vertices(g, value<GraphVertex>(nm, graph_vertex, iv[0]),
										value<GraphVertex>(nm, graph_vertex, iv[1]));
		on_edge_added(ge, nme);
		return true;
	});
	remove_attribute<NonManifoldVertex>(nm, graph_vertex);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_CONVERT_H_
