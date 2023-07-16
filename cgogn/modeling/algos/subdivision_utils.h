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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_UTILS_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_UTILS_H_

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/dart.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <set>
#include <unordered_map>

namespace cgogn
{

// forward for SFINAE
struct MapBase;

namespace modeling
{

///////////////
// MapBase:2 //
///////////////

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&> &&
												   (mesh_traits<MESH>::dimension == 2)>* = nullptr>
void hexagon_to_triangles(MESH& m, typename mesh_traits<MESH>::Face f)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	cgogn_message_assert(codegree(m, f) == 6, "hexagon_to_triangles: given face should have 6 edges");
	Dart d0 = phi1(m, f.dart);
	Dart d1 = phi<1, 1>(m, d0);
	cut_face(m, Vertex(d0), Vertex(d1));
	Dart d2 = phi<1, 1>(m, d1);
	cut_face(m, Vertex(d1), Vertex(d2));
	Dart d3 = phi<1, 1>(m, d2);
	cut_face(m, Vertex(d2), Vertex(d3));
}

//////////////
// MapBase //
//////////////

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>* = nullptr>
typename mesh_traits<MESH>::Vertex quadrangulate_face(MESH& m, typename mesh_traits<MESH>::Face f)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	cgogn_message_assert(codegree(m, f) % 2 == 0, "quadrangulate_face: given face should have a pair codegree");

	Dart d0 = phi1(m, f.dart);
	Dart d1 = phi<1, 1>(m, d0);

	cut_face(m, Vertex(d0), Vertex(d1));
	cut_edge(m, Edge(phi_1(m, d0)));

	Dart x = phi<-1, 2>(m, d0);
	Dart dd = phi<1, 1, 1, 1>(m, x);
	while (dd != x)
	{
		Dart next = phi<1, 1>(m, dd);
		cut_face(m, Vertex(dd), Vertex(phi1(m, x)));
		dd = next;
	}

	return Vertex(phi2(m, x));
}

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>* = nullptr>
typename mesh_traits<MESH>::Vertex quadrangulate_face(
	MESH& m, typename mesh_traits<MESH>::Face f,
	const CellMarker<MESH, typename mesh_traits<MESH>::Vertex>& vertex_marker)
{
	return quadrangulate_face(m, f);
}

////////////////////
// IncidenceGraph //
////////////////////

// vertices added by edge cuts are marked
inline IncidenceGraph::Vertex quadrangulate_face(IncidenceGraph& ig, IncidenceGraph::Face f,
												 CellMarker<IncidenceGraph, IncidenceGraph::Vertex>& vertex_marker)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	cgogn_message_assert(codegree(ig, f) % 2 == 0, "quadrangulate_face: given face should have a pair codegree");

	Vertex v = add_vertex(ig);
	std::unordered_map<uint32, uint32> new_edges;
	std::set<Vertex> corner_vertices;

	const std::vector<Edge>& inc_edges = (*ig.face_incident_edges_)[f.index_];
	for (uint32 i = 0; i < inc_edges.size(); ++i)
	{
		Edge e1 = inc_edges[i];
		const std::pair<Vertex, Vertex>& inc_verts1 = (*ig.edge_incident_vertices_)[e1.index_];
		// get the edge vertex
		Vertex ev = vertex_marker.is_marked(inc_verts1.first) ? inc_verts1.first : inc_verts1.second;
		// create the edge of the edge vertex if needed
		if (new_edges.count(ev.index_) == 0)
		{
			Edge ne = add_edge(ig, v, ev);
			new_edges[ev.index_] = ne.index_;
		}
		// get the corner vertex
		Vertex cv = vertex_marker.is_marked(inc_verts1.first) ? inc_verts1.second : inc_verts1.first;
		corner_vertices.insert(cv);
	}
	for (Vertex cv : corner_vertices)
	{
		std::vector<Edge> edges;
		// find the two edges that are incident to cv
		for (Edge e : inc_edges)
		{
			const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[e.index_];
			if (inc_verts.first == cv || inc_verts.second == cv)
				edges.push_back(e);
		}
		// add the two new edges
		const std::pair<Vertex, Vertex>& edge0_inc_verts = (*ig.edge_incident_vertices_)[edges[0].index_];
		Vertex ev0 = vertex_marker.is_marked(edge0_inc_verts.first) ? edge0_inc_verts.first : edge0_inc_verts.second;
		const std::pair<Vertex, Vertex>& edge1_inc_verts = (*ig.edge_incident_vertices_)[edges[1].index_];
		Vertex ev1 = vertex_marker.is_marked(edge1_inc_verts.first) ? edge1_inc_verts.first : edge1_inc_verts.second;
		edges.push_back(Edge(new_edges[ev1.index_]));
		edges.push_back(Edge(new_edges[ev0.index_]));
		// create the face
		Face f = add_face(ig, edges);
		sort_face_edges(ig, f);
	}

	remove_face(ig, f);

	return v;
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_UTILS_H_
