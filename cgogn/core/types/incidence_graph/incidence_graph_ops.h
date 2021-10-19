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

#ifndef CGOGN_CORE_INCIDENCE_GRAPH_OPS_H_
#define CGOGN_CORE_INCIDENCE_GRAPH_OPS_H_

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
namespace cgogn
{

template <typename CELL>
CELL add_cell(IncidenceGraph& ig)
{
	return CELL(ig.attribute_containers_[CELL::CELL_INDEX].new_index());
}

template <typename CELL>
void remove_cell(IncidenceGraph& ig, CELL c)
{
	ig.attribute_containers_[CELL::CELL_INDEX].release_index(c.index_);
}

inline bool sort_face_edges(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<uint8>& edges_dir = (*ig.face_incident_edges_dir_)[f.index_];

	std::vector<Edge> unordered_edges;
	unordered_edges.swap(edges);

	edges.push_back(unordered_edges.front());
	edges_dir.push_back(0);
	unordered_edges.erase(unordered_edges.begin());
	uint32 vid0 = (*ig.edge_incident_vertices_)[edges.front().index_].first.index_;
	uint32 vid1 = (*ig.edge_incident_vertices_)[edges.front().index_].second.index_;
	bool broken = false;

	while (unordered_edges.size() > 0)
	{
		uint32 i, end;
		for (i = 0, end = unordered_edges.size(); i < unordered_edges.size(); ++i)
		{
			Edge e = unordered_edges[i];
			std::pair<Vertex, Vertex>& evs = (*ig.edge_incident_vertices_)[e.index_];
			uint32 ev = (evs.first.index_ == vid1 ? evs.second.index_
												  : (evs.second.index_ == vid1 ? evs.first.index_ : INVALID_INDEX));

			if (ev != INVALID_INDEX)
			{
				vid1 = ev;
				edges.push_back(e);
				edges_dir.push_back(ev == evs.first.index_ ? 1 : 0);
				unordered_edges.erase(unordered_edges.begin() + i);
				break;
			}
		}

		broken = (vid1 == vid0) || (i == end);
		if (broken)
			break;
	}

	return (broken && unordered_edges.size() == 0);
}

inline bool same_edge(IncidenceGraph& ig, IncidenceGraph::Edge e1, IncidenceGraph::Edge e2)
{
	using Vertex = IncidenceGraph::Vertex;

	Vertex e1v1 = (*ig.edge_incident_vertices_)[e1.index_].first;
	Vertex e1v2 = (*ig.edge_incident_vertices_)[e1.index_].second;
	Vertex e2v1 = (*ig.edge_incident_vertices_)[e2.index_].first;
	Vertex e2v2 = (*ig.edge_incident_vertices_)[e2.index_].second;

	return (e1v1 == e2v1 && e1v2 == e2v2) || (e1v1 == e2v2 && e1v2 == e2v1);
}

inline void remove_edge_in_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v, IncidenceGraph::Edge edge_to_remove)
{
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.vertex_incident_edges_)[v.index_];
	auto eit = std::find(edges.begin(), edges.end(), edge_to_remove);
	if (eit != edges.end())
		edges.erase(eit);
}

inline void remove_face_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Face face_to_remove)
{
	using Face = IncidenceGraph::Face;

	std::vector<Face>& faces = (*ig.edge_incident_faces_)[e.index_];
	auto fit = std::find(faces.begin(), faces.end(), face_to_remove);
	if (fit != faces.end())
		faces.erase(fit);
}

inline void remove_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge edge_to_remove)
{
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<uint8>& edges_dir = (*ig.face_incident_edges_dir_)[f.index_];
	auto eit = edges.begin();
	auto edit = edges_dir.begin();
	for (Edge e : edges)
	{
		if (e == edge_to_remove)
			break;
		eit++;
		edit++;
	}
	if (eit != edges.end())
	{
		edges.erase(eit);
		edges_dir.erase(edit);
	}
}

inline void replace_vertex_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Vertex old_vertex,
								   IncidenceGraph::Vertex new_vertex)
{
	if ((*ig.edge_incident_vertices_)[e.index_].first == old_vertex)
		(*ig.edge_incident_vertices_)[e.index_].first = new_vertex;
	if ((*ig.edge_incident_vertices_)[e.index_].second == old_vertex)
		(*ig.edge_incident_vertices_)[e.index_].second = new_vertex;
}

inline void replace_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge old_edge,
								 IncidenceGraph::Edge new_edge)
{
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	auto eit = std::find(edges.begin(), edges.end(), old_edge);
	if (eit != edges.end())
		*eit = new_edge;
}

inline IncidenceGraph::Vertex common_vertex(IncidenceGraph& ig, IncidenceGraph::Edge e0, IncidenceGraph::Edge e1)
{
	using Vertex = IncidenceGraph::Vertex;

	Vertex v0 = (*ig.edge_incident_vertices_)[e0.index_].first;
	Vertex v1 = (*ig.edge_incident_vertices_)[e0.index_].second;
	Vertex v2 = (*ig.edge_incident_vertices_)[e1.index_].first;
	Vertex v3 = (*ig.edge_incident_vertices_)[e1.index_].second;

	if (v0 == v2 || v0 == v3)
		return v0;

	if (v1 == v2 || v1 == v3)
		return v1;

	return Vertex();
}

inline std::vector<IncidenceGraph::Vertex> sorted_face_vertices(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<Vertex> sorted_vertices;
	for (uint32 i = 0; i < edges.size(); ++i)
		sorted_vertices.push_back(common_vertex(ig, edges[i], edges[(i + 1) % edges.size()]));
	sorted_vertices.insert(sorted_vertices.begin(), sorted_vertices.back());
	sorted_vertices.pop_back();
	return sorted_vertices;
}

inline std::pair<uint32, uint32> pseudoDegree(const IncidenceGraph& ig, IncidenceGraph::Vertex v)
{
	std::pair<uint32, uint32> info;
	info.first = 0; // incident leaflets + isolated edges
	info.second = 0; // incident isolated branches  
	// uint32 degree = 0;

	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
		uint32 count = 0;
		foreach_incident_face(ig, e, [&](IncidenceGraph::Face f) -> bool {
			++count;
			return true;
		});
		switch(count)
		{
			case 0: 
				info.first += 2;
				++info.second;
				break;
			case 1: 
				info.first += 1;
				break;
			case 2:
				break;
			default:
				info.first = INVALID_INDEX;
				break;
		}

		return (info.first != INVALID_INDEX);
	});
	info.first /= 2;
	return info;
}

} // namespace cgogn

#endif // CGOGN_CORE_INCIDENCE_GRAPH_OPS_H_
