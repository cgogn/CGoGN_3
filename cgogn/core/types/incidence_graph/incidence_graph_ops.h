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

inline IncidenceGraph::Vertex add_vertex(IncidenceGraph& ig)
{
	uint32 id = ig.attribute_containers_[IncidenceGraph::Vertex::CELL_INDEX].new_index();
	(*ig.vertices_)[id] = id;
	return IncidenceGraph::Vertex(id);
}

inline IncidenceGraph::Edge add_edge(IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1)
{
	uint32 id = ig.attribute_containers_[IncidenceGraph::Edge::CELL_INDEX].new_index();
	(*ig.edges_)[id] = id;
	IncidenceGraph::Edge e(id);
	(*ig.edge_incident_vertices_)[id] = {v0, v1};
	(*ig.vertex_incident_edges_)[v0.index_][id] = e;
	(*ig.vertex_incident_edges_)[v1.index_][id] = e;

	return e;
}

inline bool sort_edges(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge>& edges)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge> unordered_edges;
	unordered_edges.swap(edges);

	edges.push_back(unordered_edges.front());
	unordered_edges.erase(unordered_edges.begin());
	uint32 vid0 = (*ig.edge_incident_vertices_)[edges.front().index_].first.index_;
	uint32 vid1 = (*ig.edge_incident_vertices_)[edges.front().index_].second.index_;
	bool broken = false;

	while (unordered_edges.size())
	{
		uint32 i, end;
		for (i = 0, end = unordered_edges.size(); i < unordered_edges.size(); ++i)
		{
			Edge e = unordered_edges[i];
			std::pair<Vertex, Vertex> evs = (*ig.edge_incident_vertices_)[e.index_];
			uint32 ev = (evs.first.index_ == vid1 ? evs.second.index_
												  : (evs.second.index_ == vid1 ? evs.first.index_ : INVALID_INDEX));

			if (ev != INVALID_INDEX)
			{
				vid1 = ev;
				edges.push_back(e);
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

inline IncidenceGraph::Face add_face(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge> edges)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	if (sort_edges(ig, edges))
	{
		uint32 id = ig.attribute_containers_[Face::CELL_INDEX].new_index();
		(*ig.faces_)[id] = id;
		(*ig.face_incident_edges_)[id] = edges;
		Face f(id);
		for (Edge e : edges)
		{
			(*ig.edge_incident_faces_)[e.index_][id] = f;
		}
		(*ig.face_incident_edges_dir_)[id] = std::vector<uint32>(edges.size());
		
		for (uint32 i = 0; i < edges.size(); ++i)
		{
			std::pair<Vertex, Vertex> evs0 = (*ig.edge_incident_vertices_)[edges[i].index_];
			std::pair<Vertex, Vertex> evs1 = (*ig.edge_incident_vertices_)[edges[(i+1)%edges.size()].index_];
			if(evs0.first.index_ == evs1.first.index_ || evs0.first.index_ == evs1.second.index_)
				(*ig.face_incident_edges_dir_)[id][i] = 1;
			else
				(*ig.face_incident_edges_dir_)[id][i] = 0;
		}

		return f;
	}
	return Face();
}

inline void remove_face(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	if (f.is_valid())
	{
		std::vector<IncidenceGraph::Edge>& edges = (*ig.face_incident_edges_)[f.index_];
		for (IncidenceGraph::Edge e : edges)
		{
			(*ig.edge_incident_faces_)[e.index_].erase(f.index_);
		}
		for (IncidenceGraph::Edge e : edges)
		{
			std::cout << (*ig.edge_incident_faces_)[e.index_].size() << std::endl;
		}

		(*ig.faces_)[f.index_] = INVALID_INDEX;
		ig.attribute_containers_[IncidenceGraph::Face::CELL_INDEX].release_index(f.index_);
	}
}

inline void remove_edge(IncidenceGraph& ig, IncidenceGraph::Edge e)
{
	while ((*ig.edge_incident_faces_)[e.index_].size())
	{
		remove_face(ig, (*ig.edge_incident_faces_)[e.index_].begin()->second);
	}

	auto [v0, v1] = (*ig.edge_incident_vertices_)[e.index_];
	// std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex> evs = (*ig.edge_incident_vertices_)[e.index_];
	// (*ig.vertex_incident_edges_)[evs.first.index_].erase(e.index_);
	// (*ig.vertex_incident_edges_)[evs.second.index_].erase(e.index_);

	(*ig.vertex_incident_edges_)[v0.index_].erase(e.index_);
	(*ig.vertex_incident_edges_)[v1.index_].erase(e.index_);

	(*ig.edges_)[e.index_] = INVALID_INDEX;
	ig.attribute_containers_[IncidenceGraph::Edge::CELL_INDEX].release_index(e.index_);
}

inline void remove_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v)
{
	while ((*ig.vertex_incident_edges_)[v.index_].size())
	{
		remove_edge(ig, (*ig.vertex_incident_edges_)[v.index_].begin()->second);
	}
	ig.attribute_containers_[IncidenceGraph::Vertex::CELL_INDEX].release_index(v.index_);
	(*ig.vertices_)[v.index_] = INVALID_INDEX;
}

inline IncidenceGraph::Vertex common_vertex(IncidenceGraph& ig, IncidenceGraph::Edge e0, IncidenceGraph::Edge e1)
{
	uint32 vid0 = (*ig.edge_incident_vertices_)[e0.index_].first.index_;
	uint32 vid1 = (*ig.edge_incident_vertices_)[e0.index_].second.index_;
	uint32 vid2 = (*ig.edge_incident_vertices_)[e1.index_].first.index_;
	uint32 vid3 = (*ig.edge_incident_vertices_)[e1.index_].second.index_;

	if (vid0 == vid2 || vid0 == vid3)
		return IncidenceGraph::Vertex(vid0);

	if (vid1 == vid2 || vid1 == vid3)
		return IncidenceGraph::Vertex(vid1);

	return IncidenceGraph::Vertex();
}

inline void sorted_face_vertices(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge> face_edges,
								 std::vector<IncidenceGraph::Vertex> sorted_vertices)
{
	for (uint32 i = 0; i < face_edges.size(); ++i)
	{
		sorted_vertices.push_back(common_vertex(ig, face_edges[i], face_edges[(i + 1) % face_edges.size()]));
	}

	sorted_vertices.insert(sorted_vertices.begin(), sorted_vertices.back());
	sorted_vertices.pop_back();
	return;
}

inline std::pair<uint32, uint32> pseudoDegree(const IncidenceGraph& ig, IncidenceGraph::Vertex v)
{
	std::pair<uint32, uint32> info;
	info.first = 0; // incident 2D connex elements 
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
				++info.first;
				break;
			case 1: 
				info.first += 0.5;
				++info.second;
				break;
			case 2:
				break;
			default:
				info.first = INVALID_INDEX;
				break;
		}

		return (info.first != INVALID_INDEX);
	});

	return info;
}

} // namespace cgogn

#endif // CGOGN_CORE_INCIDENCE_GRAPH_OPS_H_
