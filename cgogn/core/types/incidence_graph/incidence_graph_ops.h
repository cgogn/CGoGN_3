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

namespace cgogn
{

template <typename CELL>
CELL add_cell(IncidenceGraph& ig)
{
	uint32 id = ig.attribute_containers_[CELL::CELL_INDEX].new_index();
	return CELL(id);
}

template <typename CELL>
void remove_cell(IncidenceGraph& ig, CELL c)
{
	ig.attribute_containers_[CELL::CELL_INDEX].release_index(c.index_);
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

} // namespace cgogn

#endif // CGOGN_CORE_INCIDENCE_GRAPH_OPS_H_
