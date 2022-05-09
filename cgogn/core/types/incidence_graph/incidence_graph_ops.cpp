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

#include <cgogn/core/types/incidence_graph/incidence_graph.h>

#include <iostream>

namespace cgogn
{

bool sort_face_edges(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<uint8>& edges_dir = (*ig.face_incident_edges_dir_)[f.index_];
	edges_dir.clear();

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
		std::size_t i, end;
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

} // namespace cgogn
