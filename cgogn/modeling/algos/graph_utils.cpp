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

#include <cgogn/modeling/algos/graph_utils.h>

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/halfedge.h>

namespace cgogn
{

namespace modeling
{

///////////
// Graph //
///////////

bool get_graph_data(const Graph& g, GraphData& graph_data)
{
	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		if (degree(g, v) > 2)
			graph_data.intersections.push_back(v);
		return true;
	});

	CellMarker<Graph, Graph::Edge> cm(g);
	foreach_cell(g, [&](Graph::Edge e) -> bool {
		if (cm.is_marked(e))
			return true;
		cm.mark(e);
		std::vector<Graph::HalfEdge> halfedges = incident_halfedges(g, e);
		graph_data.branches.push_back({branch_extremity(g, halfedges[0], cm), branch_extremity(g, halfedges[1], cm)});
		return true;
	});

	return true;
}

////////////////////
// IncidenceGraph //
////////////////////

std::pair<uint32, uint32> pseudo_degree(const IncidenceGraph& ig, IncidenceGraph::Vertex v)
{
	std::pair<uint32, uint32> info;
	info.first = 0;	 // incident leaflets
	info.second = 0; // incident branches
	// vertices inside leaflets are (0,0)
	// vertices on the boundary of leaflets are (1,0)

	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
		uint32 nbf = degree(ig, e);
		switch (nbf)
		{
		case 0:
			info.second++;
			break;
		case 1:
			info.first++;
			break;
		case 2:
			break;
		default:
			info.first = INVALID_INDEX;
			break;
		}

		return (info.first != INVALID_INDEX);
	});
	if (info.first != INVALID_INDEX)
		info.first /= 2;
	return info;
}

uint32 get_incident_edge_id(const IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge e)
{
	const std::vector<IncidenceGraph::Edge>& inc_edges = (*ig.face_incident_edges_)[f.index_];
	for (uint32 i = 0, end = inc_edges.size(); i < end; ++i)
	{
		if (inc_edges[i] == e)
			return i;
	}
	return INVALID_INDEX;
}

std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge> branch_extremity(
	const IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Vertex v,
	CellMarker<IncidenceGraph, IncidenceGraph::Edge>& marker)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	Edge current_edge = e;
	Vertex next_vertex = v;
	while (degree(ig, next_vertex) == 2)
	{
		const std::vector<Edge>& inc_edges = (*ig.vertex_incident_edges_)[next_vertex.index_];
		current_edge = inc_edges[0] == current_edge ? inc_edges[1] : inc_edges[0];
		const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[current_edge.index_];
		next_vertex = inc_verts.first == next_vertex ? inc_verts.second : inc_verts.first;
		marker.mark(current_edge);
	}
	return {next_vertex, current_edge};
}

bool get_incidenceGraph_data(const IncidenceGraph& ig, IncidenceGraphData& igData)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	bool success = true;

	// classify vertices
	foreach_cell(ig, [&](Vertex v) -> bool {
		std::pair<uint32, uint32> info = pseudo_degree(ig, v); // { nb_leaflets, nb_edges }

		// vertex is incident to a fan edge
		if (info.first == INVALID_INDEX)
		{
			// TODO: check that it is a supported configuration
			// in case of multiple fan edges, these should be connected..
			uint32 nb_fan_edges = 0;
			foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
				uint32 d = degree(ig, e);
				if (d > 1)
					nb_fan_edges++;
				return true;
			});
			if (nb_fan_edges == 1)
				igData.leaflets_boundary_vertices_fans.push_back(v);
			return true;
		}

		if (info.first + info.second > 2)
			igData.intersections.push_back(v);
		else if (info.first == 0 && info.second == 0)
			igData.leaflets_inside_vertices.push_back(v);
		else if (info.first == 0 && info.second == 1)
			igData.extremities.push_back(v);
		else if (info.first == 0 && info.second == 2)
			igData.eejunctures.push_back(v);
		else if (info.first == 1 && info.second == 0)
			igData.leaflets_boundary_vertices_corners.push_back(v);
		else if (info.first == 1 && info.second == 1)
			igData.efjunctures.push_back(v);
		else if (info.first == 2 && info.second == 0)
			igData.ffjunctures.push_back(v);

		return true;
	});

	if (!success)
		return success;

	// get branches, fan edges & leaflet boundary edges
	CellMarker<IncidenceGraph, Edge> edge_marker(ig);
	foreach_cell(ig, [&](Edge e) -> bool {
		uint32 d = degree(ig, e);
		if (d == 1)
			igData.leaflets_boundary_edges.push_back(e);
		if (d > 2)
			igData.fan_edges.push_back(e);

		if (d > 0)
			return true;

		if (edge_marker.is_marked(e))
			return true;
		edge_marker.mark(e);

		const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[e.index_];
		igData.branches.push_back({branch_extremity(ig, e, inc_verts.first, edge_marker),
								   branch_extremity(ig, e, inc_verts.second, edge_marker)});
		return true;
	});

	// get each leaflet faces
	// for each leaflet, the orientation of the first face is propagated to all the faces of the leaflet
	CellMarker<IncidenceGraph, Face> face_marker(ig);
	foreach_cell(ig, [&](Face fa) -> bool {
		if (face_marker.is_marked(fa))
			return true;

		face_marker.mark(fa);
		std::vector<Face>& leaflet_faces = igData.leaflets.emplace_back();
		leaflet_faces.push_back(fa);
		for (uint32 i = 0; i < leaflet_faces.size(); ++i)
		{
			Face f = leaflet_faces[i];
			const std::vector<Edge>& inc_edges = (*ig.face_incident_edges_)[f.index_];
			const std::vector<uint8>& inc_edges_dir = (*ig.face_incident_edges_dir_)[f.index_];
			for (uint32 j = 0; j < inc_edges.size(); ++j)
			{
				Edge e = inc_edges[j];
				uint8 edge_dir = inc_edges_dir[j];
				const std::vector<Face>& inc_faces = (*ig.edge_incident_faces_)[e.index_];
				if (inc_faces.size() == 2)
				{
					Face af = inc_faces[0] == f ? inc_faces[1] : inc_faces[0];
					if (!face_marker.is_marked(af))
					{
						uint32 eid = get_incident_edge_id(ig, af, e);
						std::vector<Edge>& af_inc_edges = (*ig.face_incident_edges_)[af.index_];
						std::vector<uint8>& af_inc_edges_dir = (*ig.face_incident_edges_dir_)[af.index_];
						if (af_inc_edges_dir[eid] == edge_dir)
						{
							for (uint8& dir : af_inc_edges_dir)
								dir = (dir + 1) % 2;
							std::reverse(af_inc_edges.begin(), af_inc_edges.end());
							std::reverse(af_inc_edges_dir.begin(), af_inc_edges_dir.end());
						}
						face_marker.mark(af);
						leaflet_faces.push_back(af);
					}
				}
			}
		}
		return true;
	});

	return success;
}

} // namespace modeling

} // namespace cgogn
