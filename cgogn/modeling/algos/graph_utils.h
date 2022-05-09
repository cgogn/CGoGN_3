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

#ifndef CGOGN_MODELING_ALGOS_GRAPH_UTILS_H_
#define CGOGN_MODELING_ALGOS_GRAPH_UTILS_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/cell_marker.h>

namespace cgogn
{

namespace modeling
{

///////////
// Graph //
///////////

struct GraphData
{
	std::vector<std::pair<Graph::HalfEdge, Graph::HalfEdge>> branches;
	std::vector<Graph::Vertex> intersections;
};

inline Graph::HalfEdge branch_extremity(const Graph& g, Graph::HalfEdge h, CellMarker<Graph, Graph::Edge>& cm)
{
	Dart d = h.dart;
	while (degree(g, Graph::Vertex(d)) == 2)
	{
		d = alpha0(g, alpha1(g, d));
		cm.mark(Graph::Edge(d));
	}
	return Graph::HalfEdge(d);
}

bool get_graph_data(const Graph& g, GraphData& graph_data);

////////////////////
// IncidenceGraph //
////////////////////

using Branch = std::pair<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>,
						 std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>>;

struct IncidenceGraphData
{
	std::vector<Branch> branches;
	std::vector<IncidenceGraph::Vertex> extremities;
	std::vector<IncidenceGraph::Vertex> eejunctures;
	std::vector<IncidenceGraph::Vertex> efjunctures;
	std::vector<IncidenceGraph::Vertex> ffjunctures;
	std::vector<IncidenceGraph::Vertex> intersections;
	std::vector<std::vector<IncidenceGraph::Face>> leaflets;
};

std::pair<uint32, uint32> pseudo_degree(const IncidenceGraph& ig, IncidenceGraph::Vertex v);

uint32 get_incident_edge_id(const IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge e);

std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge> branch_extremity(
	const IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Vertex v,
	CellMarker<IncidenceGraph, IncidenceGraph::Edge>& marker);

bool get_incidenceGraph_data(const IncidenceGraph& ig, IncidenceGraphData& incidenceGraph_data);

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_GRAPH_UTILS_H_
