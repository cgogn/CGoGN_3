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

#include <cgogn/io/incidence_graph/incidence_graph_import.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <vector>

namespace cgogn
{

namespace io
{

void import_incidence_graph_data(IncidenceGraph& ig, IncidenceGraphImportData& incidence_graph_data)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	auto position =
		get_or_add_attribute<geometry::Vec3, Vertex>(ig, incidence_graph_data.vertex_position_attribute_name_);

	std::vector<Vertex> vertices;
	vertices.reserve(incidence_graph_data.nb_vertices_);
	for (uint32 i = 0u; i < incidence_graph_data.nb_vertices_; ++i)
	{
		Vertex v = add_vertex(ig);
		(*position)[v.index_] = incidence_graph_data.vertex_position_[i];
		vertices.push_back(v);
	}

	std::vector<Edge> edges;
	edges.reserve(incidence_graph_data.nb_edges_);
	for (uint32 i = 0; i < incidence_graph_data.nb_edges_; ++i)
	{
		Edge e = add_edge(ig, vertices[incidence_graph_data.edges_vertex_indices_[2 * i]],
						  vertices[incidence_graph_data.edges_vertex_indices_[2 * i + 1]]);
		edges.push_back(e);
	}

	uint32 faces_edge_index = 0u;
	for (uint32 i = 0; i < incidence_graph_data.nb_faces_; ++i)
	{
		uint32 nbe = incidence_graph_data.faces_nb_edges_[i];
		std::vector<Edge> face;
		face.reserve(nbe);
		for (uint32 j = 0; j < nbe; ++j)
			face.push_back(edges[incidence_graph_data.faces_edge_indices_[faces_edge_index++]]);
		Face f = add_face(ig, face);
	}
}

void import_incidence_graph_data(Graph& g, IncidenceGraphImportData& incidence_graph_data)
{
	using Vertex = Graph::Vertex;

	auto position =
		get_or_add_attribute<geometry::Vec3, Vertex>(g, incidence_graph_data.vertex_position_attribute_name_);

	std::vector<Vertex> vertices;
	vertices.reserve(incidence_graph_data.nb_vertices_);
	for (uint32 i = 0u; i < incidence_graph_data.nb_vertices_; ++i)
	{
		Vertex v = add_vertex(g);
		value<Vec3>(g, position, v) = incidence_graph_data.vertex_position_[i];
		vertices.push_back(v);
	}

	for (uint32 i = 0; i < incidence_graph_data.nb_edges_; ++i)
	{
		connect_vertices(g, vertices[incidence_graph_data.edges_vertex_indices_[2 * i]],
						 vertices[incidence_graph_data.edges_vertex_indices_[2 * i + 1]]);
	}
}

} // namespace io

} // namespace cgogn
