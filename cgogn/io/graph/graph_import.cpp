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

#include <cgogn/io/graph/graph_import.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/core/types/cmap/cmap_ops.h>

#include <vector>

namespace cgogn
{

namespace io
{

void import_graph_data(Graph& g, const GraphImportData& graph_data)
{
	using Vertex = Graph::Vertex;

	auto vertex_dart = add_attribute<Dart, Vertex>(g, "__vertex_dart");

	for (uint32 vertex_id : graph_data.vertices_id_)
	{
		Vertex v = add_vertex(g, false);
		set_index(g, v, vertex_id);
		(*vertex_dart)[vertex_id] = v.dart;
	}

	for (uint32 i = 0; i < uint32(graph_data.edges_vertex_indices_.size()); i += 2)
	{
		connect_vertices(g, Vertex((*vertex_dart)[graph_data.edges_vertex_indices_[i]]),
						 Vertex((*vertex_dart)[graph_data.edges_vertex_indices_[i + 1]]));
	}

	remove_attribute<Vertex>(g, vertex_dart);
}

void import_graph_data(IncidenceGraph& ig, const GraphImportData& graph_data)
{
	using Vertex = IncidenceGraph::Vertex;
	for (uint32 i = 0; i < uint32(graph_data.edges_vertex_indices_.size()); i += 2)
		add_edge(ig, Vertex(graph_data.edges_vertex_indices_[i]), Vertex(graph_data.edges_vertex_indices_[i + 1]));
}

} // namespace io

} // namespace cgogn
