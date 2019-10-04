/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_IO_IMPORT_GRAPH_H_
#define CGOGN_IO_IMPORT_GRAPH_H_

#include <cgogn/io/cgogn_io_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <vector>

namespace cgogn
{

namespace io
{

struct GraphImportData
{
	std::vector<uint32> vertices_id_;
	std::vector<uint32> edges_vertex_indices_;
};

void import_graph_data(Graph& g, const GraphImportData& graph_data)
{
	auto vertex_dart = add_attribute<Dart, Graph::Vertex>(g, "__vertex_dart");

	for (uint32 i : graph_data.vertices_id_)
	{
		Graph::Vertex v = add_vertex(g, false);
		g.set_embedding<Graph::Vertex>(v.dart, i);
		(*vertex_dart)[i] = v.dart;
	}

	for (uint32 i = 0; i < graph_data.edges_vertex_indices_.size(); i += 2)
		connect_vertices(
			g,
			Graph::Vertex((*vertex_dart)[graph_data.vertices_id_[graph_data.edges_vertex_indices_[i]]]),
			Graph::Vertex((*vertex_dart)[graph_data.vertices_id_[graph_data.edges_vertex_indices_[i+1]]])
		);

	remove_attribute<Graph::Vertex>(g, vertex_dart);
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_IMPORT_GRAPH_H_
