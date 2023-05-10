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

#ifndef CGOGN_IO_INCIDENCE_GRAPH_IMPORT_H_
#define CGOGN_IO_INCIDENCE_GRAPH_IMPORT_H_

#include <cgogn/io/cgogn_io_export.h>

#include <cgogn/core/types/maps/cmap/graph.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

namespace cgogn
{

namespace io
{

using geometry::Vec3;

struct IncidenceGraphImportData
{
	uint32 nb_vertices_ = 0;
	uint32 nb_edges_ = 0;
	uint32 nb_faces_ = 0;

	std::vector<Vec3> vertex_position_;
	std::string vertex_position_attribute_name_ = "position";

	std::vector<uint32> edges_vertex_indices_;
	std::vector<uint32> faces_nb_edges_;
	std::vector<uint32> faces_edge_indices_;

	inline void reserve(uint32 nb_vertices, uint32 nb_edges, uint32 nb_faces)
	{
		nb_vertices_ = nb_vertices;
		nb_edges_ = nb_edges;
		nb_faces_ = nb_faces;
		vertex_position_.reserve(nb_vertices);
		edges_vertex_indices_.reserve(nb_edges * 2u);
		faces_nb_edges_.reserve(nb_faces);
		faces_edge_indices_.reserve(nb_faces * 4u);
	}
};

void CGOGN_IO_EXPORT import_incidence_graph_data(IncidenceGraph& ig, IncidenceGraphImportData& graph_data);

void CGOGN_IO_EXPORT import_incidence_graph_data(Graph& g, IncidenceGraphImportData& graph_data);

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_INCIDENCE_GRAPH_IMPORT_H_
