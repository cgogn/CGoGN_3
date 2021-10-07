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

#ifndef CGOGN_IO_GRAPH_SKEL_H_
#define CGOGN_IO_GRAPH_SKEL_H_

#include <cgogn/io/graph/graph_import.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <fstream>
#include <set>
#include <vector>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_SKEL(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 1, "MESH dimension should be 1");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	GraphImportData graph_data;

	std::ifstream fp(filename.c_str(), std::ios::in);

	std::string line;
	line.reserve(512);
	getline_safe(fp, line); // Discard first line, it's useless
	getline_safe(fp, line); // Number of vertices

	std::stringstream issl(line);
	uint32 value;
	issl >> value;
	const uint32 nb_vertices = value;

	if (nb_vertices == 0u)
	{
		std::cerr << "File \"" << filename << " has no vertices." << std::endl;
		return false;
	}

	graph_data.reserve(nb_vertices);
	auto position = add_attribute<geometry::Vec3, Vertex>(m, "position");
	auto radius = add_attribute<geometry::Scalar, Vertex>(m, "radius");

	for (uint32 i = 0; i < nb_vertices; ++i)
	{
		uint32 id;
		float64 x, y, z, r;
		uint32 nb_neighbors;

		getline_safe(fp, line);
		std::stringstream iss(line);
		iss >> id;
		iss >> x >> y >> z;
		iss >> r;

		uint32 vertex_id = new_index<Vertex>(m);
		(*position)[vertex_id] = {x, y, z};
		(*radius)[vertex_id] = r;
		graph_data.vertices_id_.push_back(vertex_id);

		iss >> nb_neighbors;

		for (uint32 j = 0; j < nb_neighbors; ++j)
		{
			uint32 neighbor_id;
			iss >> neighbor_id;
			if (neighbor_id < id)
			{
				graph_data.edges_vertex_indices_.push_back(graph_data.vertices_id_[id]);
				graph_data.edges_vertex_indices_.push_back(graph_data.vertices_id_[neighbor_id]);
			}
		}
	}

	import_graph_data(m, graph_data);

	return true;
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_GRAPH_SKEL_H_
