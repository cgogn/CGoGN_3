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

#include <cgogn/io/graph/skel.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <vector>
#include <fstream>
#include <set>

#include <cgogn/core/types/cmap/cmap_info.h>

namespace cgogn
{

namespace io
{

bool import_SKEL(Graph& g, const std::string& filename)
{
	// GENERIC
	Scoped_C_Locale loc;
	
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

	std::vector<uint32> vertices_indices(nb_vertices);
	std::vector<geometry::Scalar> vertices_radii(nb_vertices);
	std::vector<geometry::Vec3> vertices_pos(nb_vertices);
	
	std::vector<uint32> edges_vertex_indices;
	edges_vertex_indices.reserve(nb_vertices * 3);

	for (uint32 i = 0; i < nb_vertices; ++i)
	{
		uint32 vertex_id;
		float64 x, y, z, radius;
		uint32 nb_neighbors;

		getline_safe(fp, line); 
		std::stringstream iss(line);
		iss >> vertex_id;
		iss >> x >> y >> z;
		vertices_pos[vertex_id] = {x, y, z};
		iss >> radius;
		vertices_radii[vertex_id] = radius;
		iss >> nb_neighbors;

		for (uint32 j = 0; j < nb_neighbors; ++j)
		{
			uint32 neighbor_id;
			iss >> neighbor_id;
			if (neighbor_id > vertex_id)
			{
				edges_vertex_indices.push_back(vertex_id);
				edges_vertex_indices.push_back(neighbor_id);
			}
		}
	}

	if (edges_vertex_indices.size() == 0u)
	{
		std::cerr << "File \"" << filename << " has no edges." << std::endl;
		return false;
	}

	// SPECIALIZED
	auto position = add_attribute<geometry::Vec3, Graph::Vertex>(g, "position");
	auto radius = add_attribute<geometry::Scalar, Graph::Vertex>(g, "radius");
	for (uint32 i = 0; i < nb_vertices; ++i)
	{
		uint32 vertex_id = g.attribute_containers_[Graph::Vertex::ORBIT].new_index();
		(*position)[vertex_id] = vertices_pos[i];
		(*radius)[vertex_id] = vertices_radii[i];
		vertices_indices[i] = vertex_id;
	}

	auto vertex_dart = add_attribute<Dart, Graph::Vertex>(g, "__vertex_dart");

	CMapBase::AttributeContainer& vc = g.attribute_containers_[Graph::Vertex::ORBIT];
	for (uint32 i = vc.first_index(); i != vc.last_index(); i = vc.next_index(i))
	{
		Graph::Vertex v = add_vertex(g, false);
		g.set_embedding<Graph::Vertex>(v.dart, i);
		(*vertex_dart)[i] = v.dart;
	}

	for (uint32 i = 0; i < edges_vertex_indices.size(); i += 2)
		connect_vertices(
			g,
			Graph::Vertex((*vertex_dart)[vertices_indices[edges_vertex_indices[i]]]),
			Graph::Vertex((*vertex_dart)[vertices_indices[edges_vertex_indices[i+1]]])
		);

	remove_attribute<Graph::Vertex>(g, vertex_dart);

	return true;
}

} // namespace io

} // namespace cgogn
