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

#include <cgogn/io/graph/cg.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <vector>
#include <fstream>


#include <cgogn/core/types/cmap/cmap_info.h>

namespace cgogn
{

namespace io
{

bool import_CG(Graph& g, const std::string& filename)
{
	Scoped_C_Locale loc;

	std::vector<uint32> edges_vertex_indices;

	std::ifstream fp(filename.c_str(), std::ios::in);

	std::string line;
	line.reserve(512);

	getline_safe(fp, line);
	if (line.rfind("# D") == std::string::npos)
	{
		std::cerr << "File \"" << filename << "\" is not a valid cg file." << std::endl;
		return false;
	}
	
	// read header
	std::replace(line.begin(), line.end(), ':', ' ');
	std::stringstream issl(line);
	std::string tagl;
	uint32 value;
	issl >> tagl;
	issl >> tagl;
	issl >> value; // dimension unused for now
	issl >> tagl;
	issl >> value;
	const uint32 nb_vertices = value;
	issl >> tagl;
	issl >> value;
	const uint32 nb_edges = value;

	edges_vertex_indices.reserve(nb_edges * 2);

	auto position = add_attribute<geometry::Vec3, Graph::Vertex>(g, "position");

	// read vertices
	std::vector<uint32> vertices_id;
	vertices_id.reserve(nb_vertices);

	for (uint32 i = 0; i < nb_vertices; ++i)
	{
		getline_safe(fp, line);
		std::stringstream iss(line);

		std::string tag;
		iss >> tag;

		if (tag == std::string("v"))
		{
			float64 x, y, z;
			iss >> x;
			iss >> y;
			iss >> z;

			uint32 vertex_id = g.attribute_containers_[Graph::Vertex::ORBIT].new_index();
			(*position)[vertex_id] = { x, y, z };

			vertices_id.push_back(vertex_id);
		}
	}

	// read edges
	for (uint32 i = 0; i < nb_edges; ++i)
	{
		getline_safe(fp, line);
		std::stringstream iss(line);

		std::string tag;
		iss >> tag;

		if (tag == std::string("e"))
		{
			uint32 a, b;
			iss >> a;
			iss >> b;

			edges_vertex_indices.push_back(vertices_id[a]);
			edges_vertex_indices.push_back(vertices_id[b]);
		}
	}

	if (edges_vertex_indices.size() == 0u)
	{
		std::cerr << "File \"" << filename << " has no edges." << std::endl;
		return false;
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
		connect_vertices(g, Graph::Vertex((*vertex_dart)[edges_vertex_indices[i]]), Graph::Vertex((*vertex_dart)[edges_vertex_indices[i+1]]));

	remove_attribute<Graph::Vertex>(g, vertex_dart);
    dump_map(g);
	return true;
}

} // namespace io

} // namespace cgogn
