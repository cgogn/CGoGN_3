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

#ifndef CGOGN_IO_SURFACE_OFF_H_
#define CGOGN_IO_SURFACE_OFF_H_

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <vector>
#include <fstream>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_OFF(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	SurfaceImportData surface_data;

	std::ifstream fp(filename.c_str(), std::ios::in);

	std::string line;
	line.reserve(512u);

	// read OFF header
	getline_safe(fp, line);
	if (line.rfind("OFF") == std::string::npos)
	{
		std::cerr << "File \"" << filename << "\" is not a valid off file." << std::endl;
		return false;
	}

	// read number of vertices, edges, faces
	const uint32 nb_vertices = read_uint(fp, line);
	const uint32 nb_faces = read_uint(fp, line);
	/*const uint32 nb_edges_ =*/ read_uint(fp, line);

	if (nb_vertices == 0u)
	{
		std::cerr << "File \"" << filename << " has no vertices." << std::endl;
		return false;
	}
	
	surface_data.reserve(nb_vertices, nb_faces);
	auto position = add_attribute<geometry::Vec3, CMap2::Vertex>(m, "position");

	// read vertices position
	for (uint32 i = 0u; i < nb_vertices; ++i)
	{
		float64 x = read_double(fp, line);
		float64 y = read_double(fp, line);
		float64 z = read_double(fp, line);

		uint32 vertex_id = new_cell_index<Vertex>(m);
		(*position)[vertex_id] = { x, y, z };

		surface_data.vertices_id_.push_back(vertex_id);
	}

	// read faces (vertex indices)
	for (uint32 i = 0u; i < nb_faces ; ++i)
	{
		uint32 n = read_uint(fp, line);
		
		std::vector<uint32> indices(n);
		for (uint32 j = 0u; j < n; ++j)
			indices[j] = surface_data.vertices_id_[read_uint(fp, line)];

		surface_data.faces_nb_vertices_.push_back(n);
		surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(), indices.end());
	}

	import_surface_data(m, surface_data);

	return true;
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_OFF_H_
