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

#ifndef CGOGN_IO_SURFACE_OBJ_H_
#define CGOGN_IO_SURFACE_OBJ_H_

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/utils.h>

#include <fstream>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_OBJ(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	SurfaceImportData surface_data;

	std::ifstream fp(filename.c_str(), std::ios::in);

	std::string line, tag;
	line.reserve(512u);

	// read vertices position
	do
	{
		fp >> tag;
		if (tag == std::string("v"))
		{
			surface_data.nb_vertices_++;
			float64 x = read_double(fp, line);
			float64 y = read_double(fp, line);
			float64 z = read_double(fp, line);
			surface_data.vertex_position_.push_back({x, y, z});
		}
		getline_safe(fp, line); // flush line
	} while (!fp.eof());

	if (surface_data.nb_vertices_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no vertices." << std::endl;
		return false;
	}

	// rewind
	fp.clear();
	fp.seekg(0, std::ios::beg);

	// read faces (vertex indices)
	// uint32 nb_faces = 0;
	do
	{
		fp >> tag;
		getline_safe(fp, line);
		if (tag == std::string("f"))
		{
			surface_data.nb_faces_++;
			std::vector<uint32> indices;
			std::istringstream iss(line);
			std::string str;
			while (!iss.eof())
			{
				iss >> str;
				uint32 ind = 0;
				while ((ind < str.length()) && (str[ind] != '/'))
					ind++;
				if (ind > 0)
				{
					uint32 index;
					std::stringstream iss(str.substr(0, ind));
					iss >> index;
					indices.push_back(index - 1);
				}
			}
			surface_data.faces_nb_vertices_.push_back(indices.size());
			surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
													  indices.end());
		}
	} while (!fp.eof());

	if (surface_data.nb_faces_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no faces." << std::endl;
		return false;
	}

	import_surface_data(m, surface_data);

	return true;
}

template <typename MESH>
void export_OBJ(MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
				const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	// TODO
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_OBJ_H_
