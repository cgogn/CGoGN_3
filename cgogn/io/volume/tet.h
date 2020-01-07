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

#ifndef CGOGN_IO_VOLUME_TET_H_
#define CGOGN_IO_VOLUME_TET_H_

#include <cgogn/io/utils.h>
#include <cgogn/io/volume/volume_import.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/functions/orientation.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <fstream>
#include <vector>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_TET(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 3, "MESH dimension should be 3");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	VolumeImportData volume_data;

	std::ifstream fp(filename, std::ios::in);

	std::string line;
	line.reserve(512u);

	// read number of vertices
	uint32 nb_vertices = read_uint(fp, line);
	getline_safe(fp, line);

	uint32 nb_volumes = read_uint(fp, line);
	getline_safe(fp, line);

	if (nb_vertices == 0u)
	{
		std::cerr << "File \"" << filename << " has no vertices." << std::endl;
		return false;
	}

	volume_data.reserve(nb_vertices, nb_volumes);
	auto position = add_attribute<geometry::Vec3, Vertex>(m, "position");

	// read vertices position
	for (uint32 i = 0u; i < nb_vertices; ++i)
	{
		float64 x = read_double(fp, line);
		float64 y = read_double(fp, line);
		float64 z = read_double(fp, line);

		uint32 vertex_id = new_index<Vertex>(m);
		(*position)[vertex_id] = {x, y, z};

		volume_data.vertices_id_.push_back(vertex_id);
	}

	// read volumes
	for (uint32 i = 0u; i < nb_volumes; ++i)
	{
		uint32 n = read_uint(fp, line);
		std::vector<uint32> ids(n);
		for (uint32 j = 0u; j < n; ++j)
			ids[j] = volume_data.vertices_id_[read_uint(fp, line)];

		switch (n)
		{
		case 4: {
			if (geometry::test_orientation_3D((*position)[ids[0]], (*position)[ids[1]], (*position)[ids[2]],
											  (*position)[ids[3]]) == geometry::Orientation3D::UNDER)
				std::swap(ids[1], ids[2]);
			volume_data.volumes_types_.push_back(VolumeType::Tetra);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
			break;
		}
		case 5: {
			if (geometry::test_orientation_3D((*position)[ids[4]], (*position)[ids[0]], (*position)[ids[1]],
											  (*position)[ids[2]]) == geometry::Orientation3D::OVER)
				std::swap(ids[1], ids[3]);
			volume_data.volumes_types_.push_back(VolumeType::Pyramid);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
			break;
		}
		case 6: {
			if (geometry::test_orientation_3D((*position)[ids[3]], (*position)[ids[0]], (*position)[ids[1]],
											  (*position)[ids[2]]) == geometry::Orientation3D::OVER)
			{
				std::swap(ids[1], ids[2]);
				std::swap(ids[4], ids[5]);
			}
			volume_data.volumes_types_.push_back(VolumeType::TriangularPrism);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
			break;
		}
		case 8: {
			if (geometry::test_orientation_3D((*position)[ids[4]], (*position)[ids[0]], (*position)[ids[1]],
											  (*position)[ids[2]]) == geometry::Orientation3D::OVER)
			{
				std::swap(ids[0], ids[3]);
				std::swap(ids[1], ids[2]);
				std::swap(ids[4], ids[7]);
				std::swap(ids[5], ids[6]);
			}
			volume_data.volumes_types_.push_back(VolumeType::Hexa);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
			break;
		}
		default:
			std::cout << "import_TET: Elements with " << n << " vertices are not handled. Ignoring.";
			break;
		}
	}

	import_volume_data(m, volume_data);

	return true;
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_VOLUME_TET_H_
