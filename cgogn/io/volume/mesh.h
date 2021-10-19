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

#ifndef CGOGN_IO_VOLUME_MESH_H_
#define CGOGN_IO_VOLUME_MESH_H_

#include <cgogn/io/utils.h>
#include <cgogn/io/volume/volume_import.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <fstream>
#include <vector>

namespace cgogn
{

namespace io
{

// template <typename MESH>
// bool import_MESH(MESH& m, const std::string& filename)
// {
// 	static_assert(mesh_traits<MESH>::dimension == 3, "MESH dimension should be 3");

// 	using Vertex = typename MESH::Vertex;

// 	Scoped_C_Locale loc;

// 	VolumeImportData volume_data;

// 	std::ifstream fp(filename, std::ios::in);

// 	std::string line;
// 	line.reserve(512u);

// 	getline_safe(fp, line);			  // MeshVersionFormatted 1
// 	getline_safe(fp, line);			  // Dimension
// 	uint32 dim = read_uint(fp, line); // 3

// 	getline_safe(fp, line); // Vertices
// 	uint32 nb_vertices = read_uint(fp, line);

// 	auto position = add_attribute<geometry::Vec3, Vertex>(m, "position");

// 	// read vertices position
// 	for (uint32 i = 0u; i < nb_vertices; ++i)
// 	{
// 		float64 x = read_double(fp, line);
// 		float64 y = read_double(fp, line);
// 		float64 z = read_double(fp, line);
// 		uint32 unused = read_uint(fp, line); // trailing 0

// 		uint32 vertex_id = new_index<Vertex>(m);
// 		(*position)[vertex_id] = {x, y, z};

// 		volume_data.vertices_id_.push_back(vertex_id);
// 	}

// 	getline_safe(fp, line); // Faces
// 	uint32 nb_faces = read_uint(fp, line);

// 	// read faces
// 	for (uint32 i = 0u; i < nb_faces; ++i)
// 	{
// 		getline_safe(fp, line); // unused for now
// 	}

// 	getline_safe(fp, line); // Hexahedra
// 	uint32 nb_volumes = read_uint(fp, line);

// 	// read volumes
// 	for (uint32 i = 0u; i < nb_volumes; ++i)
// 	{
// 		std::vector<uint32> ids(8);
// 		for (uint32 j = 0u; j < 8; ++j)
// 			ids[j] = volume_data.vertices_id_[read_uint(fp, line)];
// 		uint32 unused = read_uint(fp, line); // trailing 0

// 		std::vector<uint32> ordered_ids{ids[0], ids[4], ids[3], ids[2], ids[5], ids[8], ids[7], ids[6]};

// 		volume_data.volumes_types_.push_back(VolumeType::Hexa);
// 		volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ordered_ids.begin(),
// 												   ordered_ids.end());
// 	}

// 	import_volume_data(m, volume_data);

// 	return true;
// }

template <typename MESH>
void export_MESH(MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
				 const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 3, "MESH dimension should be 3");

	using Vertex = typename MESH::Vertex;
	using Face = typename MESH::Face;
	using Volume = typename MESH::Volume;

	auto vertex_id = add_attribute<uint32, Vertex>(m, "__vertex_id");

	uint32 nb_vertices = nb_cells<Vertex>(m);
	uint32 nb_faces = nb_cells<Face>(m);
	uint32 nb_volumes = nb_cells<Volume>(m);

	std::ofstream out_file;
	out_file.open(filename);
	out_file << "MeshVersionFormatted 1\n";
	out_file << "Dimension\n";
	out_file << "3\n";
	out_file << "Vertices\n";
	out_file << nb_vertices << "\n";

	uint32 id = 1;
	foreach_cell(m, [&](Vertex v) -> bool {
		const geometry::Vec3& p = value<geometry::Vec3>(m, vertex_position, v);
		value<uint32>(m, vertex_id, v) = id++;
		out_file << p[0] << " " << p[1] << " " << p[2] << " 0\n";
		return true;
	});

	out_file << "Quads\n";
	out_file << nb_faces << "\n";

	foreach_cell(m, [&](Face f) -> bool {
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			out_file << value<uint32>(m, vertex_id, v) << " ";
			return true;
		});
		out_file << "0\n";
		return true;
	});

	out_file << "Hexahedra\n";
	out_file << nb_volumes << "\n";

	foreach_cell(m, [&](Volume v) -> bool {
		Dart d1 = v.dart;
		Dart d2 = phi<2, 1, 1, 2>(m, d1);

		std::vector<uint32> indices{value<uint32>(m, vertex_id, Vertex(d1)),
									value<uint32>(m, vertex_id, Vertex(phi1(m, d2))),
									value<uint32>(m, vertex_id, Vertex(phi<1, 1>(m, d2))),
									value<uint32>(m, vertex_id, Vertex(phi_1(m, d1))),
									value<uint32>(m, vertex_id, Vertex(phi1(m, d1))),
									value<uint32>(m, vertex_id, Vertex(d2)),
									value<uint32>(m, vertex_id, Vertex(phi_1(m, d2))),
									value<uint32>(m, vertex_id, Vertex(phi<1, 1>(m, d1)))

		};

		for (uint32 i : indices)
			out_file << i << " ";

		out_file << "0\n";
		return true;
	});

	out_file << "End\n";

	remove_attribute<Vertex>(m, vertex_id);

	out_file.close();
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_VOLUME_MESH_H_
