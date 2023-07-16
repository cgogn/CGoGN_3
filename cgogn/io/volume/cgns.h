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

#ifndef CGOGN_IO_VOLUME_CGNS_H_
#define CGOGN_IO_VOLUME_CGNS_H_

namespace cgns
{
extern "C"
{
#include <thirdparty/cgns-3.4.2/src/cgnslib.h>
}
} // namespace cgns

#include <cgogn/io/utils.h>
#include <cgogn/io/volume/volume_import.h>

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <fstream>
#include <vector>

namespace cgogn
{

namespace io
{

template <typename MESH>
void export_CGNS(MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
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

	int index_file;
	cgns::cg_open(filename.c_str(), CG_MODE_WRITE, &index_file);

	int index_base;
	cgns::cg_base_write(index_file, "Base", 3, 3, &index_base);

	cgns::cg_goto(index_file, index_base, NULL);
	cgns::cg_dataclass_write(cgns::NormalizedByDimensional);

	cgns::cgsize_t isize[3];
	isize[0] = nb_vertices;
	isize[1] = nb_volumes;
	isize[2] = 0;
	int index_zone;
	cgns::cg_zone_write(index_file, index_base, "Zone1", isize, cgns::Unstructured, &index_zone);

	uint32 vert_id = 0;
	std::vector<geometry::Scalar> x(nb_vertices);
	std::vector<geometry::Scalar> y(nb_vertices);
	std::vector<geometry::Scalar> z(nb_vertices);
	foreach_cell(m, [&](Vertex v) -> bool {
		const geometry::Vec3& p = value<geometry::Vec3>(m, vertex_position, v);
		x[vert_id] = p[0];
		y[vert_id] = p[1];
		z[vert_id] = p[2];
		value<uint32>(m, vertex_id, v) = ++vert_id;
		return true;
	});

	int index_coord;
	cgns::cg_coord_write(index_file, index_base, index_zone, cgns::RealDouble, "CoordinateX", x.data(), &index_coord);
	// cgns::cg_goto(index_file, index_base, "Zone_t", index_zone, "GridCoordinates", 0, "CoordinateX", 0, NULL);
	// cgns::cg_exponents_write(cgns::RealSingle, exp);
	cgns::cg_coord_write(index_file, index_base, index_zone, cgns::RealDouble, "CoordinateY", y.data(), &index_coord);
	cgns::cg_coord_write(index_file, index_base, index_zone, cgns::RealDouble, "CoordinateZ", z.data(), &index_coord);

	// foreach_cell(m, [&](Face f) -> bool {
	// 	foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
	// 		out_file << value<uint32>(m, vertex_id, v) << " ";
	// 		return true;
	// 	});
	// 	out_file << "0\n";
	// 	return true;
	// });

	uint32 vol_id = 0;
	uint32 nb_indices = 0;
	std::vector<cgns::cgsize_t> volumes_vertex_indices(nb_volumes * 8);
	// cgns::cgsize_t volumes_vertex_indices[nb_volumes][8];
	foreach_cell(m, [&](Volume v) -> bool {
		Dart d1 = v.dart;
		Dart d2 = phi<2, 1, 1, 2>(m, d1);
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(d1));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(phi1(m, d2)));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(phi<1, 1>(m, d2)));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(phi_1(m, d1)));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(phi1(m, d1)));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(d2));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(phi_1(m, d2)));
		volumes_vertex_indices[nb_indices++] = value<uint32>(m, vertex_id, Vertex(phi<1, 1>(m, d1)));
		vol_id++;
		return true;
	});

	int index_section;
	cgns::cg_section_write(index_file, index_base, index_zone, "Elem", cgns::HEXA_8, 1, nb_volumes, 0,
						   volumes_vertex_indices.data(), &index_section);

	cgns::cg_close(index_file);

	remove_attribute<Vertex>(m, vertex_id);
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_VOLUME_CGNS_H_
