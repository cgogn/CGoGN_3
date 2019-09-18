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

#include <cgogn/io/volume/tet.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/orientation.h>

#include <vector>
#include <fstream>

namespace cgogn
{

namespace io
{

enum VolumeType
{
	Tetra,
	Pyramid,
	TriangularPrism,
	Hexa,
	Connector
};

bool import_TET(CMap3& m, const std::string& filename)
{
	Scoped_C_Locale loc;

	std::vector<VolumeType> volumes_types;
	std::vector<uint32> volumes_vertex_indices;
	
	std::ifstream fp(filename, std::ios::in);

	std::string line;
	line.reserve(512u);

	// read number of vertices
	uint32 nb_vertices = read_uint(fp, line);
	getline_safe(fp, line);

	uint32 nb_volumes = read_uint(fp, line);
	getline_safe(fp, line);

	std::cout << "nb_vertices: " << nb_vertices << std::endl;
	std::cout << "nb_volumes: " << nb_volumes << std::endl;

	volumes_types.reserve(nb_volumes);
	volumes_vertex_indices.reserve(8u * nb_volumes);

	auto position = add_attribute<geometry::Vec3, CMap3::Vertex>(m, "position");

	// read vertices position
	std::vector<uint32> vertices_id;
	vertices_id.reserve(nb_vertices);

	for (uint32 i = 0u; i < nb_vertices; ++i)
	{
		float64 x = read_double(fp, line);
		float64 y = read_double(fp, line);
		float64 z = read_double(fp, line);

		geometry::Vec3 pos{x, y, z};

		uint32 vertex_id = m.attribute_containers_[CMap3::Vertex::ORBIT].new_index();
		(*position)[vertex_id] = pos;

		vertices_id.push_back(vertex_id);
	}

	// read volumes
	for (uint32 i = 0u; i < nb_volumes; ++i)
	{
		uint32 n = read_uint(fp, line);
		std::vector<uint32> ids(n);
		for (uint32 j = 0u; j < n; ++j)
			ids[j] = vertices_id[read_uint(fp, line)];

		switch (n)
		{
			case 4: {
				if (geometry::test_orientation_3D((*position)[ids[0]], (*position)[ids[1]], (*position)[ids[2]], (*position)[ids[3]]) == geometry::Orientation3D::OVER)
					std::swap(ids[1], ids[2]);
				volumes_types.push_back(VolumeType::Tetra);
				volumes_vertex_indices.insert(volumes_vertex_indices.end(), ids.begin(), ids.end());
				break;
			}
			case 5: {
				if (geometry::test_orientation_3D((*position)[ids[4]], (*position)[ids[0]], (*position)[ids[1]], (*position)[ids[2]]) == geometry::Orientation3D::OVER)
					std::swap(ids[1], ids[3]);
				volumes_types.push_back(VolumeType::Pyramid);
				volumes_vertex_indices.insert(volumes_vertex_indices.end(), ids.begin(), ids.end());
				break;
			}
			case 6: {
				if (geometry::test_orientation_3D((*position)[ids[3]], (*position)[ids[0]], (*position)[ids[1]], (*position)[ids[2]]) == geometry::Orientation3D::OVER)
				{
					std::swap(ids[1], ids[2]);
					std::swap(ids[4], ids[5]);
				}
				volumes_types.push_back(VolumeType::TriangularPrism);
				volumes_vertex_indices.insert(volumes_vertex_indices.end(), ids.begin(), ids.end());
				break;
			}
			case 8: {
				if (geometry::test_orientation_3D((*position)[ids[4]], (*position)[ids[0]], (*position)[ids[1]], (*position)[ids[2]]) == geometry::Orientation3D::OVER)
				{
					std::swap(ids[0], ids[3]);
					std::swap(ids[1], ids[2]);
					std::swap(ids[4], ids[7]);
					std::swap(ids[5], ids[6]);
				}
				volumes_types.push_back(VolumeType::Hexa);
				volumes_vertex_indices.insert(volumes_vertex_indices.end(), ids.begin(), ids.end());
				break;
			}
			default:
				std::cout << "import_TET: Elements with " << n << " vertices are not handled. Ignoring.";
				break;
		}
	}

	if (volumes_vertex_indices.size() == 0)
		return false;
	
	auto darts_per_vertex = add_attribute<std::vector<Dart>, CMap3::Vertex>(m, "__darts_per_vertex");
	
	uint32 index = 0u;
	DartMarker dart_marker(m);
	uint32 vol_emb = 0u;

	// for each volume of table
	for (uint32 i = 0u, end = nb_volumes; i < end; ++i)
	{
		CMap3::Volume vol;
		const VolumeType vol_type = volumes_types[i];

		if (vol_type == VolumeType::Tetra) // tetrahedral case
		{
			vol = add_pyramid(static_cast<CMap2&>(m), 3u, false);

			const std::array<Dart, 4> vertices_of_tetra = {
				vol.dart,
				m.phi1(vol.dart),
				m.phi_1(vol.dart),
				m.phi_1(m.phi2(m.phi_1(vol.dart)))
			};

			for (const Dart& dv : vertices_of_tetra)
			{
				const uint32 vertex_index = volumes_vertex_indices[index++];
				static_cast<CMap2&>(m).foreach_dart_of_orbit(CMap2::Vertex(dv), [&] (Dart d) -> bool {
					m.set_embedding<CMap3::Vertex>(d, vertex_index);
					return true;
				});

				Dart dd = dv;
				do
				{
					dart_marker.mark(dd);
					(*darts_per_vertex)[vertex_index].push_back(dd);
					dd = m.phi1(m.phi2(dd));
				} while (dd != dv);
			}
		}
		else if (vol_type == VolumeType::Pyramid) // pyramidal case
		{
			vol = add_pyramid(static_cast<CMap2&>(m), 4u, false);

			const std::array<Dart, 5> vertices_of_pyramid = {
				vol.dart,
				m.phi1(vol.dart),
				m.phi1(m.phi1(vol.dart)),
				m.phi_1(vol.dart),
				m.phi_1(m.phi2(m.phi_1(vol.dart)))
			};

			for (Dart dv : vertices_of_pyramid)
			{
				const uint32 vertex_index = volumes_vertex_indices[index++];
				static_cast<CMap2&>(m).foreach_dart_of_orbit(CMap2::Vertex(dv), [&] (Dart d) -> bool {
					m.set_embedding<CMap3::Vertex>(d, vertex_index);
					return true;
				});

				Dart dd = dv;
				do
				{
					dart_marker.mark(dd);
					(*darts_per_vertex)[vertex_index].push_back(dd);
					dd = m.phi1(m.phi2(dd));
				} while (dd != dv);
			}
		}
		// else if (vol_type == VolumeType::TriangularPrism) // prism case
		// {
		// 	d = mbuild_.add_prism_topo_fp(3u);

		// 	// check if add ok (special maps)
		// 	if (d.is_nil()) break;

		// 	const std::array<Dart, 6> vertices_of_prism = {
		// 		d,
		// 		m.phi1(d),
		// 		m.phi_1(d),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(m.phi_1(d))))),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(d)))),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(m.phi1(d)))))
		// 	};

		// 	for (Dart dv : vertices_of_prism)
		// 	{
		// 		const uint32 vertex_index = this->volumes_vertex_indices_[index++];
		// 		mbuild_.template set_orbit_embedding<Vertex>(Vertex2(dv), vertex_index);

		// 		Dart dd = dv;
		// 		do
		// 		{
		// 			dart_marker.mark(dd);
		// 			(*darts_per_vertex)[vertex_index].push_back(dd);
		// 			dd = m.phi1(m.phi2(dd));
		// 		} while (dd != dv);
		// 	}
		// }
		// else if (vol_type == VolumeType::Hexa) // hexahedral case
		// {
		// 	d = mbuild_.add_prism_topo_fp(4u);

		// 	// check if add ok (special maps)
		// 	if (d.is_nil()) break;

		// 	const std::array<Dart, 8> vertices_of_hexa = {
		// 		d,
		// 		m.phi1(d),
		// 		m.phi1(m.phi1(d)),
		// 		m.phi_1(d),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(m.phi_1(d))))),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(d)))),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(m.phi1(d))))),
		// 		m.phi2(m.phi1(m.phi1(m.phi2(m.phi1(m.phi1(d))))))
		// 	};

		// 	for (Dart dv : vertices_of_hexa)
		// 	{
		// 		const uint32 vertex_index = this->volumes_vertex_indices_[index++];
		// 		mbuild_.template set_orbit_embedding<Vertex>(Vertex2(dv), vertex_index);

		// 		Dart dd = dv;
		// 		do
		// 		{
		// 			dart_marker.mark(dd);
		// 			(*darts_per_vertex)[vertex_index].push_back(dd);
		// 			dd = m.phi1(m.phi2(dd));
		// 		} while (dd != dv);
		// 	}
		// }
		// else //end of hexa
		// {
		// 	if (vol_type == VolumeType::Connector)
		// 	{
		// 		index += 4u;
		// 		// The second part of the code generates connectors automatically. We don't have to do anything here.
		// 	}
		// }

		if (m.is_embedded<CMap3::Volume>())
			set_embedding(m, vol, vol_emb++);
	}

	// reconstruct neighbourhood
	uint32 nb_boundary_faces = 0u;
	DartMarkerStore marker(m);

	m.foreach_dart([&] (Dart d) -> bool
	{
		if (m.phi3(d) == d && !marker.is_marked(d))
		{
			static_cast<CMap2&>(m).foreach_dart_of_orbit(CMap2::Face(d), [&] (Dart fd) -> bool { marker.mark(fd); return true; });

			Dart good_dart;

			// 1st step : for every dart of the face we try to find a valid phi3 candidate. If we can't it's a boundary face.
			Dart d_it = d;
			do
			{
				uint32 vindex1 = m.embedding(CMap3::Vertex(d_it));
				uint32 vindex2 = m.embedding(CMap3::Vertex(m.phi1(m.phi1(d_it))));
				const std::vector<Dart>& vec = value<std::vector<Dart>>(m, darts_per_vertex, CMap3::Vertex(m.phi1(d_it)));
				for (auto it = vec.begin(); it != vec.end() && good_dart.is_nil(); ++it)
					if (m.embedding(CMap3::Vertex(m.phi1(*it))) == vindex1 && m.embedding(CMap3::Vertex(m.phi_1(*it))) == vindex2)
						good_dart = *it;
				d_it = m.phi1(d_it);
			} while (good_dart.is_nil() && (d_it != d));
			d = m.phi_1(d_it);

			if (!good_dart.is_nil())
			{
				const uint32 degD = codegree(m, CMap3::Face(d));
				const uint32 degGD = codegree(m, CMap3::Face(good_dart));

				if (degD == degGD) // normal case : the two opposite faces have the same degree
				{
					Dart it1 = d;
					Dart it2 = good_dart;
					do
					{
						m.phi3_sew(it1, it2);
						it1 = m.phi1(it1);
						it2 = m.phi_1(it2);
					} while (it1 != d);
				}
				else
				{
					// there is one face of degree 4 and one face of degree 3
					// -> stamp volume
				}
			}
			else
				++nb_boundary_faces;
		}

		return true;
	});

	if (nb_boundary_faces > 0u)
	{
		uint32 nb_holes = m.close();
		std::cout << nb_holes << " hole(s) have been closed" << std::endl;
		std::cout << nb_boundary_faces << " boundary faces" << std::endl;
	}
	
	remove_attribute<CMap3::Vertex>(m, darts_per_vertex);

	m.foreach_dart([&] (Dart d) -> bool
	{
		if (m.phi3(d) == d)
			std::cout << "phi3 fix point Dart: " << d << std::endl;
		return true;
	});

	return true;
}

} // namespace io

} // namespace cgogn
