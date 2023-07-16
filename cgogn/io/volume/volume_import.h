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

#ifndef CGOGN_IO_VOLUME_IMPORT_H_
#define CGOGN_IO_VOLUME_IMPORT_H_

#include <cgogn/io/cgogn_io_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

namespace cgogn
{

struct CMap3;
struct GMap3;

namespace io
{

using geometry::Vec3;

enum VolumeType
{
	Tetra,
	Pyramid,
	TriangularPrism,
	Hexa,
	Connector
};

struct VolumeImportData
{
	uint32 nb_vertices_ = 0;
	uint32 nb_volumes_ = 0;

	std::vector<Vec3> vertex_position_;
	std::string vertex_position_attribute_name_ = "position";

	std::vector<VolumeType> volumes_types_;
	std::vector<uint32> volumes_vertex_indices_;

	std::vector<uint32> vertex_id_after_import_;

	inline void reserve(uint32 nb_vertices, uint32 nb_volumes)
	{
		nb_vertices_ = nb_vertices;
		nb_volumes_ = nb_volumes;
		vertex_position_.reserve(nb_vertices);
		volumes_types_.reserve(nb_volumes);
		volumes_vertex_indices_.reserve(nb_volumes * 8u);
		vertex_id_after_import_.reserve(nb_vertices);
	}
};

template <typename MESH>
auto import_volume_data(MESH& m, VolumeImportData& volume_data)
	-> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	using Vertex = typename MESH::Vertex;
	using Edge = typename MESH::Edge;
	using Face = typename MESH::Face;
	using Volume = typename MESH::Volume;
	using Vertex2 = typename MESH::Vertex2;
	using Face2 = typename MESH::Face2;

	using ParentMESH = typename MESH::Parent;

	auto position = get_or_add_attribute<geometry::Vec3, Vertex>(m, volume_data.vertex_position_attribute_name_);

	for (uint32 i = 0u; i < volume_data.nb_vertices_; ++i)
	{
		uint32 vertex_id = new_index<Vertex>(m);
		(*position)[vertex_id] = volume_data.vertex_position_[i];
		volume_data.vertex_id_after_import_.push_back(vertex_id);
	}

	auto darts_per_vertex = add_attribute<std::vector<Dart>, Vertex>(m, "__darts_per_vertex");

	uint32 index = 0u;
	DartMarker dart_marker(m);
	uint32 vol_emb = 0u;

	auto mark_n_store = [&](Dart dv, const uint32 vertex_index) {
		foreach_dart_of_orbit(m, Vertex2(dv), [&](Dart dd) {
			dart_marker.mark(dd);
			(*darts_per_vertex)[vertex_index].push_back(dd);
			return true;
		});
	};

	// for each volume of table
	for (uint32 i = 0u; i < volume_data.nb_volumes_; ++i)
	{
		Volume vol;
		const VolumeType vol_type = volume_data.volumes_types_[i];

		if (vol_type == VolumeType::Tetra) // tetrahedral case
		{
			vol = add_pyramid(static_cast<ParentMESH&>(m), 3u, false);

			const std::array<Dart, 4> vertices_of_tetra = {vol.dart, phi1(m, vol.dart), phi_1(m, vol.dart),
														   phi<-1, 2, -1>(m, vol.dart)};

			for (Dart dv : vertices_of_tetra)
			{
				const uint32 vertex_index =
					volume_data.vertex_id_after_import_[volume_data.volumes_vertex_indices_[index++]];
				foreach_dart_of_orbit(m, Vertex2(dv), [&](Dart d) -> bool {
					set_index<Vertex>(m, d, vertex_index);
					return true;
				});
				mark_n_store(dv, vertex_index);
			}
		}
		else if (vol_type == VolumeType::Pyramid) // pyramidal case
		{
			vol = add_pyramid(static_cast<ParentMESH&>(m), 4u, false);

			const std::array<Dart, 5> vertices_of_pyramid = {vol.dart, phi1(m, vol.dart), phi<1, 1>(m, vol.dart),
															 phi_1(m, vol.dart), phi<-1, 2, -1>(m, vol.dart)};

			for (Dart dv : vertices_of_pyramid)
			{
				const uint32 vertex_index =
					volume_data.vertex_id_after_import_[volume_data.volumes_vertex_indices_[index++]];
				foreach_dart_of_orbit(m, Vertex2(dv), [&](Dart d) -> bool {
					set_index<Vertex>(m, d, vertex_index);
					return true;
				});
				mark_n_store(dv, vertex_index);
			}
		}
		else if (vol_type == VolumeType::TriangularPrism) // prism case
		{
			vol = add_prism(static_cast<ParentMESH&>(m), 3u, false);

			const std::array<Dart, 6> vertices_of_prism = {vol.dart,
														   phi1(m, vol.dart),
														   phi_1(m, vol.dart),
														   phi<-1, 2, 1, 1, 2>(m, vol.dart),
														   phi<2, 1, 1, 2>(m, vol.dart),
														   phi<1, 2, 1, 1, 2>(m, vol.dart)};

			for (Dart dv : vertices_of_prism)
			{
				const uint32 vertex_index =
					volume_data.vertex_id_after_import_[volume_data.volumes_vertex_indices_[index++]];
				foreach_dart_of_orbit(m, Vertex2(dv), [&](Dart d) -> bool {
					set_index<Vertex>(m, d, vertex_index);
					return true;
				});
				mark_n_store(dv, vertex_index);
			}
		}
		else if (vol_type == VolumeType::Hexa) // hexahedral case
		{
			vol = add_prism(static_cast<ParentMESH&>(m), 4u, false);

			const std::array<Dart, 8> vertices_of_hexa = {vol.dart,
														  phi1(m, vol.dart),
														  phi<1, 1>(m, vol.dart),
														  phi_1(m, vol.dart),
														  phi<-1, 2, 1, 1, 2>(m, vol.dart),
														  phi<2, 1, 1, 2>(m, vol.dart),
														  phi<1, 2, 1, 1, 2>(m, vol.dart),
														  phi<1, 1, 2, 1, 1, 2>(m, vol.dart)};

			for (Dart dv : vertices_of_hexa)
			{
				const uint32 vertex_index =
					volume_data.vertex_id_after_import_[volume_data.volumes_vertex_indices_[index++]];
				foreach_dart_of_orbit(m, Vertex2(dv), [&](Dart d) -> bool {
					set_index<Vertex>(m, d, vertex_index);
					return true;
				});
				mark_n_store(dv, vertex_index);
			}
		}
		else // end of hexa
		{
			if (vol_type == VolumeType::Connector)
			{
				index += 4u;
				// The second part of the code generates connectors automatically. We don't have to do anything here.
			}
		}

		if (is_indexed<Volume>(m))
			set_index(m, vol, vol_emb++);
	}

	// reconstruct neighbourhood
	uint32 nb_boundary_faces = 0u;
	DartMarkerStore marker(m);

	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		if (phi3(m, d) == d && !marker.is_marked(d))
		{
			foreach_dart_of_orbit(m, Face2(d), [&](Dart fd) -> bool {
				marker.mark(fd);
				return true;
			});

			Dart good_dart;

			// 1st step : for every dart of the face we try to find a valid phi3 candidate.
			// If we can't it's a boundary face.
			Dart d_it = d;
			do
			{
				uint32 vindex1 = index_of(m, Vertex(d_it));
				uint32 vindex2 = index_of(m, Vertex(phi1(m, phi1(m, d_it))));
				const std::vector<Dart>& vec = value<std::vector<Dart>>(m, darts_per_vertex, Vertex(phi1(m, d_it)));
				for (auto it = vec.begin(); it != vec.end() && good_dart.is_nil(); ++it)
					if (index_of(m, Vertex(phi1(m, *it))) == vindex1 && index_of(m, Vertex(phi_1(m, *it))) == vindex2)
						good_dart = *it;
				d_it = phi1(m, d_it);
			} while (good_dart.is_nil() && (d_it != d));
			d = phi_1(m, d_it);

			if (!good_dart.is_nil())
			{
				// const uint32 degD = codegree(m, MESH::Face(d));
				// const uint32 degGD = codegree(m, MESH::Face(good_dart));

				// if (degD == degGD) // normal case : the two opposite faces have the same degree
				//{
				//	Dart it1 = d;
				//	Dart it2 = good_dart;
				//	do
				//	{
				//		phi3_sew(m, it1, it2);
				//		it1 = phi1(m, it1);
				//		it2 = phi_1(m, it2);
				//	} while (it1 != d);
				// }
				// else
				//{
				//	// there is one face of degree 4 and one face of degree 3
				//	// -> stamp volume
				// }

				std::vector<std::pair<Dart, Dart>> vd;
				vd.reserve(32);
				Dart it1 = d;
				Dart it2 = good_dart;
				bool stop1 = false;
				bool stop2 = false;
				do
				{
					vd.push_back(std::make_pair(it1, it2));
					it1 = phi1(m, it1);
					stop1 = (it1 == d);
					it2 = phi_1(m, it2);
					stop2 = (it2 == good_dart);
				} while (stop1 || stop2);

				if (stop1 && stop2) // normal case : the two opposite faces have the same degree
				{
					for (const auto& p : vd)
						phi3_sew(m, p.first, p.second);
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
	}

	if (nb_boundary_faces > 0u)
	{
		uint32 nb_holes = close(m);
		std::cout << nb_holes << " hole(s) have been closed" << std::endl;
		std::cout << nb_boundary_faces << " boundary faces" << std::endl;
	}

	remove_attribute<Vertex>(m, darts_per_vertex);
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_VOLUME_IMPORT_H_
