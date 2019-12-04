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

#include <cgogn/io/volume/volume_import.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

#include <vector>

namespace cgogn
{

namespace io
{

void import_volume_data(CMap3& m, const VolumeImportData& volume_data)
{
	using Vertex = CMap3::Vertex;
	using Volume = CMap3::Volume;

	auto darts_per_vertex = add_attribute<std::vector<Dart>, Vertex>(m, "__darts_per_vertex");
	
	uint32 index = 0u;
    DartMarker dart_marker(m);
	uint32 vol_emb = 0u;

	// for each volume of table
	for (uint32 i = 0u, end = volume_data.volumes_types_.size(); i < end; ++i)
	{
		Volume vol;
		const VolumeType vol_type = volume_data.volumes_types_[i];

		if (vol_type == VolumeType::Tetra) // tetrahedral case
		{
			vol = add_pyramid(static_cast<CMap2&>(m), 3u, false);

			const std::array<Dart, 4> vertices_of_tetra = {
				vol.dart,
				phi1(m,vol.dart),
				phi_1(m,vol.dart),
				phi_1(m,phi2(m,phi_1(m,vol.dart)))
			};

			for (Dart dv : vertices_of_tetra)
			{
				const uint32 vertex_index = volume_data.volumes_vertex_indices_[index++];
				foreach_dart_of_orbit(m,CMap3::Vertex2(dv), [&] (Dart d) -> bool
				{
					set_index<Vertex>(m,d, vertex_index);
					return true;
				});

				Dart dd = dv;
				do
				{
					dart_marker.mark(dd);
					(*darts_per_vertex)[vertex_index].push_back(dd);
					dd = phi1(m,phi2(m,dd));
				} while (dd != dv);
			}
		}
		else if (vol_type == VolumeType::Pyramid) // pyramidal case
		{
			vol = add_pyramid(static_cast<CMap2&>(m), 4u, false);

			const std::array<Dart, 5> vertices_of_pyramid = {
				vol.dart,
				phi1(m,vol.dart),
				phi1(m,phi1(m,vol.dart)),
				phi_1(m,vol.dart),
				phi_1(m,phi2(m,phi_1(m,vol.dart)))
			};

			for (Dart dv : vertices_of_pyramid)
			{
				const uint32 vertex_index = volume_data.volumes_vertex_indices_[index++];
				foreach_dart_of_orbit(m,CMap3::Vertex2(dv), [&] (Dart d) -> bool
				{
					set_index<Vertex>(m,d, vertex_index);
					return true;
				});

				Dart dd = dv;
				do
				{
					dart_marker.mark(dd);
					(*darts_per_vertex)[vertex_index].push_back(dd);
					dd = phi1(m,phi2(m,dd));
				} while (dd != dv);
			}
		}
		else if (vol_type == VolumeType::TriangularPrism) // prism case
		{
			vol = add_prism(static_cast<CMap2&>(m), 3u, false);

			const std::array<Dart, 6> vertices_of_prism = {
				vol.dart,
				phi1(m,vol.dart),
				phi_1(m,vol.dart),
				phi2(m,phi1(m,phi1(m,phi2(m,phi_1(m,vol.dart))))),
				phi2(m,phi1(m,phi1(m,phi2(m,vol.dart)))),
				phi2(m,phi1(m,phi1(m,phi2(m,phi1(m,vol.dart)))))
			};

			for (Dart dv : vertices_of_prism)
			{
				const uint32 vertex_index = volume_data.volumes_vertex_indices_[index++];
				foreach_dart_of_orbit(m,CMap3::Vertex2(dv), [&] (Dart d) -> bool
				{
					set_index<Vertex>(m,d, vertex_index);
					return true;
				});

				Dart dd = dv;
				do
				{
					dart_marker.mark(dd);
					(*darts_per_vertex)[vertex_index].push_back(dd);
					dd = phi1(m,phi2(m,dd));
				} while (dd != dv);
			}
		}
		else if (vol_type == VolumeType::Hexa) // hexahedral case
		{
			vol = add_prism(static_cast<CMap2&>(m), 4u,false);

			const std::array<Dart, 8> vertices_of_hexa = {
				vol.dart,
				phi1(m,vol.dart),
				phi1(m,phi1(m,vol.dart)),
				phi_1(m,vol.dart),
				phi2(m,phi1(m,phi1(m,phi2(m,phi_1(m,vol.dart))))),
				phi2(m,phi1(m,phi1(m,phi2(m,vol.dart)))),
				phi2(m,phi1(m,phi1(m,phi2(m,phi1(m,vol.dart))))),
				phi2(m,phi1(m,phi1(m,phi2(m,phi1(m,phi1(m,vol.dart))))))
			};

			for (Dart dv : vertices_of_hexa)
			{
				const uint32 vertex_index = volume_data.volumes_vertex_indices_[index++];
				auto tmp = m.mesh().attribute_containers_[Vertex::ORBIT].nb_refs(579);
				if(tmp == 0)
					std::cout << "hum" << std::endl;
				foreach_dart_of_orbit(static_cast<CMap2&>(m),CMap3::Vertex2(dv), [&] (Dart d) -> bool
				{
					set_index<Vertex>(m,d, vertex_index);
					std::cout << dv.index << " , " << d.index << std::endl;
					return true;
				});

				Dart dd = dv;
				do
				{
					dart_marker.mark(dd);
					(*darts_per_vertex)[vertex_index].push_back(dd);
					dd = phi1(m,phi2(m,dd));
				} while (dd != dv);
			}
		}
		else //end of hexa
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

    foreach_dart(m,[&] (Dart d) -> bool
	{
		if (phi3(m,d) == d && !marker.is_marked(d))
		{
			foreach_dart_of_orbit(m,CMap3::Face2(d), [&] (Dart fd) -> bool { marker.mark(fd); return true; });

			Dart good_dart;

			// 1st step : for every dart of the face we try to find a valid phi3 candidate. If we can't it's a boundary face.
			Dart d_it = d;
			do
			{
				uint32 vindex1 = index_of(m,Vertex(d_it));
				uint32 vindex2 = index_of(m,Vertex(phi1(m,phi1(m,d_it))));
				const std::vector<Dart>& vec = value<std::vector<Dart>>(m, darts_per_vertex, Vertex(phi1(m,d_it)));
				for (auto it = vec.begin(); it != vec.end() && good_dart.is_nil(); ++it)
					if (index_of(m,Vertex(phi1(m,*it))) == vindex1 && index_of(m,Vertex(phi_1(m,*it))) == vindex2)
						good_dart = *it;
				d_it = phi1(m,d_it);
			} while (good_dart.is_nil() && (d_it != d));
			d = phi_1(m,d_it);

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
						phi3_sew(m,it1, it2);
						it1 = phi1(m,it1);
						it2 = phi_1(m,it2);
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
		uint32 nb_holes = close(m);
		std::cout << nb_holes << " hole(s) have been closed" << std::endl;
		std::cout << nb_boundary_faces << " boundary faces" << std::endl;
	}
	
	remove_attribute<Vertex>(m, darts_per_vertex);
}

} // namespace io

} // namespace cgogn
