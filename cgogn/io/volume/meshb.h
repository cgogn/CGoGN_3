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

#ifndef CGOGN_IO_VOLUME_MESHB_H_
#define CGOGN_IO_VOLUME_MESHB_H_

extern "C"
{
#include <thirdparty/libMeshb/libmeshb.h>
}

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
bool import_MESHB(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 3, "MESH dimension should be 3");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	VolumeImportData volume_data;

	int32 version, dimension;
	int64 mesh_index = GmfOpenMesh(filename.c_str(), GmfRead, &version, &dimension);
	if (mesh_index == 0)
		return false;

	bool use_floats = (version == GmfFloat);

	const uint32 number_of_vertices = uint32(GmfStatKwd(mesh_index, GmfVertices));
	const uint32 number_of_tetras = uint32(GmfStatKwd(mesh_index, GmfTetrahedra));
	const uint32 number_of_hexas = uint32(GmfStatKwd(mesh_index, GmfHexahedra));
	const uint32 number_of_prisms = uint32(GmfStatKwd(mesh_index, GmfPrisms));
	const uint32 number_of_pyramids = uint32(GmfStatKwd(mesh_index, GmfPyramids));

	uint32 nb_volumes = number_of_tetras + number_of_hexas + number_of_prisms + number_of_pyramids;

	if (number_of_vertices == 0u || nb_volumes == 0u)
	{
		std::cerr << "Error while reading the file \"" << filename << "\".";
		GmfCloseMesh(mesh_index);
		return false;
	}

	volume_data.reserve(number_of_vertices, nb_volumes);
	auto position = add_attribute<geometry::Vec3, Vertex>(m, "position");
	GmfGotoKwd(mesh_index, GmfVertices);
	int32 ref;

	// read vertices position
	if (use_floats)
		for (uint32 i = 0u; i < number_of_vertices; ++i)
		{
			uint32 idx = new_index<Vertex>(m);
			std::array<float32, 3> v;
			(void)GmfGetLin(mesh_index, GmfVertices, &v[0], &v[1], &v[2], &ref);
			position->operator[](idx)[0] = v[0];
			position->operator[](idx)[1] = v[1];
			position->operator[](idx)[2] = v[2];
			volume_data.vertices_id_.push_back(idx);
		}
	else
		for (uint32 i = 0u; i < number_of_vertices; ++i)
		{
			uint32 idx = new_index<Vertex>(m);
			std::array<float64, 3> v;
			(void)GmfGetLin(mesh_index, GmfVertices, &v[0], &v[1], &v[2], &ref);
			position->operator[](idx)[0] = v[0];
			position->operator[](idx)[1] = v[1];
			position->operator[](idx)[2] = v[2];
			volume_data.vertices_id_.push_back(idx);
		}

	// read volumes
	if (number_of_tetras > 0)
	{
		GmfGotoKwd(mesh_index, GmfTetrahedra);
		std::array<int, 4> ids;
		for (uint32 i = 0; i < number_of_tetras; ++i)
		{
			(void)GmfGetLin(mesh_index, GmfTetrahedra, &ids[0], &ids[1], &ids[2], &ids[3], &ref);
			for (auto& id : ids)
			{
				--id;
				id = volume_data.vertices_id_[id];
			}
			if (geometry::test_orientation_3D((*position)[ids[0]], (*position)[ids[1]], (*position)[ids[2]],
											  (*position)[ids[3]]) == geometry::Orientation3D::UNDER)
				std::swap(ids[1], ids[2]);
			volume_data.volumes_types_.push_back(VolumeType::Tetra);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
		}
	}

	if (number_of_hexas > 0)
	{
		GmfGotoKwd(mesh_index, GmfHexahedra);
		std::array<int, 8> ids;
		for (uint32 i = 0; i < number_of_hexas; ++i)
		{
			(void)GmfGetLin(mesh_index, GmfHexahedra, &ids[0], &ids[1], &ids[2], &ids[3], &ids[4], &ids[5], &ids[6],
							&ids[7], &ref);
			for (auto& id : ids)
			{
				--id;
				id = volume_data.vertices_id_[id];
			}
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
		}
	}

	if (number_of_prisms > 0)
	{
		GmfGotoKwd(mesh_index, GmfPrisms);
		std::array<int, 6> ids;
		for (uint32 i = 0; i < number_of_prisms; ++i)
		{
			(void)GmfGetLin(mesh_index, GmfPrisms, &ids[0], &ids[1], &ids[2], &ids[3], &ids[4], &ids[5], &ref);
			for (auto& id : ids)
			{
				--id;
				id = volume_data.vertices_id_[id];
			}
			if (geometry::test_orientation_3D((*position)[ids[3]], (*position)[ids[0]], (*position)[ids[1]],
											  (*position)[ids[2]]) == geometry::Orientation3D::OVER)
			{
				std::swap(ids[1], ids[2]);
				std::swap(ids[4], ids[5]);
			}
			volume_data.volumes_types_.push_back(VolumeType::TriangularPrism);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
		}
	}

	if (number_of_pyramids > 0)
	{
		GmfGotoKwd(mesh_index, GmfPyramids);
		std::array<int, 5> ids;
		for (uint32 i = 0; i < number_of_pyramids; ++i)
		{
			(void)GmfGetLin(mesh_index, GmfPyramids, &ids[0], &ids[1], &ids[2], &ids[3], &ids[4], &ref);
			for (auto& id : ids)
			{
				--id;
				id = volume_data.vertices_id_[id];
			}
			if (geometry::test_orientation_3D((*position)[ids[4]], (*position)[ids[0]], (*position)[ids[1]],
											  (*position)[ids[2]]) == geometry::Orientation3D::OVER)
				std::swap(ids[1], ids[3]);
			volume_data.volumes_types_.push_back(VolumeType::Pyramid);
			volume_data.volumes_vertex_indices_.insert(volume_data.volumes_vertex_indices_.end(), ids.begin(),
													   ids.end());
		}
	}

	import_volume_data(m, volume_data);

	return true;
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_VOLUME_MESHB_H_
