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

#ifndef CGOGN_IO_SURFACE_PLY_H_
#define CGOGN_IO_SURFACE_PLY_H_

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <thirdparty/happly/happly.h>

#include <vector>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_PLY(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	SurfaceImportData surface_data;

	happly::PLYData plyData(filename);
	std::vector<std::array<double, 3>> position = plyData.getVertexPositions();
	std::vector<std::vector<uint32>> face_indices = plyData.getFaceIndices<uint32>();

	const uint32 nb_vertices = position.size();
	const uint32 nb_faces = face_indices.size();

	surface_data.reserve(nb_vertices, nb_faces);

	for (uint32 i = 0u; i < nb_vertices; ++i)
	{
		const std::array<double, 3>& p = position[i];
		surface_data.vertex_position_.push_back({p[0], p[1], p[2]});
	}

	for (uint32 i = 0u; i < nb_faces; ++i)
	{
		surface_data.faces_nb_vertices_.push_back(face_indices[i].size());
		surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), face_indices[i].begin(),
												  face_indices[i].end());
	}

	import_surface_data(m, surface_data);

	return true;
}

template <typename MESH>
void export_PLY(MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
				const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	// TODO
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_PLY_H_
