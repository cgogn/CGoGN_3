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

#include <cgogn/io/surface/surface_import.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/face.h>

#include <vector>

namespace cgogn
{

namespace io
{

void import_surface_data(CMap2& m, const SurfaceImportData& surface_data)
{
	using Vertex = CMap2::Vertex;

	auto darts_per_vertex = add_attribute<std::vector<Dart>, Vertex>(m, "__darts_per_vertex");

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(16u);

	for (uint32 i = 0u, end = surface_data.faces_nb_vertices_.size(); i < end; ++i)
	{
		uint32 nbv = surface_data.faces_nb_vertices_[i];

		vertices_buffer.clear();
		uint32 prev = std::numeric_limits<uint32>::max();

		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = surface_data.faces_vertex_indices_[faces_vertex_index++];
			if (idx != prev)
			{
				prev = idx;
				vertices_buffer.push_back(idx);
			}
		}
		if (vertices_buffer.front() == vertices_buffer.back())
			vertices_buffer.pop_back();

		nbv = vertices_buffer.size();
		if (nbv > 2u)
		{
			CMap1::Face f = add_face(static_cast<CMap1&>(m), nbv, false);
			Dart d = f.dart;
			for (uint32 j = 0u; j < nbv; ++j)
			{
				const uint32 vertex_index = vertices_buffer[j];
				set_index<Vertex>(m,d, vertex_index);
				(*darts_per_vertex)[vertex_index].push_back(d);
				d = phi1(m,d);
			}
		}
	}

	bool need_vertex_unicity_check = false;
	uint32 nb_boundary_edges = 0u;

    foreach_dart(*m.mesh(),[&] (Dart d) -> bool
	{
		if (phi2(m,d) == d)
		{
			uint32 vertex_index = index_of(m,Vertex(d));

			const std::vector<Dart>& next_vertex_darts = value<std::vector<Dart>>(m, darts_per_vertex, Vertex(phi1(m,d)));
			bool phi2_found = false;
			bool first_OK = true;

			for (auto it = next_vertex_darts.begin();
				 it != next_vertex_darts.end() && !phi2_found;
				 ++it)
			{
				if (index_of(m,Vertex(phi1(m,*it))) == vertex_index)
				{
					if (phi2(m,*it) == *it)
					{
						phi2_sew(m,d, *it);
						phi2_found = true;
					}
					else
						first_OK = false;
				}
			}

			if (!phi2_found)
				++nb_boundary_edges;

			if (!first_OK)
				need_vertex_unicity_check = true;
		}
		return true;
	});

	if (nb_boundary_edges > 0u)
	{
		uint32 nb_holes = close(m);
		std::cout << nb_holes << " hole(s) have been closed" << std::endl;
		std::cout << nb_boundary_edges << " boundary edges" << std::endl;
	}

	// if (need_vertex_unicity_check)
	// {
	// 	map_.template enforce_unique_orbit_embedding<Vertex::ORBIT>();
	// 	cgogn_log_warning("create_map") << "Import Surface: non manifold vertices detected and corrected";
	// }

	remove_attribute<Vertex>(m, darts_per_vertex);
}

} // namespace io

} // namespace cgogn
