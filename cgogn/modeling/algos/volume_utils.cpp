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

#include <cgogn/modeling/algos/volume_utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/io/surface/surface_import.h>

namespace cgogn
{

namespace modeling
{

//////////
// CMap //
//////////

void extract_volume_surface(CMap3& m3, CMap3::Attribute<Vec3>* m3_vertex_position, CMap2& m2,
							CMap2::Attribute<Vec3>* /*m2_vertex_position*/,
							CMap2::Attribute<CMap3::Vertex>* m2_vertex_m3_vertex,
							CMap3::Attribute<CMap2::Vertex>* m3_vertex_m2_vertex)
{
	auto m2_vertex_index = add_attribute<uint32, CMap3::Vertex>(m3, "m2_vertex_index");

	cgogn::io::SurfaceImportData surface_data;

	uint32 vertex_id = 0;
	foreach_cell(m3, [&](CMap3::Vertex v3) -> bool {
		if (is_incident_to_boundary(m3, v3))
		{
			surface_data.nb_vertices_++;
			value<uint32>(m3, m2_vertex_index, v3) = vertex_id++;
			surface_data.vertex_position_.push_back(value<Vec3>(m3, m3_vertex_position, v3));
		}
		return true;
	});

	std::vector<uint32> indices;
	indices.reserve(4);
	foreach_cell(m3, [&](CMap3::Face f) -> bool {
		if (is_incident_to_boundary(m3, f))
		{
			surface_data.nb_faces_++;
			uint32 nbv = 0;
			foreach_incident_vertex(m3, f, [&](CMap3::Vertex v3) -> bool {
				++nbv;
				indices.push_back(value<uint32>(m3, m2_vertex_index, v3));
				return true;
			});

			surface_data.faces_nb_vertices_.push_back(nbv);
			surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
													  indices.end());
			indices.clear();
		}
		return true;
	});

	import_surface_data(m2, surface_data);

	if (m2_vertex_m3_vertex || m3_vertex_m2_vertex)
	{
		foreach_cell(m3, [&](CMap3::Vertex v3) -> bool {
			if (is_incident_to_boundary(m3, v3))
			{
				uint32 vertex_id = surface_data.vertex_id_after_import_[value<uint32>(m3, m2_vertex_index, v3)];
				if (m2_vertex_m3_vertex)
					(*m2_vertex_m3_vertex)[vertex_id] = v3;
				if (m3_vertex_m2_vertex)
					value<CMap2::Vertex>(m3, m3_vertex_m2_vertex, v3) = of_index<CMap2::Vertex>(m2, vertex_id);
			}
			return true;
		});
	}

	remove_attribute<CMap3::Vertex>(m3, m2_vertex_index);
}

} // namespace modeling

} // namespace cgogn
