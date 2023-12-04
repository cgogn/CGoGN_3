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

#include <cgogn/io/surface/surface_import.h>

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/types/maps/gmap/gmap2.h>
#include <cgogn/core/types/triangle_soup/triangle_soup.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>

#include <algorithm>
#include <set>

namespace cgogn
{

namespace io
{

template <typename MESH, typename SURFDATA>
auto import_surface_data_map_tmpl(MESH& m, SURFDATA& surface_data, bool reconstruct_phi2)
	-> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	using ParentMESH = typename MESH::Parent;

	using Vertex = typename MESH::Vertex;
	using VEC = typename SURFDATA::VEC_TYPE;

	auto position = get_or_add_attribute<VEC, Vertex>(m, surface_data.vertex_position_attribute_name_);

	for (uint32 i = 0u; i < surface_data.nb_vertices_; ++i)
	{
		uint32 vertex_id = new_index<Vertex>(m);
		(*position)[vertex_id] = surface_data.vertex_position_[i];
		surface_data.vertex_id_after_import_.push_back(vertex_id);
	}

	if (!reconstruct_phi2)
	{
		uint32 faces_vertex_index = 0u;
		std::vector<uint32> vertices_buffer;
		vertices_buffer.reserve(16u);

		for (uint32 i = 0u; i < surface_data.nb_faces_; ++i)
		{
			uint32 nbv = surface_data.faces_nb_vertices_[i];
			vertices_buffer.clear();
			for (uint32 j = 0u; j < nbv; ++j)
			{
				uint32 idx =
					surface_data.vertex_id_after_import_[surface_data.faces_vertex_indices_[faces_vertex_index++]];
				vertices_buffer.push_back(idx);
			}

			nbv = uint32(vertices_buffer.size());
			typename ParentMESH::Face f = add_face(static_cast<ParentMESH&>(m), nbv, false);
			Dart d = f.dart_;
			for (uint32 j = 0u; j < nbv; ++j)
			{
				const uint32 vertex_index = vertices_buffer[j];
				set_index<Vertex>(m, d, vertex_index);
				d = phi1(m, d);
			}
		}
		return;
	}

	auto darts_per_vertex = add_attribute<std::vector<Dart>, Vertex>(m, "__darts_per_vertex");

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(16u);

	for (uint32 i = 0u; i < surface_data.nb_faces_; ++i)
	{
		uint32 nbv = surface_data.faces_nb_vertices_[i];

		vertices_buffer.clear();
		uint32 prev = std::numeric_limits<uint32>::max();

		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = surface_data.vertex_id_after_import_[surface_data.faces_vertex_indices_[faces_vertex_index++]];
			if (idx != prev)
			{
				prev = idx;
				vertices_buffer.push_back(idx);
			}
		}
		if (vertices_buffer.front() == vertices_buffer.back())
			vertices_buffer.pop_back();

		nbv = uint32(vertices_buffer.size());
		if (nbv > 2u)
		{
			typename ParentMESH::Face f = add_face(static_cast<ParentMESH&>(m), nbv, false);
			Dart d = f.dart_;
			for (uint32 j = 0u; j < nbv; ++j)
			{
				const uint32 vertex_index = vertices_buffer[j];
				set_index<Vertex>(m, d, vertex_index);
				(*darts_per_vertex)[vertex_index].push_back(d);
				d = phi1(m, d);
			}
		}
	}

	if (!reconstruct_phi2)
		return;
	
	bool need_vertex_unicity_check = false;
	uint32 nb_boundary_edges = 0u;

	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		if (phi2(m, d) == d)
		{
			uint32 vertex_index = index_of(m, Vertex(d));

			const std::vector<Dart>& next_vertex_darts =
				value<std::vector<Dart>>(m, darts_per_vertex, Vertex(phi1(m, d)));
			bool phi2_found = false;
			bool first_OK = true;

			for (auto it = next_vertex_darts.begin(); it != next_vertex_darts.end() && !phi2_found; ++it)
			{
				if (index_of(m, Vertex(phi1(m, *it))) == vertex_index)
				{
					if (phi2(m, *it) == *it)
					{
						phi2_sew(m, d, *it);
						phi2_found = true;

						if constexpr (std::is_same_v<MESH, GMap2>)
						{
							copy_index<Vertex>(m, beta2(m, d), d);
							Dart dd = beta0(m, d);
							copy_index<Vertex>(m, dd, beta2(m, dd));
						}
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
	}

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

void import_surface_data(CMap2& m, SurfaceImportData& surface_data,bool reconstruct_phi2)
{
	import_surface_data_map_tmpl<CMap2>(m, surface_data, reconstruct_phi2);
}

void import_surface_data(GMap2& m, SurfaceImportData& surface_data, bool reconstruct_phi2)
{
	import_surface_data_map_tmpl<GMap2>(m, surface_data, reconstruct_phi2);
}



void import_surface_data(CMap2& m, SurfaceImportData2D& surface_data, bool reconstruct_phi2)
{
	import_surface_data_map_tmpl<CMap2, SurfaceImportData2D>(m, surface_data,reconstruct_phi2);
}

void import_surface_data(IncidenceGraph& ig, SurfaceImportData& surface_data)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	auto position = get_or_add_attribute<geometry::Vec3, Vertex>(ig, surface_data.vertex_position_attribute_name_);

	for (uint32 i = 0u; i < surface_data.nb_vertices_; ++i)
	{
		Vertex v = add_vertex(ig);
		(*position)[v.index_] = surface_data.vertex_position_[i];
		surface_data.vertex_id_after_import_.push_back(v.index_);
	}

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(16u);
	std::vector<Edge> face_edges;
	face_edges.reserve(16u);

	std::set<std::pair<uint32, uint32>> edges;

	for (uint32 i = 0u; i < surface_data.nb_faces_; ++i)
	{
		uint32 nbv = surface_data.faces_nb_vertices_[i];

		face_edges.clear();

		vertices_buffer.clear();
		uint32 prev = std::numeric_limits<uint32>::max();

		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = surface_data.vertex_id_after_import_[surface_data.faces_vertex_indices_[faces_vertex_index++]];
			if (idx != prev)
			{
				prev = idx;
				vertices_buffer.push_back(idx);
			}
		}
		if (vertices_buffer.front() == vertices_buffer.back())
			vertices_buffer.pop_back();

		nbv = uint32(vertices_buffer.size());
		if (nbv > 2u)
		{
			for (uint32 j = 0u; j < nbv; ++j)
			{
				uint32 vertex_index_1 = vertices_buffer[j];
				uint32 vertex_index_2 = vertices_buffer[(j + 1) % nbv];

				if (vertex_index_1 > vertex_index_2)
					std::swap(vertex_index_1, vertex_index_2);

				std::pair<uint32, uint32> e(vertex_index_1, vertex_index_2);
				if (edges.find(e) == edges.end())
				{
					Edge edge = add_edge(ig, e.first, e.second);
					face_edges.push_back(edge);
				}
			}

			add_face(ig, face_edges);
		}
	}
}

void import_surface_data(TriangleSoup& ts, SurfaceImportData& surface_data)
{
	using Vertex = TriangleSoup::Vertex;
	using Face = TriangleSoup::Face;

	auto position = get_or_add_attribute<geometry::Vec3, Vertex>(ts, surface_data.vertex_position_attribute_name_);

	for (uint32 i = 0u; i < surface_data.nb_vertices_; ++i)
	{
		uint32 vertex_id = new_index<Vertex>(ts);
		(*position)[vertex_id] = surface_data.vertex_position_[i];
		surface_data.vertex_id_after_import_.push_back(vertex_id);
	}

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(8u);

	for (uint32 i = 0u; i < surface_data.nb_faces_; ++i)
	{
		uint32 nbv = surface_data.faces_nb_vertices_[i];
		vertices_buffer.clear();
		uint32 prev = std::numeric_limits<uint32>::max();
		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = surface_data.vertex_id_after_import_[surface_data.faces_vertex_indices_[faces_vertex_index++]];
			if (idx != prev)
			{
				prev = idx;
				vertices_buffer.push_back(idx);
			}
		}
		if (vertices_buffer.front() == vertices_buffer.back())
			vertices_buffer.pop_back();

		nbv = uint32(vertices_buffer.size());
		if (nbv > 2u)
		{
			for (uint32 j = 1u; j < nbv - 1u; ++j)
				add_face(ts, vertices_buffer[0], vertices_buffer[j], vertices_buffer[j + 1]);
		}
	}
}

} // namespace io

} // namespace cgogn
