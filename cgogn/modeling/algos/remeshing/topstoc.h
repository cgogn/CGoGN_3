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

#ifndef CGOGN_MODELING_ALGOS_TOPSTOC_H_
#define CGOGN_MODELING_ALGOS_TOPSTOC_H_

#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/io/surface/surface_import.h>

#include <cgogn/modeling/algos/subdivision.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;

template <typename MESH>
void debug_off(MESH& _m)
{
	using Face = typename cgogn::mesh_traits<MESH>::Face;
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;

	typename mesh_traits<MESH>::template Attribute<Vec3>* position = get_attribute<Vec3, Vertex>(_m, "position").get();
	typename mesh_traits<MESH>::template Attribute<float32>* charac =
		get_attribute<float32, Vertex>(_m, "charac").get();

	std::ofstream file;
	file.open("debug.off");
	uint nb_vertices = 0;
	foreach_cell(_m, [&](Vertex v) -> bool {
		nb_vertices++;
		return true;
	});
	uint nb_faces = 0;
	foreach_cell(_m, [&](Face f) -> bool {
		nb_faces++;
		return true;
	});

	file << "COFF" << std::endl << nb_vertices << " " << nb_faces << " 0" << std::endl;

	std::unordered_map<uint32, Vertex> vertices;
	foreach_cell(_m, [&](Vertex v) -> bool { vertices[_m.index_of(v)] = v; });

	float32 mean_charac = 0;
	float32 min_charac = value<float32>(_m, charac, vertices[0]);
	float32 max_charac = value<float32>(_m, charac, vertices[0]);
	for (int i = 0; i < nb_vertices; ++i)
	{
		float32 c = value<float32>(_m, charac, vertices[i]);
		if (c < min_charac)
			min_charac = c;
		if (c > max_charac)
			max_charac = c;
		mean_charac += c;
	}
	mean_charac /= nb_vertices;

	for (int i = 0; i < nb_vertices; ++i)
	{
		max_charac = 1;
		min_charac = 0;
		Vec3 pos = value<Vec3>(_m, position, vertices[i]);
		float32 charac_value = value<float32>(_m, charac, vertices[i]);
		uint grayscale = (255 / (max_charac - min_charac)) * (charac_value - min_charac);
		file << pos.x() << " " << pos.y() << " " << pos.z() << " ";
		file << grayscale << " " << grayscale << " " << grayscale << std::endl;
	}

	foreach_cell(_m, [&](Face f) -> bool {
		std::vector<Vertex> iv = incident_vertices(_m, f);
		file << "3 " << _m.index_of(iv[0]) << " " << _m.index_of(iv[1]) << " " << _m.index_of(iv[2]) << std::endl;
		return true;
	});

	file.close();
}

template <typename MESH, typename Vertex>
void vertex_selection(MESH& _m, CellMarkerStore<MESH, Vertex>& cm_selected, uint32 _nb_vertices_to_keep,
					  float32 adapt_coef = 0.66)
{
	typename mesh_traits<MESH>::template Attribute<Vec3>* normal = get_attribute<Vec3, Vertex>(_m, "normal").get();
	typename mesh_traits<MESH>::template Attribute<float32>* charac =
		add_attribute<float32, Vertex>(_m, "charac").get();

	std::srand(std::time(nullptr)); // init rand

	float32 mean_charac_value = 0;
	float32 max_charac_value = 0;
	uint32 nb_cells = 0;
	std::vector<Vertex> vertices;
	// compute characteristic value
	foreach_cell(_m, [&](Vertex v) -> bool {
		nb_cells++;
		uint32 count = 0;
		float32 charac_value = 0;
		foreach_adjacent_vertex_through_edge(_m, v, [&](Vertex v2) -> bool {
			count++;

			Vec3 nv = value<Vec3>(_m, normal, v);
			Vec3 nv2 = value<Vec3>(_m, normal, v2);

			float32 dist = acos(nv.x() * nv2.x() + nv.y() * nv2.y() + nv.z() * nv2.z());
			charac_value = (1.0f - dist) / 2.0f;
			return true;
		});
		charac_value /= (float)count;

		if (max_charac_value < charac_value)
			max_charac_value = charac_value;

		mean_charac_value += charac_value;
		value<float32>(_m, charac, v) = charac_value;
		vertices.push_back(v);
		return true;
	});
	mean_charac_value /= (float)nb_cells;

	std::sort(vertices.begin(), vertices.end(), [&](Vertex v, Vertex v2) -> bool {
		return value<float32>(_m, charac, v) > value<float32>(_m, charac, v2);
	});
	uint32 count = 0;

	// stochastic vertex selection
	foreach_cell(_m, [&](Vertex v) -> bool {
		float32 random_variable = std::rand() / (float)(RAND_MAX);

		float32 distrib_func =
			(1.0f + adapt_coef * (value<float32>(_m, charac, v) / mean_charac_value - 1.0f)) / max_charac_value;

		if (count++ < _nb_vertices_to_keep)
		{
			cm_selected.mark(vertices.back());
			vertices.pop_back();
		}
		// if(random_variable < distrib_func){
		// 	cm_selected.mark(v);
		// }

		return true;
	});

	remove_attribute<Vertex>(_m, charac);
}

template <typename MESH, typename Vertex>
void region_growth(MESH& _m, typename mesh_traits<MESH>::template Attribute<uint32>* _vertex_anchor,
				   CellMarkerStore<MESH, Vertex>& cm_selected)
{
	// region growth cell cache
	std::list<Vertex> rg_cache;

	//*fill the rg cache
	for (int id : cm_selected.marked_cells())
	{
		foreach_cell(_m, [&](Vertex v) -> bool {
			// if cell is selected
			if (index_of(_m, v) == id)
			{
				// add to region growth list
				rg_cache.push_back(v);
				// set its anchor to its own index
				value<uint32>(_m, _vertex_anchor, v) = id;
				return true;
			}
			return true;
		});
	}

	//*compute the region growth
	// to ignore computed cells
	CellMarker<MESH, Vertex> cm_done(_m);
	// to ignore previously anchored cells
	CellMarker<MESH, Vertex> cm_anchored(_m);
	while (!rg_cache.empty())
	{
		if (!cm_done.is_marked(rg_cache.front()))
		{
			// for each vertex in the one-ring
			foreach_adjacent_vertex_through_edge(_m, rg_cache.front(), [&](Vertex v) -> bool {
				if (!cm_selected.is_marked(v) && !cm_done.is_marked(v) && !cm_anchored.is_marked(v))
				{
					// push vertex to compute it
					rg_cache.push_back(v);
					// anchor the ring vertex to the current vertex anchor
					value<uint32>(_m, _vertex_anchor, v) = value<uint32>(_m, _vertex_anchor, rg_cache.front());
					cm_anchored.mark(v);
				}
				return true;
			});
			cm_done.mark(rg_cache.front());
			// go to next vertex
			rg_cache.pop_front();
		}
		else
		{ // already computed
			// go to next vertex
			rg_cache.pop_front();
		}
	}
}

template <typename MESH, typename Vertex, typename Face>
void compute_surface_data(MESH& _m, MESH& _new_m,
						  typename mesh_traits<MESH>::template Attribute<Vec3>* _vertex_position,
						  typename mesh_traits<MESH>::template Attribute<Vec3>* _new_vertex_position,
						  typename mesh_traits<MESH>::template Attribute<uint32>* _vertex_anchor,
						  CellMarkerStore<MESH, Vertex>& cm_selected)
{
	cgogn::io::SurfaceImportData surface_data;

	// compute the new face count
	uint32 face_count = 0;
	foreach_cell(_m, [&](Face f) -> bool {
		std::vector<Vertex> iv = incident_vertices(_m, f);
		if (value<uint32>(_m, _vertex_anchor, iv[0]) != value<uint32>(_m, _vertex_anchor, iv[1]) &&
			value<uint32>(_m, _vertex_anchor, iv[0]) != value<uint32>(_m, _vertex_anchor, iv[2]) &&
			value<uint32>(_m, _vertex_anchor, iv[1]) != value<uint32>(_m, _vertex_anchor, iv[2]))
			face_count++;
		return true;
	});

	surface_data.reserve(cm_selected.marked_cells().size(), face_count);

	// map to link the old mesh vertices to the new ones
	std::unordered_map<uint32, uint32> vmap;

	uint32 vertex_id = 0;
	for (std::vector<uint32>::const_iterator it = cm_selected.marked_cells().begin();
		 it != cm_selected.marked_cells().end(); it++)
	{
		surface_data.nb_vertices_++;
		vmap.emplace(*it, vertex_id++);
		surface_data.vertex_position_.push_back((*_vertex_position)[*it]);
	}

	// push faces for surface import
	foreach_cell(_m, [&](Face f) -> bool {
		std::vector<Vertex> iv = incident_vertices(_m, f);
		if (value<uint32>(_m, _vertex_anchor, iv[0]) != value<uint32>(_m, _vertex_anchor, iv[1]) &&
			value<uint32>(_m, _vertex_anchor, iv[0]) != value<uint32>(_m, _vertex_anchor, iv[2]) &&
			value<uint32>(_m, _vertex_anchor, iv[1]) != value<uint32>(_m, _vertex_anchor, iv[2]))
		{
			surface_data.nb_faces_++;
			surface_data.faces_nb_vertices_.push_back(3);
			surface_data.faces_vertex_indices_.push_back(vmap[value<uint32>(_m, _vertex_anchor, iv[0])]);
			surface_data.faces_vertex_indices_.push_back(vmap[value<uint32>(_m, _vertex_anchor, iv[1])]);
			surface_data.faces_vertex_indices_.push_back(vmap[value<uint32>(_m, _vertex_anchor, iv[2])]);
		}
		return true;
	});

	// computes the mesh
	io::import_surface_data(_new_m, surface_data);
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void topstoc(ui::MeshProvider<MESH>* mp, MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
			 uint32 nb_vertices_to_keep)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;

	// new mesh creation
	std::string name = "simplified_" + mp->mesh_name(m);
	MESH* new_m = mp->add_mesh(name);

	CellMarkerStore<MESH, Vertex> cm_selected(m);

	// vertex selection
	vertex_selection(m, cm_selected, nb_vertices_to_keep);

	// region growth
	auto vertex_anchor = add_attribute<uint32, Vertex>(m, "anchor");
	region_growth(m, vertex_anchor.get(), cm_selected);

	// surface data
	auto new_vertex_position = add_attribute<geometry::Vec3, CMap2::Vertex>(*new_m, "position");
	compute_surface_data<MESH, Vertex, Face>(m, *new_m, vertex_position, new_vertex_position.get(), vertex_anchor.get(),
											 cm_selected);

	// finish and cleanup
	mp->set_mesh_bb_vertex_position(*new_m, new_vertex_position);
	remove_attribute<Vertex>(m, vertex_anchor);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_TOPSTOC_H_
