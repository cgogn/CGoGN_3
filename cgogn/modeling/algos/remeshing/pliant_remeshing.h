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

#ifndef CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_
#define CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

#include <libacc/bvh_tree.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

///////////
// CMap2 //
///////////

inline void triangulate_incident_faces(CMap2& m, CMap2::Vertex v)
{
	std::vector<CMap2::Face> ifaces = incident_faces(m, v);
	for (CMap2::Face f : ifaces)
		cut_face(m, CMap2::Vertex(f.dart), CMap2::Vertex(phi<11>(m, f.dart)));
}

inline bool should_edge_flip(CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Vertex> iv = incident_vertices(m, e);
	const int32 w = degree(m, iv[0]);
	const int32 x = degree(m, iv[1]);
	const int32 y = degree(m, CMap2::Vertex(phi1(m, phi1(m, iv[0].dart))));
	const int32 z = degree(m, CMap2::Vertex(phi1(m, phi1(m, iv[1].dart))));

	// int32 flip = 0;
	// flip += w > 6 ? 1 : (w < 6 ? -1 : 0);
	// flip += x > 6 ? 1 : (x < 6 ? -1 : 0);
	// flip += y < 6 ? 1 : (y > 6 ? -1 : 0);
	// flip += z < 6 ? 1 : (z > 6 ? -1 : 0);
	// return flip > 1;

	int32 dev_pre = abs(w - 6) + abs(x - 6) + abs(y - 6) + abs(z - 6);
	int32 dev_post = abs(w - 1 - 6) + abs(x - 1 - 6) + abs(y + 1 - 6) + abs(z + 1 - 6);
	return dev_pre - dev_post > 0;
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void pliant_remeshing(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	Scalar mean_edge_length = geometry::mean_edge_length(m, vertex_position);

	const Scalar squared_min_edge_length = Scalar(0.5625) * mean_edge_length * mean_edge_length; // 0.5625 = 0.75^2
	const Scalar squared_max_edge_length = Scalar(1.5625) * mean_edge_length * mean_edge_length; // 1.5625 = 1.25^2

	CellCache<MESH> cache(m);

	auto bvh_vertex_index = add_attribute<uint32, Vertex>(m, "__bvh_vertex_index");
	std::vector<Vec3> bvh_vertex_position;
	bvh_vertex_position.reserve(nb_cells<Vertex>(m));
	uint32 idx = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, bvh_vertex_index, v) = idx++;
		bvh_vertex_position.push_back(value<Vec3>(m, vertex_position, v));
		return true;
	});
	std::vector<uint32> face_vertex_indices;
	face_vertex_indices.reserve(nb_cells<Face>(m) * 3);
	foreach_cell(m, [&](Face f) -> bool {
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			face_vertex_indices.push_back(value<uint32>(m, bvh_vertex_index, v));
			return true;
		});
		return true;
	});
	acc::BVHTree<uint32, Vec3>* surface_bvh = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, bvh_vertex_position);

	// cut long edges (and adjacent faces)
	cache.template build<Edge>();
	bool has_long_edge = false;
	do
	{
		has_long_edge = false;
		foreach_cell(cache, [&](Edge e) -> bool {
			if (geometry::squared_length(m, e, vertex_position) > squared_max_edge_length)
			{
				has_long_edge = true;
				std::vector<Vertex> iv = incident_vertices(m, e);
				Vertex v = cut_edge(m, e);
				value<Vec3>(m, vertex_position, v) =
					(value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
				triangulate_incident_faces(m, v);
			}
			return true;
		});
		cache.template build<Edge>();
	} while (has_long_edge);

	// collapse short edges
	cache.template build<Edge>();
	bool has_short_edge = false;
	do
	{
		has_short_edge = false;
		foreach_cell(cache, [&](Edge e) -> bool {
			if (geometry::squared_length(m, e, vertex_position) < squared_min_edge_length)
			{
				std::vector<Vertex> iv = incident_vertices(m, e);
				bool collapse = true;
				const Vec3& p = value<Vec3>(m, vertex_position, iv[0]);
				foreach_adjacent_vertex_through_edge(m, iv[1], [&](Vertex v) -> bool {
					const Vec3& vec = p - value<Vec3>(m, vertex_position, v);
					if (vec.squaredNorm() > squared_max_edge_length)
						collapse = false;
					return true;
				});
				if (collapse && edge_can_collapse(m, e))
				{
					has_short_edge = true;
					Vertex cv = collapse_edge(m, e);
					value<Vec3>(m, vertex_position, cv) = p;
				}
			}
			return true;
		});
		cache.template build<Edge>();
	} while (has_short_edge);

	// equalize valences with edge flips
	cache.template build<Edge>();
	foreach_cell(cache, [&](Edge e) -> bool {
		if (should_edge_flip(m, e))
			flip_edge(m, e);
		return true;
	});

	// tangential relaxation
	auto vertex_barycenter = add_attribute<Vec3, Vertex>(m, "__vertex_barycenter");
	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		Vec3 q(0, 0, 0);
		uint32 count = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
			q += value<Vec3>(m, vertex_position, av);
			++count;
			return true;
		});
		q /= Scalar(count);
		value<Vec3>(m, vertex_barycenter, v) = q;
		return true;
	});
	foreach_cell(m, [&](Vertex v) -> bool {
		Vec3 n = geometry::normal(m, v, vertex_position);
		const Vec3& q = value<Vec3>(m, vertex_barycenter, v);
		value<Vec3>(m, vertex_position, v) = q + n.dot(value<Vec3>(m, vertex_position, v) - q) * n;
		return true;
	});

	// project back on surface
	foreach_cell(m, [&](Vertex v) -> bool {
		value<Vec3>(m, vertex_position, v) = surface_bvh->closest_point(value<Vec3>(m, vertex_position, v));
		return true;
	});

	remove_attribute<Vertex>(m, bvh_vertex_index);
	remove_attribute<Vertex>(m, vertex_barycenter);
	delete surface_bvh;
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_
