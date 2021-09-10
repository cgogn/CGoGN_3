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

#include <cgogn/geometry/algos/angle.h>
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

inline bool edge_should_flip(CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Vertex> iv = incident_vertices(m, e);
	const int32 w = degree(m, iv[0]);
	const int32 x = degree(m, iv[1]);
	const int32 y = degree(m, CMap2::Vertex(phi1(m, phi1(m, iv[0].dart))));
	const int32 z = degree(m, CMap2::Vertex(phi1(m, phi1(m, iv[1].dart))));

	if (w < 4 || x < 4)
		return false;

	int32 dev_pre = abs(w - 6) + abs(x - 6) + abs(y - 6) + abs(z - 6);
	int32 dev_post = abs(w - 1 - 6) + abs(x - 1 - 6) + abs(y + 1 - 6) + abs(z + 1 - 6);
	return dev_pre - dev_post > 0;
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void pliant_remeshing(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
					  Scalar edge_length_target_ratio = 1.0)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	Scalar edge_length_target = geometry::mean_edge_length(m, vertex_position) * edge_length_target_ratio;

	const Scalar squared_min_edge_length = Scalar(0.5625) * edge_length_target * edge_length_target; // 0.5625 = 0.75^2
	const Scalar squared_max_edge_length = Scalar(1.5625) * edge_length_target * edge_length_target; // 1.5625 = 1.25^2

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

	auto feature_edge = add_attribute<bool, Edge>(m, "__feature_edge");
	feature_edge->fill(false);
	Scalar angle_threshold = 44.0 * M_PI / 180.0;
	foreach_cell(m, [&](Edge e) -> bool {
		if (std::fabs(geometry::angle(m, e, vertex_position)) > angle_threshold)
			value<bool>(m, feature_edge, e) = true;
		return true;
	});

	for (uint32 i = 0; i < 5; ++i)
	{
		// cut long edges (and adjacent faces)
		bool has_long_edge = false;
		do
		{
			cache.template build<Edge>();
			has_long_edge = false;
			foreach_cell(cache, [&](Edge e) -> bool {
				if (geometry::squared_length(m, e, vertex_position) > squared_max_edge_length)
				{
					has_long_edge = true;
					std::vector<Vertex> iv = incident_vertices(m, e);
					Vertex v = cut_edge(m, e);
					if (value<bool>(m, feature_edge, e))
					{
						foreach_incident_edge(m, v, [&](Edge ie) -> bool {
							value<bool>(m, feature_edge, ie) = true;
							return true;
						});
					}
					value<Vec3>(m, vertex_position, v) =
						(value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
					triangulate_incident_faces(m, v);
					foreach_incident_edge(m, v, [&](Edge ie) -> bool {
						if (std::fabs(geometry::angle(m, ie, vertex_position)) > angle_threshold)
							value<bool>(m, feature_edge, ie) = true;
						return true;
					});
				}
				return true;
			});
		} while (has_long_edge);

		// collapse short edges
		bool has_short_edge = false;
		do
		{
			has_short_edge = false;
			foreach_cell(m, [&](Edge e) -> bool {
				if (value<bool>(m, feature_edge, e))
					return true;
				if (geometry::squared_length(m, e, vertex_position) < squared_min_edge_length)
				{
					std::vector<Vertex> iv = incident_vertices(m, e);
					bool collapse = true;
					const Vec3& p = value<Vec3>(m, vertex_position, iv[0]);
					foreach_adjacent_vertex_through_edge(m, iv[1], [&](Vertex v) -> bool {
						const Vec3& vec = p - value<Vec3>(m, vertex_position, v);
						if (vec.squaredNorm() > squared_max_edge_length)
							collapse = false;
						return collapse;
					});
					foreach_incident_edge(m, iv[0], [&](Edge ie) -> bool {
						if (value<bool>(m, feature_edge, ie))
							collapse = false;
						return collapse;
					});
					if (collapse && edge_can_collapse(m, e))
					{
						has_short_edge = true;
						Vec3 mp =
							(value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
						Vertex cv = collapse_edge(m, e);
						value<Vec3>(m, vertex_position, cv) = mp;
					}
				}
				return true;
			});
		} while (has_short_edge);

		// equalize valences with edge flips
		foreach_cell(m, [&](Edge e) -> bool {
			if (!value<bool>(m, feature_edge, e) && edge_should_flip(m, e) && edge_can_flip(m, e))
				flip_edge(m, e);
			return true;
		});

		// tangential relaxation
		// + project back on surface
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			bool move_vertex = true;
			foreach_incident_edge(m, v, [&](Edge ie) -> bool {
				if (value<bool>(m, feature_edge, ie))
					move_vertex = false;
				return move_vertex;
			});
			if (move_vertex && !is_incident_to_boundary(m, v))
			{
				Vec3 q(0, 0, 0);
				uint32 count = 0;
				foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
					q += value<Vec3>(m, vertex_position, av);
					++count;
					return true;
				});
				q /= Scalar(count);
				Vec3 n = geometry::normal(m, v, vertex_position);
				Vec3 r = q + n.dot(value<Vec3>(m, vertex_position, v) - q) * n;
				value<Vec3>(m, vertex_position, v) = surface_bvh->closest_point(r);
			}
			return true;
		});
	}

	remove_attribute<Edge>(m, feature_edge);
	remove_attribute<Vertex>(m, bvh_vertex_index);
	delete surface_bvh;
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_
