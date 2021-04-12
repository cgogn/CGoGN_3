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

/////////////
// GENERIC //
/////////////

template <typename MESH>
void pliant_remeshing(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	Scalar mean_edge_length = geometry::mean_edge_length(m, vertex_position);

	const Scalar squared_min_edge_length = Scalar(0.5625) * mean_edge_length * mean_edge_length; // 0.5625 = 0.75^2
	const Scalar squared_max_edge_length = Scalar(1.5625) * mean_edge_length * mean_edge_length; // 1.5625 = 1.25^2

	CellCache<MESH> cache(m);

	// cut long edges (and adjacent faces)
	auto long_edges_selection = [&](Edge e) -> bool {
		return geometry::squared_length(m, e, vertex_position) > squared_max_edge_length;
	};

	cache.template build<Edge>(long_edges_selection);
	while (cache.template size<Edge>() > 0)
	{
		foreach_cell(cache, [&](Edge e) -> bool {
			std::vector<Vertex> iv = incident_vertices(m, e);
			Vertex v = cut_edge(m, e);
			value<Vec3>(m, vertex_position, v) =
				(value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
			triangulate_incident_faces(m, v);
			return true;
		});
		cache.template build<Edge>(long_edges_selection);
	}

	// collapse short edges
	auto short_edges_selection = [&](Edge e) -> bool {
		bool collapse = geometry::squared_length(m, e, vertex_position) < squared_min_edge_length;
		if (collapse)
		{
			std::vector<Vertex> iv = incident_vertices(m, e);
			const Vec3& p = value<Vec3>(m, vertex_position, iv[0]);
			foreach_adjacent_vertex_through_edge(m, iv[1], [&](Vertex v) -> bool {
				const Vec3& vec = p - value<Vec3>(m, vertex_position, v);
				if (vec.squaredNorm() > squared_max_edge_length)
					collapse = false;
				return true;
			});
			if (collapse)
			{
				if (!edge_can_collapse(m, e))
					collapse = false;
			}
		}
		return collapse;
	};
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_
