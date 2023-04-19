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

#ifndef CGOGN_GEOMETRY_ALGOS_SELECTION_H_
#define CGOGN_GEOMETRY_ALGOS_SELECTION_H_

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/inclusion.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

namespace cgogn
{

namespace geometry
{

CellCache<CMap2> within_sphere(const CMap2& m, typename CMap2::Vertex center, geometry::Scalar radius,
							   const typename CMap2::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename CMap2::Vertex;
	using HalfEdge = typename CMap2::HalfEdge;
	using Edge = typename CMap2::Edge;
	using Face = typename CMap2::Face;

	CellCache<CMap2> cache(m);

	const Vec3& center_position = value<Vec3>(m, vertex_position, center);

	DartMarkerStore dm(m);

	auto mark_vertex = [&](Vertex v) {
		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
			// mark a dart of the vertex
			dm.mark(d);

			// check if the edge of d is now completely marked
			// (which means all the vertices of the edge are in the sphere)
			Edge e(d);
			bool all_in = true;
			foreach_dart_of_orbit(m, e, [&](Dart dd) -> bool {
				if (!dm.is_marked(dd))
					all_in = false;
				return all_in;
			});
			if (all_in)
				cache.add(e);

			// check if the face of d is now completely marked
			// (which means all the vertices of the face are in the sphere)
			Face f(d);
			all_in = true;
			foreach_dart_of_orbit(m, f, [&](Dart dd) -> bool {
				if (!dm.is_marked(dd))
					all_in = false;
				return all_in;
			});
			if (all_in)
				cache.add(f);

			return true;
		});
	};

	cache.add(center);
	mark_vertex(center);

	uint32 i = 0;
	while (i < cache.cell_vector<Vertex>().size())
	{
		Vertex v = cache.cell_vector<Vertex>()[i];
		foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
			const Vec3& p = value<Vec3>(m, vertex_position, av);
			if (in_sphere(p, center_position, radius))
			{
				if (!dm.is_marked(av.dart))
				{
					cache.add(av);
					mark_vertex(av);
				}
			}
			else
				cache.add(HalfEdge(phi2(m, av.dart)));
			return true;
		});
		++i;
	}

	return cache;
}

template <typename CELL, typename MESH>
std::vector<CELL> within_normal_angle_threshold(
	const MESH& m, CELL start, geometry::Scalar angle_threshold,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "MESH dimension should be >= 2");
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");

	Vec3 n = normal(m, start, vertex_position);
	CellMarker<CMap2, CELL> marker(m);
	std::vector<CELL> cells;
	cells.push_back(start);

	auto check_cell = [&](CELL c) -> bool {
		if (!marker.is_marked(c))
		{
			Vec3 cn = normal(m, c, vertex_position);
			if (angle(n, cn) < angle_threshold)
			{
				marker.mark(c);
				cells.push_back(c);
			}
		}
		return true;
	};

	for (uint32 i = 0; i < cells.size(); ++i)
	{
		CELL c = cells[i];
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Vertex>)
			foreach_adjacent_vertex_through_edge(m, c, check_cell);
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
			foreach_adjacent_edge_through_face(m, c, check_cell);
		if constexpr (std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
			foreach_adjacent_face_through_edge(m, c, check_cell);
	}
	return cells;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_SELECTION_H_
