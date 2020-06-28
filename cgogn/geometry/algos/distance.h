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

#ifndef CGOGN_GEOMETRY_ALGOS_DISTANCE_H_
#define CGOGN_GEOMETRY_ALGOS_DISTANCE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/thread_pool.h>

#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
Vec3 closest_point_on_surface(const MESH& m,
							  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							  const Vec3& p)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	Vec3 closest(0, 0, 0);
	Scalar min_dist = std::numeric_limits<Scalar>::max();

	foreach_cell(m, [&](Face f) -> bool {
		std::vector<Vertex> vertices = incident_vertices(m, f);
		// std::vector<const Vec3*> vertices_position;
		// std::transform(vertices.begin(), vertices.end(), std::back_inserter(vertices_position),
		// 			   [&](Vertex v) -> const Vec3* { return &value<Vec3>(m, vertex_position, v); });
		Scalar u, v, w;
		// assume triangle faces
		const Vec3& a = value<Vec3>(m, vertex_position, vertices[0]);
		const Vec3& b = value<Vec3>(m, vertex_position, vertices[1]);
		const Vec3& c = value<Vec3>(m, vertex_position, vertices[2]);
		closest_point_in_triangle(p, a, b, c, u, v, w);
		Vec3 pos = u * a + v * b + w * c;
		Scalar dist = (pos - p).squaredNorm();
		if (dist < min_dist)
		{
			closest = pos;
			min_dist = dist;
		}
		return true;
	});

	return closest;
}

template <typename MESH, typename GRID>
Vec3 closest_point_on_surface(const MESH& m,
							  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							  const GRID& g, const Vec3& p)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	Vec3 closest(0, 0, 0);
	Scalar min_dist = std::numeric_limits<Scalar>::max();

	g.foreach_face_around(p, [&](Face f) {
		std::vector<Vertex> vertices = incident_vertices(m, f);
		// std::vector<const Vec3*> vertices_position;
		// std::transform(vertices.begin(), vertices.end(), std::back_inserter(vertices_position),
		// 			   [&](Vertex v) -> const Vec3* { return &value<Vec3>(m, vertex_position, v); });
		Scalar u, v, w;
		// assume triangle faces
		const Vec3& a = value<Vec3>(m, vertex_position, vertices[0]);
		const Vec3& b = value<Vec3>(m, vertex_position, vertices[1]);
		const Vec3& c = value<Vec3>(m, vertex_position, vertices[2]);
		closest_point_in_triangle(p, a, b, c, u, v, w);
		Vec3 pos = u * a + v * b + w * c;
		Scalar dist = (pos - p).squaredNorm();
		if (dist < min_dist)
		{
			closest = pos;
			min_dist = dist;
		}
	});

	return closest;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_DISTANCE_H_
