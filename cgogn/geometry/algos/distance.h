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

	std::vector<std::tuple<Vec3, Scalar>> closest_per_thread(
		thread_pool()->nb_workers(), {{0, 0, 0}, std::numeric_limits<Scalar>::max()}); // position, squaredDistance

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 worker_index = current_worker_index();
		std::vector<Vertex> vertices = incident_vertices(m, f);
		std::vector<const Vec3*> vertices_position;
		std::transform(vertices.begin(), vertices.end(), std::back_inserter(vertices_position),
					   [&](Vertex v) -> const Vec3* { return &value<Vec3>(m, vertex_position, v); });
		Scalar u, v, w;
		// assume triangle faces
		closest_point_in_triangle(p, *(vertices_position[0]), *(vertices_position[1]), *(vertices_position[2]), u, v,
								  w);
		Vec3 pos = u * *(vertices_position[0]) + v * *(vertices_position[1]) + w * *(vertices_position[2]);
		Scalar dist = (pos - p).squaredNorm();
		if (dist < std::get<1>(closest_per_thread[worker_index]))
			closest_per_thread[worker_index] = {pos, dist};
		return true;
	});

	Vec3 closest_p(0, 0, 0);
	Scalar min_dist = std::numeric_limits<Scalar>::max();
	for (auto& [pos, dist] : closest_per_thread)
	{
		if (dist < min_dist)
		{
			min_dist = dist;
			closest_p = pos;
		}
	}

	return closest_p;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_DISTANCE_H_
