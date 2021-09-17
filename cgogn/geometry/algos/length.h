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

#ifndef CGOGN_GEOMETRY_ALGOS_LENGTH_H_
#define CGOGN_GEOMETRY_ALGOS_LENGTH_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
Scalar length(const MESH& m, typename mesh_traits<MESH>::Edge e,
			  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, e);
	return (value<Vec3>(m, vertex_position, vertices[0]) - value<Vec3>(m, vertex_position, vertices[1])).norm();
}

template <typename MESH>
Scalar squared_length(const MESH& m, typename mesh_traits<MESH>::Edge e,
					  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, e);
	return (value<Vec3>(m, vertex_position, vertices[0]) - value<Vec3>(m, vertex_position, vertices[1])).squaredNorm();
}

template <typename MESH>
Scalar mean_edge_length(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Edge = typename mesh_traits<MESH>::Edge;

	thread_pool()->execute_all([]() {
		float64_value() = 0.0;
		uint32_value() = 0;
	});

	parallel_foreach_cell(m, [&](Edge e) -> bool {
		float64_value() += length(m, e, vertex_position);
		++uint32_value();
		return true;
	});

	Scalar length_sum = 0.0;
	uint32 nbe = 0;

	std::mutex mutex;
	thread_pool()->execute_all([&]() {
		std::lock_guard<std::mutex> lock(mutex);
		length_sum += float64_value();
		nbe += uint32_value();
	});

	return length_sum / Scalar(nbe);
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_LENGTH_H_
