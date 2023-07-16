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

#ifndef CGOGN_GEOMETRY_ALGOS_CENTROID_H_
#define CGOGN_GEOMETRY_ALGOS_CENTROID_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

template <typename MESH>
struct mesh_traits;

namespace geometry
{

template <typename VEC, typename CELL, typename MESH>
VEC centroid(const MESH& m, CELL c, const typename mesh_traits<MESH>::template Attribute<VEC>* vertex_attribute)
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Scalar = typename vector_traits<VEC>::Scalar;
	VEC result;
	if constexpr (vector_traits<VEC>::SIZE == 1)
		result = 0;
	else
		result.setZero();
	uint32 count = 0;
	foreach_incident_vertex(m, c, [&](Vertex v) -> bool {
		result += value<VEC>(m, vertex_attribute, v);
		++count;
		return true;
	});
	result /= Scalar(count);
	return result;
}

template <typename VEC, typename MESH>
VEC centroid(const MESH& m, const typename mesh_traits<MESH>::template Attribute<VEC>* vertex_attribute)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Scalar = typename vector_traits<VEC>::Scalar;
	VEC result;
	result.setZero();
	uint32 count = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		result += value<VEC>(m, vertex_attribute, v);
		++count;
		return true;
	});
	result /= Scalar(count);
	return result;
}

template <typename VEC, typename CELL, typename MESH>
void compute_centroid(const MESH& m, const typename mesh_traits<MESH>::template Attribute<VEC>* vertex_attribute,
					  typename mesh_traits<MESH>::template Attribute<VEC>* cell_centroid)
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");

	parallel_foreach_cell(m, [&](CELL c) -> bool {
		value<VEC>(m, cell_centroid, c) = centroid<VEC>(m, c, vertex_attribute);
		return true;
	});
}

template <typename VEC, typename MESH>
typename mesh_traits<MESH>::Vertex central_vertex(const MESH& m,
												  const typename mesh_traits<MESH>::template Attribute<VEC>* attribute)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Scalar = typename vector_traits<VEC>::Scalar;
	VEC center = centroid<VEC>(m, attribute);
	Scalar min_dist = std::numeric_limits<Scalar>::max();
	Vertex min_vertex;
	foreach_cell(m, [&](Vertex v) -> bool {
		Scalar distance = (value(m, attribute, v) - center).squaredNorm();
		if (distance < min_dist)
		{
			min_dist = distance;
			min_vertex = v;
		}
		return true;
	});
	return min_vertex;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_CENTROID_H_
