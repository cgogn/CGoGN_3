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

#ifndef CGOGN_GEOMETRY_ALGOS_LENGTH_H_
#define CGOGN_GEOMETRY_ALGOS_LENGTH_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/vector_ops.h>

namespace cgogn
{

namespace geometry
{

template <typename VEC, typename MESH>
typename vector_traits<VEC>::Scalar
length(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	const typename mesh_traits<MESH>::template AttributePtr<VEC> vertex_position
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, e);
	return norm(value<VEC>(m, vertex_position, vertices[0]) - value<VEC>(m, vertex_position, vertices[1]));
}

template <typename VEC, typename MESH>
typename vector_traits<VEC>::Scalar
mean_edge_length(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<VEC> vertex_position
)
{
	using Scalar = typename vector_traits<VEC>::Scalar;
	using Edge = typename mesh_traits<MESH>::Edge;
	Scalar length_sum = 0;
	uint32 nbe = 0;
	foreach_cell(m, [&] (Edge e) -> bool
	{
		length_sum += length<VEC>(m, e, vertex_position);
		++nbe;
		return true;
	});
	return length_sum / Scalar(nbe);
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_LENGTH_H_
