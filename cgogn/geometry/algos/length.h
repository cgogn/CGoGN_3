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

namespace cgogn
{

namespace geometry
{

template <typename MESH>
double
length(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3> vertex_position
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, e);
	return (value<Vec3>(m, vertex_position, vertices[0]) - value<Vec3>(m, vertex_position, vertices[1])).norm();
}

template <typename MESH>
double
mean_edge_length(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3> vertex_position
)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	double length_sum = 0.0;
	uint32 nbe = 0;
	foreach_cell(m, [&] (Edge e) -> bool
	{
		length_sum += length(m, e, vertex_position);
		++nbe;
		return true;
	});
	return length_sum / double(nbe);
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_LENGTH_H_
