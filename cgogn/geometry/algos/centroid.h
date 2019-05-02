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

#ifndef CGOGN_GEOMETRY_ALGOS_CENTROID_H_
#define CGOGN_GEOMETRY_ALGOS_CENTROID_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/vector_ops.h>

namespace cgogn
{

namespace geometry
{

template <typename VEC, typename CELL, typename MESH,
		  typename = typename std::enable_if<is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value>::type>
VEC
centroid(
	const MESH& m,
	CELL c,
	const typename mesh_traits<MESH>::template AttributePtr<VEC> attribute
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Scalar = typename vector_traits<VEC>::Scalar;
	VEC result;
	set_zero(result);
	uint32 count = 0;
	foreach_incident_vertex(m, c, [&] (Vertex v)
	{
		result += value<VEC>(m, attribute, v);
		++count;
	});
	result /= Scalar(count);
	return result;
}

template <typename VEC, typename MESH>
void
centroid(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<VEC> attribute
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Scalar = typename vector_traits<VEC>::Scalar;
	VEC result;
	set_zero(result);
	uint32 count = 0;
	foreach_cell(m, [&] (Vertex v)
	{
		result += value<VEC>(m, attribute, v);
		++count;
	});
	result /= Scalar(count);
	return result;
}

template <typename VEC, typename CELL, typename MESH,
		  typename = typename std::enable_if<is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value>::type>
void
compute_centroid(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<VEC> attribute,
	typename mesh_traits<MESH>::template AttributePtr<VEC> cell_centroid
)
{
	foreach_cell(m, [&] (CELL c)
	{
		value<VEC>(m, cell_centroid, c) = centroid<VEC>(m, c, attribute);
	});
}

template <typename VEC, typename MESH>
typename mesh_traits<MESH>::Vertex
central_vertex(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<VEC> attribute
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Scalar = typename vector_traits<VEC>::Scalar;
	VEC center = centroid<VEC>(m, attribute);
	Scalar min_dist = std::numeric_limits<Scalar>::max();
	Vertex min_vertex;
	foreach_cell(m, [&] (Vertex v)
	{
		Scalar distance = square_norm(value(m, attribute, v) - center);
		if (distance < min_dist)
		{
			min_dist = distance;
			min_vertex = v;
		}
	});
	return min_vertex;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_CENTROID_H_
