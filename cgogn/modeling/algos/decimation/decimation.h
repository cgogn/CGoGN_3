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

#ifndef CGOGN_GEOMETRY_ALGOS_DECIMATION_H_
#define CGOGN_GEOMETRY_ALGOS_DECIMATION_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/mesh_ops/edge.h>

#include <cgogn/modeling/algos/decimation/edge_approximator.h>
#include <cgogn/modeling/algos/decimation/edge_queue_edge_length.h>

namespace cgogn
{

namespace modeling
{

/////////////
// GENERIC //
/////////////

template <typename VEC, typename MESH>
void decimate(MESH& m, typename mesh_traits<MESH>::template AttributePtr<VEC> vertex_position, uint32 nb_vertices_to_remove)
{
	using Scalar = typename geometry::vector_traits<VEC>::Scalar;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	CellQueue<Edge> edge_queue;
	using EdgeInfo = typename CellQueue<Edge>::CellInfo;
	auto edge_info = add_attribute<EdgeInfo, Edge>(m, "__decimate_edge_info");
	auto edge_cost = [&] (Edge e) -> Scalar
	{
		return geometry::length<VEC>(m, e, vertex_position);
	};

	foreach_cell(m, [&] (Edge e) -> bool
	{
		update_edge_queue(m, e, edge_queue, edge_info, edge_cost);
		return true;
	});

	uint32 count = 0;
	for (auto it = edge_queue.begin(); it != edge_queue.end(); ++it)
	{
		VEC newpos = mid_edge<VEC>(m, vertex_position, *it);

		Edge e1, e2;
		pre_collapse_edge_length(m, *it, e1, e2, edge_queue, edge_info);
		Vertex v = collapse_edge(m, *it);
		value<VEC>(m, vertex_position, v) = newpos;
		post_collapse_edge_length(m, e1, e2, edge_queue, edge_info, edge_cost);

		++count;
		if (count >= nb_vertices_to_remove)
			break;
	}

	remove_attribute<Edge>(m, edge_info);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_DECIMATION_H_
