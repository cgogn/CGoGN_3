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

#ifndef CGOGN_GEOMETRY_ALGOS_DECIMATION_H_
#define CGOGN_GEOMETRY_ALGOS_DECIMATION_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/decimation/QEM_helper.h>
#include <cgogn/modeling/algos/decimation/edge_approximator.h>
#include <cgogn/modeling/algos/decimation/edge_queue_update.h>

namespace cgogn
{

namespace modeling
{

using geometry::Vec3;
using geometry::Scalar;

/////////////
// GENERIC //
/////////////

template <typename MESH>
void decimate(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
			  uint32 nb_vertices_to_remove)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	CellQueue<Edge> edge_queue;
	using EdgeQueueInfo = typename CellQueue<Edge>::CellQueueInfo;
	auto edge_queue_info = add_attribute<EdgeQueueInfo, Edge>(m, "__decimate_edge_queue_info");

	DecimationQEM_Helper<MESH> helper(m, vertex_position);

	auto before = [&](Edge e) { helper.before_collapse(e); };
	auto approx = [&](Edge e) -> Vec3 { return mid_point(m, e, vertex_position); }; // helper.edge_optimal(e); };
	auto edge_cost = [&](Edge e) -> Scalar { return helper.edge_cost(e, approx(e)); };
	auto after = [&](Vertex v) { helper.after_collapse(v); };

	// auto before = [](Edge e) {};
	// auto approx = [&](Edge e) -> Vec3 { return mid_edge(m, e, vertex_position); };
	// auto edge_cost = [&](Edge e) -> Scalar { return geometry::length(m, e, vertex_position); };
	// auto after = [](Vertex v) {};

	foreach_cell(m, [&](Edge e) -> bool {
		update_edge_queue(m, e, edge_queue, edge_queue_info.get(), edge_cost);
		return true;
	});

	uint32 count = 0;
	for (auto it = edge_queue.begin(); it != edge_queue.end(); ++it)
	{
		auto it_e = *it;
		Vec3 newpos = approx(it_e);

		Edge e1, e2;
		pre_collapse(m, it_e, e1, e2, edge_queue, edge_queue_info.get());
		before(it_e);
		Vertex v = collapse_edge(m, it_e);
		value<Vec3>(m, vertex_position, v) = newpos;
		after(v);
		post_collapse(m, e1, e2, edge_queue, edge_queue_info.get(), edge_cost);

		++count;
		if (count >= nb_vertices_to_remove)
			break;
	}

	remove_attribute<Edge>(m, edge_queue_info);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_DECIMATION_H_
