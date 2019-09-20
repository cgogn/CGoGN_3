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

#ifndef CGOGN_MODELING_DECIMATION_EDGE_QUEUE_EDGE_LENGTH_H_
#define CGOGN_MODELING_DECIMATION_EDGE_QUEUE_EDGE_LENGTH_H_

#include <cgogn/modeling/algos/decimation/cell_queue.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/length.h>

#include <map>

namespace cgogn
{

namespace modeling
{

/////////////
// GENERIC //
/////////////

template <typename MESH, typename FUNC>
void
update_edge_queue(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	CellQueue<typename mesh_traits<MESH>::Edge>& edge_queue,
	typename mesh_traits<MESH>::template Attribute<typename CellQueue<typename mesh_traits<MESH>::Edge>::CellQueueInfo>* edge_queue_info,
	const FUNC& edge_cost
)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using EdgeQueueInfo = typename CellQueue<Edge>::CellQueueInfo;

	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Given function should take an Edge as parameter");
	static_assert(std::is_floating_point<func_return_type<FUNC>>::value, "Given function should return a floating point value");

	EdgeQueueInfo& ei = value<EdgeQueueInfo>(m, edge_queue_info, e);
	if (edge_can_collapse(m, e))
	{
		if (ei.valid_)
		{
			edge_queue.cells_.erase(ei.it_);
			ei.it_ = edge_queue.cells_.insert(std::make_pair(edge_cost(e), e));
		}
		else
		{
			ei.it_ = edge_queue.cells_.insert(std::make_pair(edge_cost(e), e));
			ei.valid_ = true;
		}
	}
	else
	{
		if (ei.valid_)
		{
			edge_queue.cells_.erase(ei.it_);
			ei.valid_ = false;
		}
	}
}

//////////////
// CMapBase //
//////////////

template <typename MESH,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
void
pre_collapse_edge_length(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	typename mesh_traits<MESH>::Edge& e1,
	typename mesh_traits<MESH>::Edge& e2,
	CellQueue<typename mesh_traits<MESH>::Edge>& edge_queue,
	typename mesh_traits<MESH>::template Attribute<typename CellQueue<typename mesh_traits<MESH>::Edge>::CellQueueInfo>* edge_queue_info
)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using EdgeQueueInfo = typename CellQueue<Edge>::CellQueueInfo;

	EdgeQueueInfo& ei = value<EdgeQueueInfo>(m, edge_queue_info, e);
	if (ei.valid_) { edge_queue.cells_.erase(ei.it_); ei.valid_ = false; }

	e1 = Edge(m.phi2(m.phi_1(e.dart)));
	e2 = Edge(m.phi2(m.phi_1(m.phi2(e.dart))));

	Dart ed1 = e.dart;
	Dart ed2 = m.phi2(ed1);

	EdgeQueueInfo& ei1 = value<EdgeQueueInfo>(m, edge_queue_info, Edge(m.phi1(ed1)));
	if (ei1.valid_) { edge_queue.cells_.erase(ei1.it_); ei1.valid_ = false; }

	EdgeQueueInfo& ei2 = value<EdgeQueueInfo>(m, edge_queue_info, Edge(m.phi_1(ed1)));
	if (ei2.valid_) { edge_queue.cells_.erase(ei2.it_); ei2.valid_ = false; }

	EdgeQueueInfo& ei3 = value<EdgeQueueInfo>(m, edge_queue_info, Edge(m.phi1(ed2)));
	if (ei3.valid_) { edge_queue.cells_.erase(ei3.it_); ei3.valid_ = false; }

	EdgeQueueInfo& ei4 = value<EdgeQueueInfo>(m, edge_queue_info, Edge(m.phi_1(ed2)));
	if (ei4.valid_) { edge_queue.cells_.erase(ei4.it_); ei4.valid_ = false; }
}

template <typename MESH, typename FUNC,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
void
post_collapse_edge_length(
	const MESH& m,
	typename mesh_traits<MESH>::Edge& e1,
	typename mesh_traits<MESH>::Edge& e2,
	CellQueue<typename mesh_traits<MESH>::Edge>& edge_queue,
	typename mesh_traits<MESH>::template Attribute<typename CellQueue<typename mesh_traits<MESH>::Edge>::CellQueueInfo>* edge_queue_info,
	const FUNC& edge_cost
)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using EdgeQueueInfo = typename CellQueue<Edge>::CellQueueInfo;

	Dart vit = e1.dart;
	do
	{
		update_edge_queue(m, Edge(m.phi1(vit)), edge_queue, edge_queue_info, edge_cost);
		if(vit == e1.dart || vit == e2.dart)
		{
			update_edge_queue(m, Edge(vit), edge_queue, edge_queue_info, edge_cost);

			Dart vit2 = m.template phi<121>(vit);
			Dart stop = m.phi2(vit);
			do
			{
				update_edge_queue(m, Edge(vit2), edge_queue, edge_queue_info, edge_cost);
				update_edge_queue(m, Edge(m.phi1(vit2)), edge_queue, edge_queue_info, edge_cost);
				vit2 = m.phi1(m.phi2(vit2));
			} while(vit2 != stop);
		}
		else
			update_edge_queue(m, Edge(vit), edge_queue, edge_queue_info, edge_cost);

		vit = m.phi2(m.phi_1(vit));
	} while(vit != e1.dart);
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
pre_collapse_edge_length(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	typename mesh_traits<MESH>::Edge& e1,
	typename mesh_traits<MESH>::Edge& e2,
	CellQueue<typename mesh_traits<MESH>::Edge>& edge_queue,
	typename mesh_traits<MESH>::template Attribute<typename CellQueue<typename mesh_traits<MESH>::Edge>::CellQueueInfo>* edge_queue_info
)
{
	return pre_collapse_edge_length(m.mesh(), e, e1, e2, edge_queue, edge_queue_info);
}

template <typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
post_collapse_edge_length(
	const MESH& m,
	typename mesh_traits<MESH>::Edge& e1,
	typename mesh_traits<MESH>::Edge& e2,
	CellQueue<typename mesh_traits<MESH>::Edge>& edge_queue,
	typename mesh_traits<MESH>::template Attribute<typename CellQueue<typename mesh_traits<MESH>::Edge>::CellQueueInfo>* edge_queue_info,
	const FUNC& edge_cost
)
{
	return post_collapse_edge_length(m.mesh(), e1, e2, edge_queue, edge_queue_info, edge_cost);
}

////////////////
// CellFilter //
////////////////

template <typename MESH, typename FUNC,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
void
post_collapse_edge_length(
	const CellFilter<MESH>& cf,
	typename mesh_traits<MESH>::Edge& e1,
	typename mesh_traits<MESH>::Edge& e2,
	CellQueue<typename mesh_traits<MESH>::Edge>& edge_queue,
	typename mesh_traits<MESH>::template Attribute<typename CellQueue<typename mesh_traits<MESH>::Edge>::CellQueueInfo>* edge_queue_info,
	const FUNC& edge_cost
)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using EdgeQueueInfo = typename CellQueue<Edge>::CellQueueInfo;

	const MESH& m = cf.mesh();

	Dart vit = e1.dart;
	do
	{
		Edge e(m.phi1(vit));
		if (cf.filter(e))
			update_edge_queue(m, e, edge_queue, edge_queue_info, edge_cost);
		if(vit == e1.dart || vit == e2.dart)
		{
			e = Edge(vit);
			if (cf.filter(e))
				update_edge_queue(m, e, edge_queue, edge_queue_info, edge_cost);

			Dart vit2 = m.template phi<121>(vit);
			Dart stop = m.phi2(vit);
			do
			{
				e = Edge(vit2);
				if (cf.filter(e))
					update_edge_queue(m, e, edge_queue, edge_queue_info, edge_cost);
				e = Edge(m.phi1(vit2));
				if (cf.filter(e))
					update_edge_queue(m, e, edge_queue, edge_queue_info, edge_cost);
				vit2 = m.phi1(m.phi2(vit2));
			} while(vit2 != stop);
		}
		else
		{
			e = Edge(vit);
			if (cf.filter(e))
				update_edge_queue(m, e, edge_queue, edge_queue_info, edge_cost);
		}

		vit = m.phi2(m.phi_1(vit));
	} while(vit != e1.dart);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_QUEUE_EDGE_LENGTH_H_
