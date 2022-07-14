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

#ifndef CGOGN_MODELING_ALGOS_REMESHING_BLOSSOM_QUAD_REMESHING_H_
#define CGOGN_MODELING_ALGOS_REMESHING_BLOSSOM_QUAD_REMESHING_H_

#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

#include <libacc/bvh_tree.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

///////////
// CMap2 //
///////////

inline std::vector<CMap2::Vertex> edge_quad_vertices(CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Vertex> eqv(4);
	eqv[0] = CMap2::Vertex(phi1(m, e.dart));
	eqv[1] = CMap2::Vertex(phi_1(m, e.dart));
	eqv[2] = CMap2::Vertex(e.dart);
	eqv[3] = CMap2::Vertex(phi<2, -1>(m, e.dart));
	return eqv;
}

inline std::vector<CMap2::Edge> butterfly_neighbour_edges(CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Edge> edges;

	Dart d = e.dart;
	if (!is_boundary(m, d))
	{
		Dart d12 = phi<1, 2>(m, d);
		if (!is_boundary(m, d12))
		{
			Dart d1212 = phi<1, 2>(m, d12);
			if (!is_boundary(m, d1212))
				edges.push_back(CMap2::Edge(d1212));
			Dart d12_12 = phi<-1, 2>(m, d12);
			if (!is_boundary(m, d12_12))
				edges.push_back(CMap2::Edge(d12_12));
		}
		Dart d_12 = phi<-1, 2>(m, d);
		if (!is_boundary(m, d_12))
		{
			Dart d_1212 = phi<1, 2>(m, d_12);
			if (!is_boundary(m, d_1212))
				edges.push_back(CMap2::Edge(d_1212));
			Dart d_12_12 = phi<-1, 2>(m, d_12);
			if (!is_boundary(m, d_12_12))
				edges.push_back(CMap2::Edge(d_12_12));
		}
	}
	Dart d2 = phi2(m, e.dart);
	if (!is_boundary(m, d2))
	{
		Dart d212 = phi<1, 2>(m, d2);
		if (!is_boundary(m, d212))
		{
			Dart d21212 = phi<1, 2>(m, d212);
			if (!is_boundary(m, d21212))
				edges.push_back(CMap2::Edge(d21212));
			Dart d212_12 = phi<-1, 2>(m, d212);
			if (!is_boundary(m, d212_12))
				edges.push_back(CMap2::Edge(d212_12));
		}
		Dart d2_12 = phi<-1, 2>(m, d2);
		if (!is_boundary(m, d2_12))
		{
			Dart d2_1212 = phi<1, 2>(m, d2_12);
			if (!is_boundary(m, d2_1212))
				edges.push_back(CMap2::Edge(d2_1212));
			Dart d2_12_12 = phi<-1, 2>(m, d2_12);
			if (!is_boundary(m, d2_12_12))
				edges.push_back(CMap2::Edge(d2_12_12));
		}
	}

	return edges;
}

/////////////
// GENERIC //
/////////////

inline Scalar quad_quality(const std::vector<Vec3>& points)
{
	const auto ps = points.size();

	std::vector<Vec3> edges(ps, {0, 0, 0});
	std::vector<Scalar> edgeSqrNorm(ps, -1);

	for (uint32 i = 0; i < ps; ++i)
	{
		uint32 j = (i + 1) % ps;
		edges[i] = points[j] - points[i];
		edgeSqrNorm[i] = edges[i].squaredNorm();
	}

	auto determinant = [](const Vec3& p1, const Vec3& p2) { return p1[0] * p2[1] - p1[1] * p2[0]; };

	std::vector<Scalar> detJ(ps, -1);
	for (uint32 i = 0; i < ps; ++i)
	{
		uint32 il = (i - 1 + ps) % ps;
		detJ[i] = determinant(edges[i], -edges[il]);
	}

	std::vector<Scalar> c(ps, -1);

	for (uint32 i = 0; i < ps; ++i)
	{
		uint32 il = (i - 1 + ps) % ps;
		c[i] = detJ[i] / (edgeSqrNorm[i] + edgeSqrNorm[il]);
	}

	auto detMin = std::min_element(detJ.begin(), detJ.end());

	if (*detMin > 0)
	{
		auto cMin = std::min_element(c.begin(), c.end());
		return 2 * *cMin;
	}
	else
		return *detMin;
}

template <typename MESH>
void blossom_quad_remeshing(MESH& m,
							std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	struct prio_edge
	{
		Scalar w;
		Edge e;
		int rank = -1;
	};

	uint32 nbe = nb_cells<Edge>(m);

	std::vector<prio_edge> edges;
	edges.reserve(nbe);

	foreach_cell(m, [&](Edge e) -> bool {
		if (is_incident_to_boundary(m, e))
			return true;

		std::vector<Face> ifaces = incident_faces(m, e);
		if (codegree(m, ifaces[0]) > 3 || codegree(m, ifaces[1]) > 3)
			return true;

		std::vector<Vertex> eqv = edge_quad_vertices(m, e);
		std::vector<Vec3> eqv_pos;
		eqv_pos.reserve(4);
		for (Vertex v : eqv)
			eqv_pos.push_back(value<Vec3>(m, vertex_position, v));
		Scalar weight = quad_quality(eqv_pos);

		edges.push_back({weight, e});

		return true;
	});

	auto deleted_faces = add_attribute<bool, Face>(m, "__deleted_faces");
	deleted_faces->fill(false);

	std::vector<Edge> edges_to_delete;
	edges_to_delete.reserve(nbe / 2);

	std::sort(edges.begin(), edges.end(), [](const prio_edge& a, const prio_edge& b) { return a.w > b.w; });

	// auto edge_rank = add_attribute<int, Edge>(m, "__edge_rank");
	// edge_rank->fill(std::numeric_limits<int>::max());
	// for (uint32 i = 0; i < edges.size(); ++i)
	// 	value<int>(m, edge_rank, edges[i].e) = i;

	for (const auto& [w, e, rank] : edges)
	{
		// std::queue<Edge> queue;
		// queue.push(e);
		// while (!queue.empty())
		// {
		// Edge e = queue.front();
		// queue.pop();

		// int rank = value<int>(m, edge_rank, e);
		// if (rank == std::numeric_limits<int>::max())
		// 	continue;

		std::vector<Face> ifaces = incident_faces(m, e);

		if (value<bool>(m, deleted_faces, ifaces[0]) || value<bool>(m, deleted_faces, ifaces[1]))
			continue;
		if (w < 0.1)
			continue;

		value<bool>(m, deleted_faces, ifaces[0]) = true;
		value<bool>(m, deleted_faces, ifaces[1]) = true;
		edges_to_delete.push_back(e);

		// 	std::vector<Edge> next_edges = butterfly_neighbour_edges(m, e);

		// 	std::sort(next_edges.begin(), next_edges.end(),
		// 			  [&](Edge a, Edge b) { return value<int>(m, edge_rank, a) < value<int>(m, edge_rank, b); });

		// 	for (Edge ne : next_edges)
		// 		queue.push(ne);
		// }
	}

	remove_attribute<Face>(m, deleted_faces);
	// remove_attribute<Edge>(m, edge_rank);

	for (Edge e : edges_to_delete)
		merge_incident_faces(m, e);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_REMESHING_BLOSSOM_QUAD_REMESHING_H_
