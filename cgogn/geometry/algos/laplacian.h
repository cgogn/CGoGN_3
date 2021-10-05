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

#ifndef CGOGN_GEOMETRY_ALGOS_LAPLACIAN_H_
#define CGOGN_GEOMETRY_ALGOS_LAPLACIAN_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>

// #include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Sparse>

namespace cgogn
{

namespace geometry
{

///////////
// CMap2 //
///////////

Scalar edge_cotan_weight(const CMap2& m, CMap2::Edge e, const CMap2::Attribute<Vec3>* vertex_position)
{
	Scalar result = 0.0;

	Dart d1 = e.dart;
	Dart d2 = phi2(m, d1);

	const Vec3& p1 = value<Vec3>(m, vertex_position, CMap2::Vertex(d1));
	const Vec3& p2 = value<Vec3>(m, vertex_position, CMap2::Vertex(d2));

	const Vec3& p3 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi_1(m, d1)));
	Vec3 vecR = p1 - p3;
	Vec3 vecL = p2 - p3;
	Scalar e1value = vecR.dot(vecL) / vecR.cross(vecL).norm();

	result += e1value / 2.0;

	if (!is_boundary(m, d2))
	{
		const Vec3& p4 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi_1(m, d2)));
		Vec3 vecR = p2 - p4;
		Vec3 vecL = p1 - p4;
		Scalar e2value = vecR.dot(vecL) / vecR.cross(vecL).norm();

		result += e2value / 2.0;
	}

	return std::clamp(result, 0.0, 1.0);
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void compute_edge_cotan_weight(const MESH& m,
							   const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							   typename mesh_traits<MESH>::template Attribute<Scalar>* edge_weight)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	parallel_foreach_cell(m, [&](Edge e) -> bool {
		value<Scalar>(m, edge_weight, e) = edge_cotan_weight(m, e, vertex_position);
		return true;
	});
}

template <typename MESH>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> cotan_laplacian_matrix(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<uint32>* vertex_index,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Scalar>* edge_weight)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	// setup matrix
	uint32 nb_vertices = nb_cells<Vertex>(m);
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL(nb_vertices, nb_vertices);
	std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
	LAPLcoeffs.reserve(nb_vertices * 10);
	foreach_cell(m, [&](Edge e) -> bool {
		Scalar w = value<Scalar>(m, edge_weight, e);
		auto vertices = incident_vertices(m, e);
		uint32 vidx1 = value<uint32>(m, vertex_index, vertices[0]);
		uint32 vidx2 = value<uint32>(m, vertex_index, vertices[1]);
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w));
		return true;
	});
	LAPL.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());

	return LAPL;
}

template <typename MESH>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> cotan_laplacian_matrix(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<uint32>* vertex_index,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	auto edge_weight = add_attribute<Scalar, Edge>(m, "__edge_weight");
	compute_edge_cotan_weight(m, vertex_position, edge_weight.get());
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL =
		cotan_laplacian_matrix(m, vertex_index, vertex_position, edge_weight.get());
	remove_attribute<Edge>(m, edge_weight);
	return LAPL;
}

template <typename MESH>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> topo_laplacian_matrix(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<uint32>* vertex_index)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	uint32 nb_vertices = nb_cells<Vertex>(m);
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL(nb_vertices, nb_vertices);
	std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
	LAPLcoeffs.reserve(nb_vertices * 10);
	foreach_cell(m, [&](Edge e) -> bool {
		auto vertices = incident_vertices(m, e);
		uint32 vidx1 = value<uint32>(m, vertex_index, vertices[0]);
		uint32 vidx2 = value<uint32>(m, vertex_index, vertices[1]);
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), 1));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), 1));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -1));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -1));
		return true;
	});
	LAPL.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());

	return LAPL;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_LAPLACIAN_H_
