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

#ifndef CGOGN_MODELING_ALGOS_SKELETON_H_
#define CGOGN_MODELING_ALGOS_SKELETON_H_

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/medial_axis.h>
#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>

#include <Eigen/Sparse>

#include <array>
#include <vector>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

template <typename MESH>
void mean_curvature_skeleton(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							 typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	Scalar wH = 1.0;
	Scalar wL = 1.0;

	// modeling::pliant_remeshing(m, vertex_position, 0.9);
	// geometry::compute_normal(m, vertex_position, vertex_normal);

	// auto sbc = add_attribute<Vec3, Vertex>(m, "__shrinking_ball_centers");
	// modeling::shrinking_ball_centers(m, vertex_position, vertex_normal, sbc.get());

	auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");
	uint32 nb_vertices = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, vertex_index, v) = nb_vertices++;
		return true;
	});

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * nb_vertices, nb_vertices);

	auto edge_weight = add_attribute<Scalar, Edge>(m, "__edge_weight");
	parallel_foreach_cell(m, [&](Edge e) -> bool {
		std::vector<Scalar> angles = geometry::opposite_angles(m, e, vertex_position);
		Scalar& weight = value<Scalar>(m, edge_weight, e);
		for (Scalar a : angles)
			weight += std::tan(M_PI_2 - a);
		weight /= uint32(angles.size());
		return true;
	});

	std::vector<Eigen::Triplet<Scalar>> Acoeffs;
	Acoeffs.reserve(nb_vertices * 10);

	foreach_cell(m, [&](Edge e) -> bool {
		Scalar w = wL * value<Scalar>(m, edge_weight, e);
		auto vertices = incident_vertices(m, e);
		uint32 vidx1 = value<uint32>(m, vertex_index, vertices[0]);
		uint32 vidx2 = value<uint32>(m, vertex_index, vertices[1]);
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w));
		return true;
	});

	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), wH));
		return true;
	});

	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);

	Eigen::MatrixXd x(nb_vertices, 3);
	Eigen::MatrixXd b(2 * nb_vertices, 3);

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		b(vidx, 0) = 0;
		b(vidx, 1) = 0;
		b(vidx, 2) = 0;
		const Vec3& pos = value<Vec3>(m, vertex_position, v);
		b(nb_vertices + vidx, 0) = wH * pos[0];
		b(nb_vertices + vidx, 1) = wH * pos[1];
		b(nb_vertices + vidx, 2) = wH * pos[2];
		x(vidx, 0) = pos[0];
		x(vidx, 1) = pos[1];
		x(vidx, 2) = pos[2];
		return true;
	});

	x = solver.solveWithGuess(b, x);

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		Vec3& pos = value<Vec3>(m, vertex_position, v);
		pos[0] = x(vidx, 0);
		pos[1] = x(vidx, 1);
		pos[2] = x(vidx, 2);
		return true;
	});

	geometry::compute_normal(m, vertex_position, vertex_normal);

	// remove_attribute<Vertex>(m, sbc);
	remove_attribute<Vertex>(m, vertex_index);
	remove_attribute<Edge>(m, edge_weight);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SKELETON_H_
