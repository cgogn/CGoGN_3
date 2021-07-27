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

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/medial_axis.h>
#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>

#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <vector>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

///////////
// CMap2 //
///////////

inline void cut_incident_faces(CMap2& m, CMap2::Vertex v)
{
	std::vector<CMap2::Face> ifaces = incident_faces(m, v);
	for (CMap2::Face f : ifaces)
		cut_face(m, CMap2::Vertex(f.dart), CMap2::Vertex(phi<11>(m, f.dart)));
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void mean_curvature_skeleton(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							 typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
							 typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_medial_point)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	std::vector<Vec3> medial_points;
	medial_points.reserve(nb_cells<Vertex>(m));
	foreach_cell(m, [&](Vertex v) -> bool {
		medial_points.push_back(value<Vec3>(m, vertex_medial_point, v));
		return true;
	});
	acc::KDTree<3, uint32>* medial_points_kdt = new acc::KDTree<3, uint32>(medial_points);

	Vec3 bb_min, bb_max;
	for (uint32 i = 0; i < 3; ++i)
	{
		bb_min[i] = std::numeric_limits<float64>::max();
		bb_max[i] = std::numeric_limits<float64>::lowest();
	}
	for (const Vec3& p : *vertex_position)
	{
		for (uint32 i = 0; i < 3; ++i)
		{
			if (p[i] < bb_min[i])
				bb_min[i] = p[i];
			if (p[i] > bb_max[i])
				bb_max[i] = p[i];
		}
	}
	Scalar bb_diag = (bb_max - bb_min).norm();
	Scalar edge_collapse_threshold = 0.01 * bb_diag;

	for (uint32 i = 0; i < 15; ++i)
	{
		Scalar wH = 20.0;
		Scalar wL = 10.0;
		Scalar wM = 1.0;

		auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");
		uint32 nb_vertices = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			value<uint32>(m, vertex_index, v) = nb_vertices++;
			return true;
		});

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(3 * nb_vertices, nb_vertices);

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

		// smooth
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

		// velocity
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), wH));
			return true;
		});

		// medial
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(2 * nb_vertices + vidx), int(vidx), wM));
			return true;
		});

		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(3 * nb_vertices, 3);

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			b(vidx, 0) = 0;
			b(vidx, 1) = 0;
			b(vidx, 2) = 0;
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			b(nb_vertices + vidx, 0) = wH * pos[0];
			b(nb_vertices + vidx, 1) = wH * pos[1];
			b(nb_vertices + vidx, 2) = wH * pos[2];
			const Vec3& medp = value<Vec3>(m, vertex_medial_point, v);
			b(2 * nb_vertices + vidx, 0) = wH * medp[0];
			b(2 * nb_vertices + vidx, 1) = wH * medp[1];
			b(2 * nb_vertices + vidx, 2) = wH * medp[2];
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

		uint32 nb_cut_edges = 0;
		foreach_cell(m, [&](Edge e) -> bool {
			std::vector<Scalar> angles = geometry::opposite_angles(m, e, vertex_position);
			bool should_cut_edge =
				std::any_of(angles.begin(), angles.end(), [](Scalar s) { return s > 110 * M_PI / 180.0; });
			if (should_cut_edge)
			{
				std::vector<Vertex> iv = incident_vertices(m, e);
				Vec3 mp = (value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
				Vertex cv = cut_edge(m, e);
				cut_incident_faces(m, cv);
				value<Vec3>(m, vertex_position, cv) = mp;
				std::pair<uint32, Scalar> k_res;
				medial_points_kdt->find_nn(mp, &k_res);
				value<Vec3>(m, vertex_medial_point, cv) = medial_points_kdt->vertex(k_res.first);
				++nb_cut_edges;
			}
		});

		uint32 nb_collapsed_edges = 0;
		foreach_cell(m, [&](Edge e) -> bool {
			Scalar l = geometry::length(m, e, vertex_position);
			if (l < edge_collapse_threshold)
			{
				if (edge_can_collapse(m, e))
				{
					std::vector<Vertex> iv = incident_vertices(m, e);
					Vec3 mp = (value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
					Vertex cv = collapse_edge(m, e);
					value<Vec3>(m, vertex_position, cv) = mp;
					std::pair<uint32, Scalar> k_res;
					medial_points_kdt->find_nn(mp, &k_res);
					value<Vec3>(m, vertex_medial_point, cv) = medial_points_kdt->vertex(k_res.first);
					++nb_collapsed_edges;
				}
			}
		});

		std::cout << "iteration " << i << ": " << std::endl;
		std::cout << "  nb cut edges " << nb_cut_edges << ": " << std::endl;
		std::cout << "  nb collapsed edges " << nb_collapsed_edges << ": " << std::endl;

		geometry::compute_normal(m, vertex_position, vertex_normal);

		remove_attribute<Vertex>(m, vertex_index);
		remove_attribute<Edge>(m, edge_weight);
	}
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SKELETON_H_
