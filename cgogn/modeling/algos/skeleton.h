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
#include <cgogn/core/functions/traversals/halfedge.h>
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

inline CMap2::Vertex opposite_vertex(CMap2& m, CMap2::HalfEdge he)
{
	return CMap2::Vertex(phi_1(m, he.dart));
}

inline void compute_halfedges_opposite_angle(CMap2& m, CMap2::Face f, CMap2::Attribute<Vec3>* vertex_position,
											 CMap2::Attribute<Scalar>* halfedge_opposite_angle)
{
	Scalar zero_threshold = 1e-5;

	Dart ha = f.dart;
	Dart hb = phi1(m, ha);
	Dart hc = phi1(m, hb);
	Scalar a = geometry::length(m, CMap2::Edge(ha), vertex_position);
	Scalar a2 = a * a;
	Scalar b = geometry::length(m, CMap2::Edge(hb), vertex_position);
	Scalar b2 = b * b;
	Scalar c = geometry::length(m, CMap2::Edge(hc), vertex_position);
	Scalar c2 = c * c;
	if (a < zero_threshold || b < zero_threshold || c < zero_threshold)
	{
		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(ha)) = -1;
		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hb)) = -1;
		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hc)) = -1;
	}
	else
	{
		/// Opposite angles (from law of cosines)
		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(ha)) =
			acos(std::clamp((-a2 + b2 + c2) / (2 * b * c), -1.0, 1.0));
		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hb)) =
			acos(std::clamp((+a2 - b2 + c2) / (2 * a * c), -1.0, 1.0));
		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hc)) =
			acos(std::clamp((+a2 + b2 - c2) / (2 * a * b), -1.0, 1.0));
	}
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void mean_curvature_skeleton(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							 typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_medial_point)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;

	// std::vector<Vec3> medial_points;
	// medial_points.reserve(nb_cells<Vertex>(m));
	// foreach_cell(m, [&](Vertex v) -> bool {
	// 	medial_points.push_back(value<Vec3>(m, vertex_medial_point, v));
	// 	return true;
	// });
	// acc::KDTree<3, uint32>* medial_points_kdt = new acc::KDTree<3, uint32>(medial_points);

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

	Scalar wL = 1.0;
	Scalar wH = 0.1;
	Scalar wM = 0.2;
	Scalar edge_collapse_threshold = 0.004 * bb_diag;
	Scalar zero_threshold = 1e-5;

	auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");
	auto vertex_wL = add_attribute<Scalar, Vertex>(m, "__vertex_wL");
	auto vertex_wH = add_attribute<Scalar, Vertex>(m, "__vertex_wH");
	auto vertex_wM = add_attribute<Scalar, Vertex>(m, "__vertex_wM");
	auto halfedge_opposite_angle = add_attribute<Scalar, HalfEdge>(m, "__halfedge_opposite_angle");
	auto edge_weight = add_attribute<Scalar, Edge>(m, "__edge_weight");
	auto vertex_is_fixed = add_attribute<bool, Vertex>(m, "__vertex_is_fixed");
	auto vertex_is_split = add_attribute<bool, Vertex>(m, "__vertex_is_split");

	vertex_wL->fill(wL);
	vertex_wH->fill(wH);
	vertex_wM->fill(wM);
	vertex_is_fixed->fill(false);
	vertex_is_split->fill(false);

	CellCache<MESH> cache(m);

	for (uint32 i = 0; i < 2; ++i)
	{
		std::cout << "index vertices" << std::endl;
		uint32 nb_vertices = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			value<uint32>(m, vertex_index, v) = nb_vertices++;
			return true;
		});

		std::cout << "compute laplacian edge weights" << std::endl;
		foreach_cell(m, [&](Edge e) -> bool {
			std::vector<Scalar> angles = geometry::opposite_angles(m, e, vertex_position);
			Scalar weight = 0.0;
			for (Scalar a : angles)
				weight += 1.0 / std::tan(std::clamp(a, Scalar(-0.99), Scalar(0.99)));
			// weight += std::tan(M_PI_2 - a);
			if (weight < 0.0)
				weight = 0.0;
			value<Scalar>(m, edge_weight, e) = weight;
			return true;
		});

		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_vertices * 10);

		// smooth
		std::cout << "setup matrix laplacian coeffs" << std::endl;
		foreach_cell(m, [&](Edge e) -> bool {
			Scalar w = value<Scalar>(m, edge_weight, e);
			auto iv = incident_vertices(m, e);
			uint32 vidx1 = value<uint32>(m, vertex_index, iv[0]);
			uint32 vidx2 = value<uint32>(m, vertex_index, iv[1]);
			Scalar v1_wL = value<Scalar>(m, vertex_wL, iv[0]);
			Scalar v2_wL = value<Scalar>(m, vertex_wL, iv[1]);
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w* v1_wL));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w* v2_wL));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w* v1_wL));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w* v2_wL));
			return true;
		});

		// velocity
		std::cout << "setup matrix velocity coeffs" << std::endl;
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			Acoeffs.push_back(
				Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), value<Scalar>(m, vertex_wH, v)));
			return true;
		});

		// medial
		std::cout << "setup matrix medial coeffs" << std::endl;
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			Acoeffs.push_back(
				Eigen::Triplet<Scalar>(int(2 * nb_vertices + vidx), int(vidx), value<Scalar>(m, vertex_wM, v)));
			return true;
		});

		std::cout << "setup matrix" << std::endl;
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(3 * nb_vertices, nb_vertices);
		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(3 * nb_vertices, 3);

		std::cout << "setup right hand side" << std::endl;
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			b(vidx, 0) = 0;
			b(vidx, 1) = 0;
			b(vidx, 2) = 0;
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			Scalar v_wH = value<Scalar>(m, vertex_wH, v);
			b(nb_vertices + vidx, 0) = v_wH * pos[0];
			b(nb_vertices + vidx, 1) = v_wH * pos[1];
			b(nb_vertices + vidx, 2) = v_wH * pos[2];
			const Vec3& medp = value<Vec3>(m, vertex_medial_point, v);
			Scalar v_wM = value<Scalar>(m, vertex_wM, v);
			b(2 * nb_vertices + vidx, 0) = v_wM * medp[0];
			b(2 * nb_vertices + vidx, 1) = v_wM * medp[1];
			b(2 * nb_vertices + vidx, 2) = v_wM * medp[2];
			x(vidx, 0) = pos[0];
			x(vidx, 1) = pos[1];
			x(vidx, 2) = pos[2];
			return true;
		});

		std::cout << "solve least squares problem" << std::endl;
		// Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);
		// x = solver.solveWithGuess(b, x);
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At = A.transpose();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(At * A);
		x = solver.solve(At * b);

		std::cout << "store solution" << std::endl;
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			Vec3& pos = value<Vec3>(m, vertex_position, v);
			pos[0] = x(vidx, 0);
			pos[1] = x(vidx, 1);
			pos[2] = x(vidx, 2);
			return true;
		});

		std::cout << "update vertices constraints" << std::endl;
		foreach_cell(m, [&](Vertex v) -> bool {
			if (value<bool>(m, vertex_is_fixed, v))
			{
				value<Scalar>(m, vertex_wL, v) = 0.0;
				value<Scalar>(m, vertex_wH, v) = 1.0 / zero_threshold;
				value<Scalar>(m, vertex_wM, v) = 0.0;
				return true;
			}
			value<Scalar>(m, vertex_wL, v) = wL;
			value<Scalar>(m, vertex_wH, v) = wH;
			value<Scalar>(m, vertex_wM, v) = wM;
			if (value<bool>(m, vertex_is_split, v))
			{
				value<Scalar>(m, vertex_wL, v) = wL;
				value<Scalar>(m, vertex_wH, v) = wH;
				value<Scalar>(m, vertex_wM, v) = 0.0;
			}
			return true;
		});

		std::cout << "collapse short edges" << std::endl;
		uint32 nb_collapsed_edges = 0;
		cache.template build<Edge>();
		bool has_short_edge = false;
		do
		{
			has_short_edge = false;
			foreach_cell(cache, [&](Edge e) -> bool {
				if (geometry::length(m, e, vertex_position) < edge_collapse_threshold)
				{
					if (edge_can_collapse(m, e))
					{
						has_short_edge = true;
						std::vector<Vertex> iv = incident_vertices(m, e);
						Vec3 newp =
							(value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
						const Vec3& mp0 = value<Vec3>(m, vertex_medial_point, iv[0]);
						const Vec3& mp1 = value<Vec3>(m, vertex_medial_point, iv[1]);
						Scalar d0 = (mp0 - newp).squaredNorm();
						Scalar d1 = (mp1 - newp).squaredNorm();
						Vec3 newmp = d0 < d1 ? mp0 : mp1;
						Vertex cv = collapse_edge(m, e);
						value<Vec3>(m, vertex_position, cv) = newp;
						value<Vec3>(m, vertex_medial_point, cv) = newmp;
						// std::pair<uint32, Scalar> k_res;
						// medial_points_kdt->find_nn(mp, &k_res);
						// value<Vec3>(m, vertex_medial_point, cv) = medial_points_kdt->vertex(k_res.first);
						++nb_collapsed_edges;
					}
				}
				return true;
			});
			cache.template build<Edge>();
		} while (has_short_edge);
		std::cout << "  nb collapsed edges: " << nb_collapsed_edges << std::endl;

		std::cout << "compute angles" << std::endl;
		foreach_cell(m, [&](Face f) -> bool {
			compute_halfedges_opposite_angle(m, f, vertex_position, halfedge_opposite_angle.get());
			return true;
		});

		std::cout << "cut long edges" << std::endl;
		uint32 nb_cut_edges = 0;
		cache.template build<Edge>();
		bool has_long_edge = false;
		do
		{
			has_long_edge = false;
			foreach_cell(cache, [&](Edge e) -> bool {
				// std::vector<Scalar> op_angles = geometry::opposite_angles(m, e, vertex_position);
				std::vector<HalfEdge> ihe = incident_halfedges(m, e);
				Scalar alpha0 = value<Scalar>(m, halfedge_opposite_angle, ihe[0]);
				Scalar alpha1 = value<Scalar>(m, halfedge_opposite_angle, ihe[1]);
				Scalar cut_threshold = 110.0 * M_PI / 180.0;

				if (alpha0 < cut_threshold || alpha1 < cut_threshold)
					return true;

				has_long_edge = true;
				std::vector<Vertex> iv = incident_vertices(m, e);
				const Vec3& p0 = value<Vec3>(m, vertex_position, iv[0]);
				const Vec3& p1 = value<Vec3>(m, vertex_position, iv[1]);

				uint32 largest = alpha0 > alpha1 ? 0 : 1;
				Vertex op_vertex = opposite_vertex(m, ihe[largest]);
				Vec3 vec = (p1 - p0).normalized();
				Vec3 proj = value<Vec3>(m, vertex_position, op_vertex) - p0;
				Scalar t = vec.dot(proj);
				Vec3 newp = p0 + t * vec;

				const Vec3& mp0 = value<Vec3>(m, vertex_medial_point, iv[0]);
				const Vec3& mp1 = value<Vec3>(m, vertex_medial_point, iv[1]);
				vec = (mp1 - mp0).normalized();
				Vec3 newmp = mp0 + t * vec;

				Vertex cv = cut_edge(m, e);
				cut_incident_faces(m, cv);
				value<Vec3>(m, vertex_position, cv) = newp;
				value<Vec3>(m, vertex_medial_point, cv) = newmp;
				// std::pair<uint32, Scalar> k_res;
				// medial_points_kdt->find_nn(mp, &k_res);
				// value<Vec3>(m, vertex_medial_point, cv) = medial_points_kdt->vertex(k_res.first);
				value<bool>(m, vertex_is_split, cv) = true;
				value<bool>(m, vertex_is_fixed, cv) = false;

				++nb_cut_edges;
				return true;
			});
			cache.template build<Edge>();
		} while (has_long_edge);
		std::cout << "  nb cut edges: " << nb_cut_edges << std::endl;

		std::cout << "detect & mark degeneracies" << std::endl;
		uint32 nb_fixed_vertices = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			if (value<bool>(m, vertex_is_fixed, v))
				return true;
			uint32 count = 0;
			foreach_incident_edge(m, v, [&](Edge ie) -> bool {
				Scalar l = geometry::length(m, ie, vertex_position);
				if (l < edge_collapse_threshold / 10.0 && !edge_can_collapse(m, ie))
					++count;
				return true;
			});
			bool is_fixed = count >= 2;
			if (is_fixed)
			{
				value<bool>(m, vertex_is_fixed, v) = true;
				++nb_fixed_vertices;
			}
			else
				value<bool>(m, vertex_is_fixed, v) = false;
			return true;
		});
		std::cout << "  nb fixed vertices: " << nb_fixed_vertices << std::endl;
	}

	remove_attribute<Vertex>(m, vertex_index);
	remove_attribute<Vertex>(m, vertex_wL);
	remove_attribute<Vertex>(m, vertex_wH);
	remove_attribute<Vertex>(m, vertex_wM);
	remove_attribute<Edge>(m, edge_weight);
	remove_attribute<HalfEdge>(m, halfedge_opposite_angle);
	remove_attribute<Vertex>(m, vertex_is_fixed);
	remove_attribute<Vertex>(m, vertex_is_split);

	// delete medial_points_kdt;
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SKELETON_H_
