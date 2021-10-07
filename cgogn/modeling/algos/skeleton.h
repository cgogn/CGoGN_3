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

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/vector_traits.h>

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

inline std::vector<CMap2::Vertex> opposite_vertices(CMap2& m, CMap2::Edge e)
{
	return {CMap2::Vertex(phi_1(m, e.dart)), CMap2::Vertex(phi_1(m, phi2(m, e.dart)))};
}

// inline void compute_halfedges_opposite_angle(CMap2& m, CMap2::Face f, CMap2::Attribute<Vec3>* vertex_position,
// 											 CMap2::Attribute<Scalar>* halfedge_opposite_angle)
// {
// 	Scalar zero_threshold = 1e-5;

// 	Dart ha = f.dart;
// 	Dart hb = phi1(m, ha);
// 	Dart hc = phi1(m, hb);
// 	Scalar a = geometry::length(m, CMap2::Edge(ha), vertex_position);
// 	Scalar a2 = a * a;
// 	Scalar b = geometry::length(m, CMap2::Edge(hb), vertex_position);
// 	Scalar b2 = b * b;
// 	Scalar c = geometry::length(m, CMap2::Edge(hc), vertex_position);
// 	Scalar c2 = c * c;
// 	if (a < zero_threshold || b < zero_threshold || c < zero_threshold)
// 	{
// 		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(ha)) = -1;
// 		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hb)) = -1;
// 		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hc)) = -1;
// 	}
// 	else
// 	{
// 		/// Opposite angles (from law of cosines)
// 		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(ha)) =
// 			acos(std::clamp((-a2 + b2 + c2) / (2 * b * c), -1.0, 1.0));
// 		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hb)) =
// 			acos(std::clamp((+a2 - b2 + c2) / (2 * a * c), -1.0, 1.0));
// 		value<Scalar>(m, halfedge_opposite_angle, CMap2::HalfEdge(hc)) =
// 			acos(std::clamp((+a2 + b2 - c2) / (2 * a * b), -1.0, 1.0));
// 	}
// }

/////////////
// GENERIC //
/////////////

template <typename MESH>
struct MeanCurvatureSkeleton_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	MeanCurvatureSkeleton_Helper(MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position,
								 Scalar surface_resampling_ratio)
		: m_(m), vertex_position_(vertex_position)
	{
		modeling::pliant_remeshing(m_, vertex_position_, surface_resampling_ratio, false, true);

		vertex_normal_ = add_attribute<Vec3, Vertex>(m_, "__vertex_normal");
		geometry::compute_normal(m_, vertex_position_.get(), vertex_normal_.get());

		vertex_medial_point_ = add_attribute<Vec3, Vertex>(m_, "__vertex_medial_point");
		geometry::shrinking_ball_centers(m_, vertex_position_.get(), vertex_normal_.get(), vertex_medial_point_.get());

		vertex_is_fixed_ = add_attribute<bool, Vertex>(m_, "__vertex_is_fixed");
		vertex_is_fixed_color_ = add_attribute<Vec3, Vertex>(m_, "__vertex_is_fixed_color");
		vertex_is_fixed_->fill(false);
		vertex_is_fixed_color_->fill({1, 1, 1});

		auto [bb_min, bb_max] = geometry::bounding_box(*vertex_position);
		Scalar bb_diag = (bb_max - bb_min).norm();

		edge_collapse_threshold_ = 0.004 * bb_diag;

		vertex_index_ = add_attribute<uint32, Vertex>(m_, "__vertex_index");
		edge_weight_ = add_attribute<Scalar, Edge>(m_, "__edge_weight");
	}
	~MeanCurvatureSkeleton_Helper()
	{
		remove_attribute<Vertex>(m_, vertex_normal_);
		remove_attribute<Vertex>(m_, vertex_medial_point_);
		remove_attribute<Vertex>(m_, vertex_is_fixed_);
		remove_attribute<Vertex>(m_, vertex_index_);
		remove_attribute<Edge>(m_, edge_weight_);
	}

	MESH& m_;

	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	std::shared_ptr<Attribute<Vec3>> vertex_normal_;
	std::shared_ptr<Attribute<Vec3>> vertex_medial_point_;
	std::shared_ptr<Attribute<bool>> vertex_is_fixed_;
	std::shared_ptr<Attribute<Vec3>> vertex_is_fixed_color_;
	std::shared_ptr<Attribute<uint32>> vertex_index_;
	std::shared_ptr<Attribute<Scalar>> edge_weight_;

	Scalar wL_, wH_, wM_, edge_collapse_threshold_;
};

template <typename MESH>
void mean_curvature_skeleton(MESH& m,
							 std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& vertex_position,
							 Scalar wL, Scalar wH, Scalar wM, Scalar surface_resampling_ratio = 0.9)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;

	// static map to store helpers associated to meshes
	// allows to store context without polluting outer context and function api
	static std::unordered_map<MESH*, MeanCurvatureSkeleton_Helper<MESH>> helpers_;
	auto [it, inserted] = helpers_.try_emplace(&m, m, vertex_position, surface_resampling_ratio);
	MeanCurvatureSkeleton_Helper<MESH>& helper = it->second;

	uint32 nb_vertices = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, helper.vertex_index_, v) = nb_vertices++;
		return true;
	});

	geometry::compute_edge_cotan_weight(m, helper.vertex_position_.get(), helper.edge_weight_.get());

	std::vector<Eigen::Triplet<Scalar>> Acoeffs;
	Acoeffs.reserve(nb_vertices * 10);

	// smoothness
	foreach_cell(m, [&](Edge e) -> bool {
		Scalar w = value<Scalar>(m, helper.edge_weight_, e);
		auto iv = incident_vertices(m, e);
		uint32 vidx1 = value<uint32>(m, helper.vertex_index_, iv[0]);
		uint32 vidx2 = value<uint32>(m, helper.vertex_index_, iv[1]);
		Scalar v_wL1 = value<bool>(m, helper.vertex_is_fixed_, iv[0]) ? 0.0 : wL;
		Scalar v_wL2 = value<bool>(m, helper.vertex_is_fixed_, iv[1]) ? 0.0 : wL;
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w* v_wL1));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w* v_wL1));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w* v_wL2));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w* v_wL2));
		return true;
	});

	// velocity
	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, helper.vertex_index_, v);
		Scalar v_wH = value<bool>(m, helper.vertex_is_fixed_, v) ? 1e6 : wH;
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), v_wH));
		return true;
	});

	// medial attraction
	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, helper.vertex_index_, v);
		Scalar v_wM = value<bool>(m, helper.vertex_is_fixed_, v) ? 0.0 : wM;
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(2 * nb_vertices + vidx), int(vidx), v_wM));
		return true;
	});

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(3 * nb_vertices, nb_vertices);
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	Eigen::MatrixXd x(nb_vertices, 3);
	Eigen::MatrixXd b(3 * nb_vertices, 3);

	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, helper.vertex_index_, v);
		b(vidx, 0) = 0;
		b(vidx, 1) = 0;
		b(vidx, 2) = 0;
		const Vec3& pos = value<Vec3>(m, helper.vertex_position_, v);
		Scalar v_wH = value<bool>(m, helper.vertex_is_fixed_, v) ? 1e6 : wH;
		b(nb_vertices + vidx, 0) = v_wH * pos[0];
		b(nb_vertices + vidx, 1) = v_wH * pos[1];
		b(nb_vertices + vidx, 2) = v_wH * pos[2];
		const Vec3& medp = value<Vec3>(m, helper.vertex_medial_point_, v);
		Scalar v_wM = value<bool>(m, helper.vertex_is_fixed_, v) ? 0.0 : wM;
		b(2 * nb_vertices + vidx, 0) = v_wM * medp[0];
		b(2 * nb_vertices + vidx, 1) = v_wM * medp[1];
		b(2 * nb_vertices + vidx, 2) = v_wM * medp[2];
		return true;
	});

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At = A.transpose();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(At * A);
	x = solver.solve(At * b);

	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, helper.vertex_index_, v);
		Vec3& pos = value<Vec3>(m, helper.vertex_position_, v);
		pos[0] = x(vidx, 0);
		pos[1] = x(vidx, 1);
		pos[2] = x(vidx, 2);
		return true;
	});

	bool has_flat_edge = false;
	do
	{
		has_flat_edge = false;
		foreach_cell(m, [&](Edge e) -> bool {
			std::vector<Vertex> iv = incident_vertices(m, e);
			if (degree(m, iv[0]) < 5 || degree(m, iv[1]) < 5)
				return true;

			std::vector<Scalar> op_angles = geometry::opposite_angles(m, e, helper.vertex_position_.get());
			Scalar flip_threshold_low = 140.0 * M_PI / 180.0;
			if (op_angles[0] > flip_threshold_low && op_angles[1] > flip_threshold_low)
			{
				if (edge_can_flip(m, e))
				{
					if (flip_edge(m, e))
						has_flat_edge = true;
				}
			}
			return true;
		});
	} while (has_flat_edge);

	// uint32 nb_cut_edges = 0;
	// bool has_long_edge = false;
	// do
	// {
	// 	cache.template build<Edge>();
	// 	has_long_edge = false;
	// 	foreach_cell(cache, [&](Edge e) -> bool {
	// 		std::vector<Vertex> op_vertices = opposite_vertices(m, e);
	// 		std::vector<Scalar> op_angles = geometry::opposite_angles(m, e, helper.vertex_position_);
	// 		Scalar alpha0 = op_angles[0];
	// 		Scalar alpha1 = op_angles[1];
	// 		// std::vector<HalfEdge> ihe = incident_halfedges(m, e);
	// 		// Scalar alpha0 = value<Scalar>(m, halfedge_opposite_angle, ihe[0]);
	// 		// Scalar alpha1 = value<Scalar>(m, halfedge_opposite_angle, ihe[1]);
	// 		Scalar cut_threshold_low = 120.0 * M_PI / 180.0;
	// 		Scalar cut_threshold_high = 178.0 * M_PI / 180.0;

	// 		if (alpha0 < cut_threshold_low || alpha1 < cut_threshold_low || alpha0 > cut_threshold_high ||
	// 			alpha1 > cut_threshold_high)
	// 			return true;

	// 		has_long_edge = true;
	// 		std::vector<Vertex> iv = incident_vertices(m, e);
	// 		const Vec3& p0 = value<Vec3>(m, helper.vertex_position_, iv[0]);
	// 		const Vec3& p1 = value<Vec3>(m, helper.vertex_position_, iv[1]);

	// 		// uint32 largest = alpha0 > alpha1 ? 0 : 1;
	// 		// Vertex op_vertex = op_vertices[largest]; // opposite_vertex(m, ihe[largest]);
	// 		// Vec3 vec = (p1 - p0).normalized();
	// 		// Vec3 proj = value<Vec3>(m, helper.vertex_position_, op_vertex) - p0;
	// 		// Scalar t = vec.dot(proj);
	// 		// Vec3 newp = p0 + t * vec;
	// 		Vec3 newp = (p0 + p1) * 0.5;

	// 		const Vec3& mp0 = value<Vec3>(m, helper.vertex_medial_point_, iv[0]);
	// 		const Vec3& mp1 = value<Vec3>(m, helper.vertex_medial_point_, iv[1]);
	// 		// vec = (mp1 - mp0).normalized();
	// 		// Vec3 newmp = mp0 + t * vec;
	// 		Vec3 newmp = (mp0 + mp1) * 0.5;

	// 		Vertex cv = cut_edge(m, e);
	// 		cut_incident_faces(m, cv);
	// 		value<Vec3>(m, helper.vertex_position_, cv) = newp;
	// 		value<Vec3>(m, helper.vertex_medial_point_, cv) = newmp;
	// 		// std::pair<uint32, Scalar> k_res;
	// 		// medial_points_kdt->find_nn(mp, &k_res);
	// 		// value<Vec3>(m, helper.vertex_medial_point_, cv) = medial_points_kdt->vertex(k_res.first);

	// 		// value<bool>(m, vertex_is_split, cv) = true;
	// 		// value<bool>(m, vertex_is_fixed, cv) = false;

	// 		++nb_cut_edges;
	// 		return true;
	// 	});
	// } while (has_long_edge);

	bool has_short_edge = false;
	do
	{
		has_short_edge = false;
		foreach_cell(m, [&](Edge e) -> bool {
			std::vector<Vertex> iv = incident_vertices(m, e);
			if (geometry::length(m, e, helper.vertex_position_.get()) < helper.edge_collapse_threshold_)
			{
				if (edge_can_collapse(m, e))
				{
					has_short_edge = true;
					Vec3 newp = (value<Vec3>(m, helper.vertex_position_, iv[0]) +
								 value<Vec3>(m, helper.vertex_position_, iv[1])) *
								0.5;
					const Vec3& mp0 = value<Vec3>(m, helper.vertex_medial_point_, iv[0]);
					const Vec3& mp1 = value<Vec3>(m, helper.vertex_medial_point_, iv[1]);
					Scalar d0 = (mp0 - newp).squaredNorm();
					Scalar d1 = (mp1 - newp).squaredNorm();
					Vec3 newmp = d0 < d1 ? mp0 : mp1;
					Vertex cv = collapse_edge(m, e);
					value<Vec3>(m, helper.vertex_position_, cv) = newp;
					value<Vec3>(m, helper.vertex_medial_point_, cv) = newmp;
				}
			}
			return true;
		});
	} while (has_short_edge);

	foreach_cell(m, [&](Vertex v) -> bool {
		if (value<bool>(m, helper.vertex_is_fixed_, v))
			return true;
		uint32 count = 0;
		foreach_incident_edge(m, v, [&](Edge ie) -> bool {
			Scalar l = geometry::length(m, ie, helper.vertex_position_.get());
			if (l < helper.edge_collapse_threshold_ && !edge_can_collapse(m, ie))
				++count;
			return true;
		});
		bool is_fixed = count > 2;
		if (is_fixed)
		{
			value<Vec3>(m, helper.vertex_is_fixed_color_, v) = {1, 0, 0};
			value<bool>(m, helper.vertex_is_fixed_, v) = is_fixed;
		}
		return true;
	});
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SKELETON_H_
