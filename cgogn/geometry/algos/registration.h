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

#ifndef CGOGN_GEOMETRY_ALGOS_REGISTRATION_H_
#define CGOGN_GEOMETRY_ALGOS_REGISTRATION_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>
#include <simpleICP/simpleicp.h>

namespace cgogn
{

namespace geometry
{

template <typename MESHS, typename MESHT>
Mat4 rigid_register_mesh(MESHS& source, typename mesh_traits<MESHS>::template Attribute<Vec3>* source_vertex_position,
						 MESHT& target,
						 const typename mesh_traits<MESHT>::template Attribute<Vec3>* target_vertex_position)
{
	using SVertex = typename mesh_traits<MESHS>::Vertex;
	using TVertex = typename mesh_traits<MESHT>::Vertex;

	uint32 source_nbv = nb_cells<SVertex>(source);
	Eigen::MatrixXd X_source(source_nbv, 3);
	auto source_vertex_index = add_attribute<uint32, SVertex>(source, "__vertex_index");
	uint32 source_vertex_idx = 0;
	foreach_cell(source, [&](SVertex v) -> bool {
		const Vec3& p = value<Vec3>(source, source_vertex_position, v);
		X_source(source_vertex_idx, 0) = p[0];
		X_source(source_vertex_idx, 1) = p[1];
		X_source(source_vertex_idx, 2) = p[2];
		value<uint32>(source, source_vertex_index, v) = source_vertex_idx++;
		return true;
	});

	uint32 target_nbv = nb_cells<TVertex>(target);
	Eigen::MatrixXd X_target(target_nbv, 3);
	auto target_vertex_index = add_attribute<uint32, TVertex>(target, "__vertex_index");
	uint32 target_vertex_idx = 0;
	foreach_cell(target, [&](TVertex v) -> bool {
		const Vec3& p = value<Vec3>(target, target_vertex_position, v);
		X_target(target_vertex_idx, 0) = p[0];
		X_target(target_vertex_idx, 1) = p[1];
		X_target(target_vertex_idx, 2) = p[2];
		value<uint32>(target, target_vertex_index, v) = target_vertex_idx++;
		return true;
	});

	Mat4 t = SimpleICP(X_target, X_source);

	Eigen::MatrixXd X_sourceH(4, source_nbv);
	foreach_cell(source, [&](SVertex v) -> bool {
		uint32 idx = value<uint32>(source, source_vertex_index, v);
		const Vec3& p = value<Vec3>(source, source_vertex_position, v);
		X_sourceH(0, idx) = p[0];
		X_sourceH(1, idx) = p[1];
		X_sourceH(2, idx) = p[2];
		X_sourceH(3, idx) = 1.0;
		return true;
	});
	Eigen::MatrixXd res = t * X_sourceH;
	foreach_cell(source, [&](SVertex v) -> bool {
		uint32 idx = value<uint32>(source, source_vertex_index, v);
		Vec3& p = value<Vec3>(source, source_vertex_position, v);
		p[0] = res(0, idx);
		p[1] = res(1, idx);
		p[2] = res(2, idx);
		return true;
	});

	remove_attribute<SVertex>(source, source_vertex_index);
	remove_attribute<TVertex>(target, target_vertex_index);

	return t;
}

template <typename MESH>
struct NonRigidRegistration_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	struct BVH_Hit
	{
		bool hit = false;
		Face face;
		Vec3 bcoords;
		Scalar dist;
		Vec3 pos;
	};

	NonRigidRegistration_Helper(MESH& source, const std::shared_ptr<Attribute<Vec3>>& source_vertex_position,
								MESH& target, const Attribute<Vec3>* target_vertex_position, Scalar fit_to_target)
		: source_(source), source_vertex_position_(source_vertex_position), target_(&target),
		  target_vertex_position_(target_vertex_position), fit_to_target_(fit_to_target)
	{
		source_vertex_position_init_ = add_attribute<Vec3, Vertex>(source, "__nrrh_vertex_position_init");
		source_vertex_position_init_->copy(source_vertex_position_.get());

		source_vertex_rotation_matrix_ = add_attribute<Mat3, Vertex>(source, "__nrrh_vertex_rotation_matrix");
		Mat3 rm;
		rm.setZero();
		source_vertex_rotation_matrix_->fill(rm);

		build_target_bvh(target, target_vertex_position);

		// index source vertices
		source_vertex_index_ = add_attribute<uint32, Vertex>(source_, "__nrrh_vertex_index");
		source_nb_vertices_ = 0;
		foreach_cell(source_, [&](Vertex v) -> bool {
			value<uint32>(source_, source_vertex_index_, v) = source_nb_vertices_++;
			return true;
		});

		// compute source topo laplacian
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL(source_nb_vertices_, source_nb_vertices_);
		Eigen::MatrixXd vpos(source_nb_vertices_, 3);
		std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
		LAPLcoeffs.reserve(source_nb_vertices_ * 10);
		foreach_cell(source, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(source_, source_vertex_index_, v);
			const Vec3& pos = value<Vec3>(source_, source_vertex_position_, v);
			vpos(vidx, 0) = pos[0];
			vpos(vidx, 1) = pos[1];
			vpos(vidx, 2) = pos[2];
			uint32 nbv = 0;
			foreach_adjacent_vertex_through_edge(source_, v, [&](Vertex av) -> bool {
				uint32 avidx = value<uint32>(source_, source_vertex_index_, av);
				LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
				++nbv;
				return true;
			});
			LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(nbv)));
			return true;
		});
		LAPL.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());
		lapl_ = LAPL * vpos;
		rlapl_ = lapl_;

		// build solver
		build_solver(fit_to_target_);
	}
	~NonRigidRegistration_Helper()
	{
		remove_attribute<Vertex>(source_, source_vertex_index_);
		remove_attribute<Vertex>(source_, source_vertex_position_init_);
		remove_attribute<Vertex>(source_, source_vertex_rotation_matrix_);
		delete target_bvh_;
	}

	void init_source_steady_pos()
	{
		source_vertex_position_init_->copy(source_vertex_position_.get());
	}

	void build_target_bvh(MESH& target, const Attribute<Vec3>* target_vertex_position)
	{
		target_ = &target;
		target_vertex_position_ = target_vertex_position;

		auto target_vertex_index = add_attribute<uint32, Vertex>(target, "__nrrh_vertex_index");
		std::vector<Vec3> target_vertex_positions;
		uint32 target_nb_vertices = 0;
		foreach_cell(target, [&](Vertex v) -> bool {
			value<uint32>(target, target_vertex_index, v) = target_nb_vertices++;
			target_vertex_positions.push_back(value<Vec3>(*target_, target_vertex_position_, v));
			return true;
		});
		target_faces_.clear();
		target_faces_.reserve(target_nb_vertices * 2);
		std::vector<uint32> target_face_vertex_indices;
		target_face_vertex_indices.reserve(target_nb_vertices * 6);
		foreach_cell(*target_, [&](Face f) -> bool {
			target_faces_.push_back(f);
			foreach_incident_vertex(target, f, [&](Vertex v) -> bool {
				target_face_vertex_indices.push_back(value<uint32>(*target_, target_vertex_index, v));
				return true;
			});
			return true;
		});
		if (target_bvh_)
			delete target_bvh_;
		target_bvh_ = new acc::BVHTree<uint32, Vec3>(target_face_vertex_indices, target_vertex_positions);
		remove_attribute<Vertex>(target, target_vertex_index);
	}

	bool is_inside_target(const Vec3& p)
	{
		std::pair<uint32, Vec3> cp;
		target_bvh_->closest_point(p, &cp);
		Vec3 dir = (cp.second - p).normalized();
		Vec3 n = geometry::normal(*target_, target_faces_[cp.first], target_vertex_position_);
		return dir.dot(n) >= 0.0;
	}

	BVH_Hit intersect_target_bvh(const acc::Ray<Vec3>& r)
	{
		acc::BVHTree<uint32, Vec3>::Hit h;
		if (target_bvh_->intersect(r, &h))
		{
			Face f = target_faces_[h.idx];
			std::vector<Vertex> vertices = incident_vertices(*target_, f);
			Vec3 p = h.bcoords[0] * value<Vec3>(*target_, target_vertex_position_, vertices[0]) +
					 h.bcoords[1] * value<Vec3>(*target_, target_vertex_position_, vertices[1]) +
					 h.bcoords[2] * value<Vec3>(*target_, target_vertex_position_, vertices[2]);
			return {true, f, {h.bcoords[0], h.bcoords[1], h.bcoords[2]}, h.t, p};
		}
		else
			return BVH_Hit();
	}

	void build_solver(Scalar fit_to_target)
	{
		fit_to_target_ = fit_to_target;
		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(source_nb_vertices_ * 10);
		foreach_cell(source_, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(source_, source_vertex_index_, v);
			uint32 nbv = 0;
			foreach_adjacent_vertex_through_edge(source_, v, [&](Vertex av) -> bool {
				uint32 avidx = value<uint32>(source_, source_vertex_index_, av);
				Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
				++nbv;
				return true;
			});
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(nbv)));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(source_nb_vertices_ + vidx), int(vidx), fit_to_target_));
			return true;
		});
		A_.resize(2 * source_nb_vertices_, source_nb_vertices_);
		A_.setFromTriplets(Acoeffs.begin(), Acoeffs.end());
		At_ = A_.transpose();
		solver_.compute(At_ * A_);
	}

	MESH& source_;
	uint32 source_nb_vertices_;
	std::shared_ptr<Attribute<Vec3>> source_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> source_vertex_position_init_;
	std::shared_ptr<Attribute<Mat3>> source_vertex_rotation_matrix_;
	std::shared_ptr<Attribute<uint32>> source_vertex_index_;
	Eigen::MatrixXd lapl_;
	Eigen::MatrixXd rlapl_;

	MESH* target_ = nullptr;
	const Attribute<Vec3>* target_vertex_position_ = nullptr;

	Scalar fit_to_target_;
	acc::BVHTree<uint32, Vec3>* target_bvh_ = nullptr;
	std::vector<Face> target_faces_;

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A_;
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At_;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver_;
};

enum ProximityPolicy : uint32
{
	NEAREST_POINT,
	NORMAL_RAY
};

template <typename MESH>
void non_rigid_register_mesh(
	MESH& source, std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& source_vertex_position,
	MESH& target, const typename mesh_traits<MESH>::template Attribute<Vec3>* target_vertex_position,
	Scalar fit_to_target, bool relax, bool init_source_steady_pos, ProximityPolicy prox = NEAREST_POINT)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	// static map to store helpers associated to meshes
	// allows to store context without polluting outer context and function api
	static std::unordered_map<MESH*, NonRigidRegistration_Helper<MESH>> helpers_;
	auto [it, inserted] =
		helpers_.try_emplace(&source, source, source_vertex_position, target, target_vertex_position, fit_to_target);
	NonRigidRegistration_Helper<MESH>& helper = it->second;
	if (&target != helper.target_ || target_vertex_position != helper.target_vertex_position_)
		helper.build_target_bvh(target, target_vertex_position);
	if (fit_to_target != helper.fit_to_target_)
		helper.build_solver(fit_to_target);
	if (init_source_steady_pos)
		helper.init_source_steady_pos();

	// rotate laplacian coordinates
	parallel_foreach_cell(source, [&](Vertex v) -> bool {
		Mat3 cov;
		cov.setZero();
		const Vec3& pos = value<Vec3>(source, helper.source_vertex_position_, v);
		const Vec3& pos_i = value<Vec3>(source, helper.source_vertex_position_init_, v);
		foreach_adjacent_vertex_through_edge(source, v, [&](Vertex av) -> bool {
			Vec3 vec = (value<Vec3>(source, helper.source_vertex_position_, av) - pos).normalized();
			Vec3 vec_i = (value<Vec3>(source, helper.source_vertex_position_init_, av) - pos_i).normalized();
			for (uint32 i = 0; i < 3; ++i)
				for (uint32 j = 0; j < 3; ++j)
					cov(i, j) += vec[i] * vec_i[j];
			return true;
		});
		Eigen::JacobiSVD<Mat3> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Mat3 R = svd.matrixU() * svd.matrixV().transpose();
		if (R.determinant() < 0)
		{
			Mat3 U = svd.matrixU();
			for (uint32 i = 0; i < 3; ++i)
				U(i, 2) *= -1;
			R = U * svd.matrixV().transpose();
		}
		value<Mat3>(source, helper.source_vertex_rotation_matrix_, v) = R;
		return true;
	});
	parallel_foreach_cell(source, [&](Vertex v) -> bool {
		// uint32 degree = 0;
		// Mat3 r;
		// r.setZero();
		// foreach_adjacent_vertex_through_edge(source, v, [&](Vertex av) -> bool {
		// 	r += value<Mat3>(source, helper.source_vertex_rotation_matrix_, av);
		// 	++degree;
		// 	return true;
		// });
		// r += value<Mat3>(source, helper.source_vertex_rotation_matrix_, v);
		// r /= degree + 1;
		const Mat3& r = value<Mat3>(source, helper.source_vertex_rotation_matrix_, v);
		uint32 vidx = value<uint32>(source, helper.source_vertex_index_, v);
		Vec3 l;
		l[0] = helper.lapl_(vidx, 0);
		l[1] = helper.lapl_(vidx, 1);
		l[2] = helper.lapl_(vidx, 2);
		Vec3 rl = r * l;
		helper.rlapl_(vidx, 0) = rl[0];
		helper.rlapl_(vidx, 1) = rl[1];
		helper.rlapl_(vidx, 2) = rl[2];
		return true;
	});

	// setup RHS
	Eigen::MatrixXd b(2 * helper.source_nb_vertices_, 3);
	Vec3 pos;
	parallel_foreach_cell(source, [&](Vertex v) -> bool {
		Vec3 pos;
		switch (prox)
		{
		case NEAREST_POINT: {
			pos = helper.target_bvh_->closest_point(value<Vec3>(source, helper.source_vertex_position_, v));
		}
		break;
		case NORMAL_RAY: {
			const Vec3& p = value<Vec3>(source, helper.source_vertex_position_, v);
			Vec3 n{0, 0, 0};
			foreach_incident_face(source, v, [&](Face f) -> bool {
				Vec3 nf = geometry::normal(source, f, helper.source_vertex_position_.get());
				Vec3 cf = geometry::centroid<Vec3>(source, f, helper.source_vertex_position_.get());
				bool inside = helper.is_inside_target(cf);
				if (!inside)
					nf *= -1;
				typename NonRigidRegistration_Helper<MESH>::BVH_Hit h =
					helper.intersect_target_bvh({cf, nf, 0, acc::inf});
				if (h.hit)
					n += inside ? h.pos - cf : cf - h.pos;
				return true;
			});
			n.normalize();

			if (!helper.is_inside_target(p))
				n *= -1;

			typename NonRigidRegistration_Helper<MESH>::BVH_Hit h = helper.intersect_target_bvh({p, n, 0, acc::inf});
			if (h.hit)
				pos = h.pos;
			else
				pos = helper.target_bvh_->closest_point(p);
		}
		break;
		};

		uint32 vidx = value<uint32>(source, helper.source_vertex_index_, v);

		b(vidx, 0) = relax ? 0.0 : helper.rlapl_(vidx, 0);
		b(vidx, 1) = relax ? 0.0 : helper.rlapl_(vidx, 1);
		b(vidx, 2) = relax ? 0.0 : helper.rlapl_(vidx, 2);
		b(helper.source_nb_vertices_ + vidx, 0) = fit_to_target * pos[0];
		b(helper.source_nb_vertices_ + vidx, 1) = fit_to_target * pos[1];
		b(helper.source_nb_vertices_ + vidx, 2) = fit_to_target * pos[2];

		return true;
	});

	// solve
	Eigen::MatrixXd vpos = helper.solver_.solve(helper.At_ * b);

	// store result
	parallel_foreach_cell(source, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(source, helper.source_vertex_index_, v);
		Vec3& pos = value<Vec3>(source, helper.source_vertex_position_, v);
		pos[0] = vpos(vidx, 0);
		pos[1] = vpos(vidx, 1);
		pos[2] = vpos(vidx, 2);
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_REGISTRATION_H_
