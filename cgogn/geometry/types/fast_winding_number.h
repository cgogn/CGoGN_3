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

#ifndef CGOGN_GEOMETRY_TYPES_FAST_WINDING_NUMBER_H_
#define CGOGN_GEOMETRY_TYPES_FAST_WINDING_NUMBER_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <libacc/bvh_tree.h>

namespace cgogn
{

namespace geometry
{

inline Scalar frobenius_product(const Mat3& A, const Mat3& B)
{
	return A.cwiseProduct(B).sum();
}
//tensor 3rd order
struct Tensor3
{
	using Mat3 = Eigen::Matrix<Scalar, 3, 3>;
	std::array<Mat3, 3> M; // M[k](i,j) = S_ccn(i,j,k)

	Tensor3()
	{
		for (auto& A : M)
			A.setZero();
	}

	inline void add(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& n, Scalar a)
	{
		Mat3 C = (1.0 / 3.0) * (v1 * v1.transpose() + v2 * v2.transpose() + v3 * v3.transpose());
		M[0].noalias() += (a * n[0]) * C;
		M[1].noalias() += (a * n[1]) * C;
		M[2].noalias() += (a * n[2]) * C;
	}

	Tensor3& operator+=(const Tensor3& o)
	{
		for (int k = 0; k < 3; ++k)
			M[k] += o.M[k];
		return *this;
	}

	
	Scalar operator*(const Tensor3& T) const
	{
		Scalar sum = 0;
		for (int k = 0; k < 3; ++k)
			sum += frobenius_product(M[k], T.M[k]);
		return sum;
	}
	
};
template <int ORDER, typename SCALAR, typename VEC3>
struct Fast_Winding_Number_Coeff;

template <typename SCALAR, typename VEC3>
struct Fast_Winding_Number_Coeff<1, SCALAR, VEC3>
{
	using Fwn_Coeff = Fast_Winding_Number_Coeff<1, SCALAR, VEC3>;
	inline void add(const VEC3& c, const VEC3& n, SCALAR a)
	{
		//do nothing
	}                                                                                         

	Fwn_Coeff& operator+=(Fwn_Coeff const& other)
	{
		sum_area += other.sum_area;
		weighted_centroid += other.weighted_centroid;
		weighted_normal += other.weighted_normal;
		return *this;
	}

	friend Fwn_Coeff operator+(Fwn_Coeff lhs, const Fwn_Coeff& rhs)
	{
		lhs += rhs;
		return lhs;
	}
	Scalar sum_area = 0;
	VEC3 weighted_centroid = VEC3(0, 0, 0);
	VEC3 weighted_normal = VEC3(0, 0, 0);
};

template <typename SCALAR, typename VEC3>
struct Fast_Winding_Number_Coeff<2, SCALAR, VEC3>
{
	using Fwn_Coeff = Fast_Winding_Number_Coeff<2, SCALAR, VEC3>;
	inline void add(const VEC3& c, const VEC3& n, SCALAR a)
	{
		Q += a * (c * n.transpose());
	}

	Fwn_Coeff& operator+=(Fwn_Coeff const& other)
	{
		sum_area += other.sum_area;
		weighted_centroid += other.weighted_centroid;
		weighted_normal += other.weighted_normal;
		Q += other.Q;

		return *this;
	}

	friend Fwn_Coeff operator+(Fwn_Coeff lhs, const Fwn_Coeff& rhs)
	{
		lhs += rhs;
		return lhs;
	}
	Scalar sum_area = 0;
	VEC3 weighted_centroid = VEC3(0, 0, 0);
	VEC3 weighted_normal = VEC3(0, 0, 0);
	Eigen::Matrix<Scalar, 3, 3> Q = Eigen::Matrix<Scalar, 3, 3>::Zero();
};

template <typename Scalar, typename VEC3>
struct Fast_Winding_Number_Coeff<3, Scalar, VEC3>
{
	using Fwn_Coeff = Fast_Winding_Number_Coeff<3, Scalar, VEC3>;
	inline void add(const VEC3& v1, const VEC3& v2, const VEC3& v3, const VEC3& c, const VEC3& n, Scalar a)
	{
		Q += a * (c * n.transpose());
		T.add(v1, v2, v3, n, a);
	}

	Fwn_Coeff& operator+=(Fwn_Coeff const& other)
	{
		sum_area += other.sum_area;
		weighted_centroid += other.weighted_centroid;
		weighted_normal += other.weighted_normal;
		Q += other.Q;
		T += other.T;
		return *this;
	}

	friend Fwn_Coeff operator+(Fwn_Coeff lhs, const Fwn_Coeff& rhs)
	{
		lhs += rhs;
		return lhs;
	}

	Scalar sum_area = 0;
	VEC3 weighted_centroid = VEC3(0, 0, 0);
	VEC3 weighted_normal = VEC3(0, 0, 0);
	Eigen::Matrix<Scalar, 3, 3> Q = Eigen::Matrix<Scalar, 3, 3>::Zero();
	Tensor3 T;
};

template <typename MESH, int ORDER = 3>
class Fast_Winding_Number
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using BVH = acc::BVHTree<std::size_t, Vec3>;
	using Coeff = Fast_Winding_Number_Coeff<ORDER, Scalar, Vec3>;
	using Mat3 = Eigen::Matrix<Scalar, 3, 3>;

public:
	Fast_Winding_Number(MESH& m, const BVH& bvh_tree, Attribute<Vec3>* vertex_position, Attribute<Vec3>* face_normal,
						const std::vector<Face>& bvh_faces, Scalar beta)
		: m_(m), bvh_tree_(&bvh_tree), vertex_position_(vertex_position), face_normal_(face_normal),
		  bvh_faces_(bvh_faces), beta_(beta)
	{
		face_area_ = get_or_add_attribute<Scalar, Face>(m_, "f_area");
		face_centroid_ = get_or_add_attribute<Vec3, Face>(m_, "f_centroid");
		uint32 nb_node = bvh_tree_->node_count();
		coeffs_.resize(nb_node);
		p_tilde_.resize(nb_node, Vec3(0, 0, 0));
		node_radius_.resize(nb_node, 0);
		compute_area<Face>(m_, vertex_position_, face_area_.get());
		compute_centroid<Vec3, Face, MESH>(m_, vertex_position_, face_centroid_.get());
		precompute_coeffs();
	}
	~Fast_Winding_Number()
	{
		remove_attribute<Face>(m_, face_area_);
		remove_attribute<Face>(m_, face_centroid_);
	}
	
	Scalar evaluate_fast_winding_number(const Vec3& q) const
	{
		Scalar w = 0;
		std::stack<std::size_t> st;
		st.push(0);
		while (!st.empty())
		{
			std::size_t node_id = st.top();
			st.pop();
			Vec3 r = p_tilde_[node_id] - q;
			if (r.norm() > beta_ * node_radius_[node_id])
			{
				w += direct_eval(q, node_id);
			}
			else
			{
				if (bvh_tree_->is_leaf(node_id))
				{
					w += solid_angle_leaf(q, node_id);
				}
				else
				{
					auto [l_id, r_id] = bvh_tree_->children(node_id);
					st.push(l_id);
					st.push(r_id);
				}
			}
		}
		return w;
	}

	// Exact (ground truth) winding number by summing solid angles of all faces.
	// Returns value already normalized (divided by 4*pi) so comparable with evaluate_fast_winding_number.
	Scalar exact_winding_number(const Vec3& q) const
	{
		Scalar w = 0;
		for (const Face& f : bvh_faces_)
		{
			w += solid_angle(q, f) / (4 * M_PI);

		}
		return w;
	}

private:
	void precompute_coeffs()
	{
		std::stack<std::size_t> st;
		st.push(0);
		while (!st.empty())
		{
			std::size_t node_id = st.top();
			Coeff& coeff = coeffs_[node_id];
			st.pop();
			std::vector<std::size_t> primitives;
			primitives.reserve(nb_cells<Face>(m_));
			bvh_tree_->collect_primitives_id(node_id, primitives);
			Scalar sum_area = 0;
			Vec3 sum_ac(0, 0, 0); // sum of weighted centroid
			Vec3 sum_an(0, 0, 0); // sum of weighted normal
			// compute p_tilde
			for (std::size_t& fid : primitives)
			{
				Face f = bvh_faces_[fid];
				Scalar a = value<Scalar>(m_, face_area_, f);
				Vec3 c = value<Vec3>(m_, face_centroid_, f);
				Vec3 n = value<Vec3>(m_, face_normal_, f);
				sum_area += a;
				sum_ac += a * c;
				sum_an += a * n;
			}
			if (sum_area > 0)
			{
				p_tilde_[node_id] = sum_ac / sum_area;
			}
			Scalar max_norm = 0;
			for (std::size_t& fid : primitives)
			{
				Face f = bvh_faces_[fid];
				Scalar a = value<Scalar>(m_, face_area_, f);
				Vec3 c = value<Vec3>(m_, face_centroid_, f);
				Vec3 n = value<Vec3>(m_, face_normal_, f);
				Vec3 r = c - p_tilde_[node_id];
				auto iv = incident_vertices(m_, f);
				for (Vertex v : iv)
				{
					max_norm = std::max(max_norm, (value<Vec3>(m_, vertex_position_, v) - p_tilde_[node_id]).norm());
				}
				if constexpr (ORDER < 3)
					coeff.add(r, n, a);
				else
				{
					auto iv = incident_vertices(m_, f);
					Vec3 p1 = value<Vec3>(m_, vertex_position_, iv[0]);
					Vec3 p2 = value<Vec3>(m_, vertex_position_, iv[1]);
					Vec3 p3 = value<Vec3>(m_, vertex_position_, iv[2]);
					Vec3 v1 = 0.5 * (p1 + p2) - p_tilde_[node_id];
					Vec3 v2 = 0.5 * (p2 + p3) - p_tilde_[node_id];
					Vec3 v3 = 0.5 * (p3 + p1) - p_tilde_[node_id];

					coeff.add(v1, v2, v3, r, n, a);
				}
			}
			//node_radius_[node_id] = aabb_radius(node_id);
			node_radius_[node_id] = max_norm;
			/*std::cout << "Node computed radius: " << max_norm << ", AABB computed radius: " << aabb_radius(node_id)
					  << std::endl;*/
			coeff.sum_area = sum_area;
			coeff.weighted_centroid = sum_ac;
			coeff.weighted_normal = sum_an;
			if (!bvh_tree_->is_leaf(node_id))
			{
				auto [l_id, r_id] = bvh_tree_->children(node_id);
				st.push(l_id);
				st.push(r_id);
			}
			
		}
	}

	Tensor3 third_derivative(const Vec3& R, Scalar r) const 
	{
		Tensor3 K;

		const Scalar c1 = 15.0 / (4 * M_PI * std::pow(r, 7));
		const Scalar c2 = 3.0 / (4 * M_PI * std::pow(r, 5));
		for (int k = 0; k < 3; ++k)
		{
			Mat3 M = Mat3::Zero();
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					Scalar term1 = c1 * R[i] * R[j] * R[k];
					Scalar term2 = c2 * ((i == j ? R[k] : 0) + (i == k ? R[j] : 0) + (j == k ? R[i] : 0));
					M(i, j) = term1 - term2;
				}
			}
			K.M[k] = M;
		}
		return K;
	}
	
	Scalar solid_angle_leaf(const Vec3& q, const std::size_t node_id) const
	{
		Scalar w = 0.0;
		auto [first, last] = bvh_tree_->range(node_id);
		for (auto i = first; i < last; i++)
		{
			std::size_t fid = bvh_tree_->get_primitive_index(i);
			Face f = bvh_faces_[fid];

			w += solid_angle(q, f) / (4 * M_PI);
		}
		return w;
	}

	Scalar direct_eval(const Vec3& q, const std::size_t node_id) const
	{
		Scalar w = 0.0;
		const Vec3 R = p_tilde_[node_id] - q;
		Scalar r = std::max(R.norm(), 1e-20);
		if constexpr (ORDER >= 1)
		{
			Vec3 G1 = R / (4 * M_PI * r * r * r);
			w += coeffs_[node_id].weighted_normal.dot(G1);
		}
		if constexpr (ORDER >= 2)
		{
			Scalar r3 = r * r * r;
			Scalar r5 = r3 * r * r;
			Mat3 I = Mat3::Identity();
			Mat3 G2 = I / (4 * M_PI * r3) - 3 * R * R.transpose() / (4 * M_PI * r5);
			w += frobenius_product(coeffs_[node_id].Q, G2);
		}
		if constexpr (ORDER >= 3)
		{
			auto G3 = third_derivative(R, r);
			w += 0.5 * (coeffs_[node_id].T * G3);
		}
		return w;
	}

	Scalar solid_angle(const Vec3& q, const Face& f) const
	{
		auto iv = incident_vertices(m_, f);
		Vec3 v1 = value<Vec3>(m_, vertex_position_, iv[0]) - q;
		Vec3 v2 = value<Vec3>(m_, vertex_position_, iv[1]) - q;
		Vec3 v3 = value<Vec3>(m_, vertex_position_, iv[2]) - q;
		Scalar l1 = v1.norm();
		Scalar l2 = v2.norm();
		Scalar l3 = v3.norm();

		if (l1 == 0 || l2 == 0 || l3 == 0)
			return 0;
		v1 /= l1;
		v2 /= l2;
		v3 /= l3;

		const Scalar numerator = v1.dot((v2 - v1).cross(v3 - v1));
		if (numerator == 0)
			return 0;
		const Scalar denominator = 1 + v1.dot(v2) + v2.dot(v3) + v3.dot(v1);
		return 2 * std::atan2(numerator, denominator);
	}
	Scalar aabb_radius(const std::size_t node_id) const
	{
		Scalar r = 0.0;
		const auto box = bvh_tree_->node_aabb(node_id);
		Vec3 bb_min = box.min;
		Vec3 bb_max = box.max;
		Vec3 corner[8] = {Vec3(bb_min.x(), bb_min.y(), bb_min.z()), Vec3(bb_min.x(), bb_min.y(), bb_max.z()),
						  Vec3(bb_min.x(), bb_max.y(), bb_min.z()), Vec3(bb_min.x(), bb_max.y(), bb_max.z()),
						  Vec3(bb_max.x(), bb_min.y(), bb_min.z()), Vec3(bb_max.x(), bb_min.y(), bb_max.z()),
						  Vec3(bb_max.x(), bb_max.y(), bb_min.z()), Vec3(bb_max.x(), bb_max.y(), bb_max.z())};
		for (int i = 0; i < 8; i++)
		{
			Scalar d = (corner[i] - p_tilde_[node_id]).norm();
			if (d > r)
				r = d;
		}
		return r;
	}

private:
	Scalar beta_;
	const BVH* bvh_tree_;
	MESH& m_;
	const Attribute<Vec3>* vertex_position_;
	const Attribute<Vec3>* face_normal_;
	std::shared_ptr<Attribute<Scalar>> face_area_;
	std::shared_ptr<Attribute<Vec3>> face_centroid_;

	const std::vector<Face>& bvh_faces_;
	std::vector<Coeff> coeffs_;
	std::vector<Scalar> node_radius_; // radius of the aabb of the node
	std::vector<Vec3> p_tilde_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_QUADRIC_H_