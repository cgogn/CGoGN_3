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

#ifndef CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_
#define CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
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

inline void triangulate_incident_faces(CMap2& m, CMap2::Vertex v)
{
	std::vector<CMap2::Face> ifaces = incident_faces(m, v);
	for (CMap2::Face f : ifaces)
		cut_face(m, CMap2::Vertex(f.dart), CMap2::Vertex(phi<1, 1>(m, f.dart)));
}

inline bool edge_should_flip(CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Vertex> iv = incident_vertices(m, e);
	const int32 w = degree(m, iv[0]);
	const int32 x = degree(m, iv[1]);
	const int32 y = degree(m, CMap2::Vertex(phi<1, 1>(m, iv[0].dart)));
	const int32 z = degree(m, CMap2::Vertex(phi<1, 1>(m, iv[1].dart)));

	if (w < 4 || x < 4)
		return false;

	int32 dev_pre = abs(w - 6) + abs(x - 6) + abs(y - 6) + abs(z - 6);
	int32 dev_post = abs(w - 1 - 6) + abs(x - 1 - 6) + abs(y + 1 - 6) + abs(z + 1 - 6);
	return dev_pre - dev_post > 0;
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
struct PliantRemeshing_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	PliantRemeshing_Helper(MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
		: m_(m), vertex_position_(vertex_position), surface_bvh_(nullptr)
	{
		compute_bvh();
	}
	~PliantRemeshing_Helper()
	{
		if (surface_bvh_)
			delete surface_bvh_;
		if (feature_edge_)
			remove_attribute<Edge>(m_, feature_edge_);
		if (vertex_lfs_)
			remove_attribute<Vertex>(m_, vertex_lfs_);
	}

	void compute_bvh()
	{
		if (surface_bvh_)
			delete surface_bvh_;
		auto bvh_vertex_index = add_attribute<uint32, Vertex>(m_, "__bvh_vertex_index");
		std::vector<Vec3> bvh_vertex_position;
		bvh_vertex_position.reserve(nb_cells<Vertex>(m_));
		uint32 idx = 0;
		foreach_cell(m_, [&](Vertex v) -> bool {
			value<uint32>(m_, bvh_vertex_index, v) = idx++;
			bvh_vertex_position.push_back(value<Vec3>(m_, vertex_position_, v));
			return true;
		});
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_cells<Face>(m_) * 3);
		foreach_cell(m_, [&](Face f) -> bool {
			foreach_incident_vertex(m_, f, [&](Vertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(m_, bvh_vertex_index, v));
				return true;
			});
			return true;
		});
		surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, bvh_vertex_position);
		remove_attribute<Vertex>(m_, bvh_vertex_index);
	}

	void compute_lfs()
	{
		auto vertex_normal = add_attribute<Vec3, Vertex>(m_, "__vertex_normal");
		geometry::compute_normal<Vertex>(m_, vertex_position_.get(), vertex_normal.get());
		auto vertex_medial_point = add_attribute<Vec3, Vertex>(m_, "__vertex_medial_point");
		auto vertex_medial_point_radius = add_attribute<Scalar, Vertex>(m_, "__vertex_medial_point_radius");
		geometry::shrinking_ball_centers(m_, vertex_position_.get(), vertex_normal.get(), vertex_medial_point.get(),
										 vertex_medial_point_radius.get());

		vertex_lfs_ = get_or_add_attribute<Scalar, Vertex>(m_, "__vertex_lfs");
		lfs_min_ = std::numeric_limits<float64>::max();
		lfs_max_ = std::numeric_limits<float64>::lowest();
		lfs_mean_ = 0.0;
		uint32 nbv = 0;
		parallel_foreach_cell(m_, [&](Vertex v) -> bool {
			++nbv;
			uint32 vidx = index_of(m_, v);
			Scalar lfs = ((*vertex_medial_point)[vidx] - (*vertex_position_)[vidx]).norm();
			if (lfs < lfs_min_)
				lfs_min_ = lfs;
			if (lfs > lfs_max_)
				lfs_max_ = lfs;
			lfs_mean_ += lfs;
			(*vertex_lfs_)[vidx] = lfs;
			return true;
		});
		lfs_mean_ /= Scalar(nbv);

		remove_attribute<Vertex>(m_, vertex_normal);
		remove_attribute<Vertex>(m_, vertex_medial_point);
		remove_attribute<Vertex>(m_, vertex_medial_point_radius);

		std::cout << "lfs min: " << lfs_min_ << std::endl;
		std::cout << "lfs max: " << lfs_max_ << std::endl;
		std::cout << "lfs mean: " << lfs_mean_ << std::endl;
	}

	void detect_features()
	{
		Scalar angle_threshold = 60.0 * M_PI / 180.0;
		if (!feature_edge_)
			feature_edge_ = add_attribute<bool, Edge>(m_, "__feature_edge");
		if (!feature_vertex_)
			feature_vertex_ = add_attribute<bool, Vertex>(m_, "__feature_vertex");
		if (!feature_corner_)
			feature_corner_ = add_attribute<bool, Vertex>(m_, "__feature_corner");
		feature_edge_->fill(false);
		feature_vertex_->fill(false);
		feature_corner_->fill(false);
		parallel_foreach_cell(m_, [&](Edge e) -> bool {
			if (std::fabs(geometry::angle(m_, e, vertex_position_.get())) > angle_threshold)
			{
				value<bool>(m_, feature_edge_, e) = true;
				std::vector<Vertex> iv = incident_vertices(m_, e);
				value<bool>(m_, feature_vertex_, iv[0]) = true;
				value<bool>(m_, feature_vertex_, iv[1]) = true;
			}
			return true;
		});
		parallel_foreach_cell(m_, [&](Vertex v) -> bool {
			bool is_corner = false;
			uint32 nb_incident_feature_edge = 0;
			foreach_incident_edge(m_, v, [&](Edge e) -> bool {
				if (value<bool>(m_, feature_edge_, e))
					++nb_incident_feature_edge;
				if (nb_incident_feature_edge > 2)
					is_corner = true;
				return !is_corner;
			});
			value<bool>(m_, feature_corner_, v) = is_corner;
			return true;
		});
	}

	MESH& m_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	std::shared_ptr<Attribute<bool>> feature_edge_;
	std::shared_ptr<Attribute<bool>> feature_vertex_;
	std::shared_ptr<Attribute<bool>> feature_corner_;
	std::shared_ptr<Attribute<Scalar>> vertex_lfs_;
	Scalar lfs_min_, lfs_max_, lfs_mean_;
	acc::BVHTree<uint32, Vec3>* surface_bvh_;
};

template <typename MESH>
void pliant_remeshing(MESH& m, std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& vertex_position,
					  Scalar edge_length_target_ratio = 1.0, bool preserve_features = false, bool lfs_adaptive = false,
					  bool recompute_bvh = false)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	// static map to store helpers associated to meshes
	// allows to store context without polluting outer context and function api
	static std::unordered_map<MESH*, PliantRemeshing_Helper<MESH>> helpers_;
	auto [it, inserted] = helpers_.try_emplace(&m, m, vertex_position);
	PliantRemeshing_Helper<MESH>& helper = it->second;
	if (recompute_bvh || vertex_position != helper.vertex_position_)
		helper.compute_bvh();
	if (preserve_features)
		helper.detect_features();
	if (lfs_adaptive)
		helper.compute_lfs();

	Scalar edge_length_target = geometry::mean_edge_length(m, vertex_position.get()) * edge_length_target_ratio;

	const Scalar squared_min_edge_length = Scalar(0.5625) * edge_length_target * edge_length_target; // 0.5625 = 0.75^2
	const Scalar squared_max_edge_length = Scalar(1.5625) * edge_length_target * edge_length_target; // 1.5625 = 1.25^2

	CellCache<MESH> cache(m);

	for (uint32 i = 0; i < 3; ++i)
	{
		// cut long edges (and adjacent faces)
		bool has_long_edge = false;
		do
		{
			cache.template build<Edge>();
			has_long_edge = false;
			foreach_cell(cache, [&](Edge e) -> bool {
				std::vector<Vertex> iv = incident_vertices(m, e);
				Scalar lfs = 0.0; // init to zero for warning remove
				Scalar coeff = 1.0;
				if (lfs_adaptive)
				{
					lfs = (value<Scalar>(m, helper.vertex_lfs_, iv[0]) + value<Scalar>(m, helper.vertex_lfs_, iv[1])) *
						  0.5;
					if (lfs < helper.lfs_mean_)
						coeff = 0.25 + ((lfs - helper.lfs_min_) / (helper.lfs_mean_ - helper.lfs_min_) * 0.75);
					else
						coeff = 1 + ((lfs - helper.lfs_mean_) / (helper.lfs_max_ - helper.lfs_mean_) * 3.0);
				}
				Scalar threshold = squared_max_edge_length * coeff;
				if (geometry::squared_length(m, e, vertex_position.get()) > threshold)
				{
					has_long_edge = true;
					Vertex v = cut_edge(m, e);
					if (preserve_features)
					{
						if (value<bool>(m, helper.feature_edge_, e))
						{
							foreach_incident_edge(m, v, [&](Edge ie) -> bool {
								value<bool>(m, helper.feature_edge_, ie) = true;
								return true;
							});
						}
					}
					value<Vec3>(m, vertex_position, v) =
						(value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
					if (lfs_adaptive)
						value<Scalar>(m, helper.vertex_lfs_, v) = lfs;
					if (preserve_features)
					{
						value<bool>(m, helper.feature_corner_, v) = false;
						if (value<bool>(m, helper.feature_edge_, e))
							value<bool>(m, helper.feature_vertex_, v) = true;
					}
					triangulate_incident_faces(m, v);
				}
				return true;
			});
		} while (has_long_edge);

		// collapse short edges
		bool has_short_edge = false;
		do
		{
			has_short_edge = false;
			foreach_cell(m, [&](Edge e) -> bool {
				std::vector<Vertex> iv = incident_vertices(m, e);
				Scalar lfs;
				Scalar coeff = 1.0;
				if (lfs_adaptive)
				{
					lfs = (value<Scalar>(m, helper.vertex_lfs_, iv[0]) + value<Scalar>(m, helper.vertex_lfs_, iv[1])) *
						  0.5;
					if (lfs < helper.lfs_mean_)
						coeff = 0.25 + ((lfs - helper.lfs_min_) / (helper.lfs_mean_ - helper.lfs_min_) * 0.75);
					else
						coeff = 1 + ((lfs - helper.lfs_mean_) / (helper.lfs_max_ - helper.lfs_mean_) * 3.0);
				}
				Scalar threshold = squared_min_edge_length * coeff;
				if (geometry::squared_length(m, e, vertex_position.get()) < threshold)
				{
					bool collapse = true;
					const Vec3& p = value<Vec3>(m, vertex_position, iv[0]);
					foreach_adjacent_vertex_through_edge(m, iv[1], [&](Vertex v) -> bool {
						const Vec3& vec = p - value<Vec3>(m, vertex_position, v);
						if (vec.squaredNorm() > threshold)
							collapse = false;
						return collapse;
					});
					if (preserve_features)
					{
						if (value<bool>(m, helper.feature_corner_, iv[0]) ||
							value<bool>(m, helper.feature_corner_, iv[1]))
							collapse = false;
						if ((value<bool>(m, helper.feature_vertex_, iv[0]) &&
							 !value<bool>(m, helper.feature_vertex_, iv[1])) ||
							(!value<bool>(m, helper.feature_vertex_, iv[0]) &&
							 value<bool>(m, helper.feature_vertex_, iv[1])))
							collapse = false;
					}
					if (collapse && edge_can_collapse(m, e))
					{
						has_short_edge = true;
						Vec3 mp = value<Vec3>(m, vertex_position, iv[1]);
						// (value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1])) * 0.5;
						Vertex cv = collapse_edge(m, e);
						value<Vec3>(m, vertex_position, cv) = mp;
					}
				}
				return true;
			});
		} while (has_short_edge);

		// equalize valences with edge flips
		foreach_cell(m, [&](Edge e) -> bool {
			if (preserve_features && value<bool>(m, helper.feature_edge_, e))
				return true;

			if (!edge_can_flip(m, e))
				return true;

			if (edge_should_flip(m, e))
				flip_edge(m, e);
			else
			{
				// Delaunay flips
				std::vector<Vertex> iv = incident_vertices(m, e);
				if (degree(m, iv[0]) > 4 && degree(m, iv[1]) > 4)
				{
					std::vector<Scalar> op_angles = geometry::opposite_angles(m, e, vertex_position.get());
					if (op_angles[0] + op_angles[1] > M_PI)
						flip_edge(m, e);
				}
			}

			return true;
		});

		auto vertex_area = add_attribute<Scalar, Vertex>(m, "vertex_area__");
		geometry::compute_area<Vertex>(m, vertex_position.get(), vertex_area.get());

		// tangential relaxation
		// + project back on surface
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			if (is_incident_to_boundary(m, v))
				return true;
			Vec3 new_pos = value<Vec3>(m, vertex_position, v);
			if (preserve_features)
			{
				if (!value<bool>(m, helper.feature_corner_, v))
				{
					if (value<bool>(m, helper.feature_vertex_, v))
					{
						// Vec3 q(0, 0, 0);
						// uint32 count = 0;
						// foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
						// 	if (value<bool>(m, helper.feature_vertex_, av))
						// 	{
						// 		q += value<Vec3>(m, vertex_position, av);
						// 		++count;
						// 	}
						// 	return true;
						// });
						// if (count == 2)
						// {
						// 	q /= Scalar(count);
						// 	Vec3 n = geometry::normal(m, v, vertex_position.get());
						// 	new_pos = q + n.dot(value<Vec3>(m, vertex_position, v) - q) * n;
						// }
					}
					else
					{
						Vec3 q(0, 0, 0);
						Scalar total_area = 0.0;
						foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
							Scalar a = value<Scalar>(m, vertex_area, av);
							q += a * value<Vec3>(m, vertex_position, av);
							total_area += a;
							return true;
						});
						q /= Scalar(total_area);
						Vec3 n = geometry::normal(m, v, vertex_position.get());
						new_pos = q + n.dot(value<Vec3>(m, vertex_position, v) - q) * n;
						new_pos = helper.surface_bvh_->closest_point(new_pos);
					}
				}
			}
			else
			{
				Vec3 q(0, 0, 0);
				Scalar total_area = 0.0;
				foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
					Scalar a = value<Scalar>(m, vertex_area, av);
					q += a * value<Vec3>(m, vertex_position, av);
					total_area += a;
					return true;
				});
				q /= Scalar(total_area);
				Vec3 n = geometry::normal(m, v, vertex_position.get());
				new_pos = q + n.dot(value<Vec3>(m, vertex_position, v) - q) * n;
				new_pos = helper.surface_bvh_->closest_point(new_pos);
			}
			value<Vec3>(m, vertex_position, v) = new_pos;
			return true;
		});

		remove_attribute<Vertex>(m, vertex_area);
	}
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_REMESHING_PLIANT_REMESHING_H_
