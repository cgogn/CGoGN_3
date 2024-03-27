/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_SPHERE_FITTING_H_
#define CGOGN_MODULE_SPHERE_FITTING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>

#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>
#include <libacc/kd_tree.h>

#include <GLFW/glfw3.h>

#include <boost/synapse/connect.hpp>
#include <set>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD>
class SphereFitting : public ViewModule
{
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SEdge = typename mesh_traits<SURFACE>::Edge;
	using SFace = typename mesh_traits<SURFACE>::Face;

	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	template <typename T>
	using NMAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	enum UpdateMethod : uint32
	{
		MEAN,
		FIT
	};

	struct SurfaceParameters
	{
		bool initialized_ = false;

		SURFACE* surface_;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_normal_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_face_normal_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> medial_axis_position_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> medial_axis_radius_ = nullptr;
		std::shared_ptr<SAttribute<SVertex>> medial_axis_secondary_vertex_ = nullptr;
		std::shared_ptr<SAttribute<bool>> medial_axis_selected_ = nullptr;
		std::shared_ptr<SAttribute<PVertex>> surface_vertex_cluster_ = nullptr;

		acc::BVHTree<uint32, Vec3>* surface_bvh_ = nullptr;
		std::vector<SFace> surface_bvh_faces_;
		acc::KDTree<3, uint32>* surface_kdt_ = nullptr;
		std::vector<SVertex> surface_kdt_vertices_;

		POINTS* spheres_;
		std::shared_ptr<PAttribute<Vec3>> spheres_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_color_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<SVertex>>> spheres_cluster_ = nullptr;
		std::shared_ptr<PAttribute<std::set<PVertex>>> spheres_neighbor_clusters_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_deviation_ = nullptr;
		PAttribute<Scalar>* selected_spheres_error_ = nullptr;

		NONMANIFOLD* skeleton_;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;

		float32 filter_radius_threshold_ = 0.01f;
		float32 filter_angle_threshold_ = M_PI / 4.0f;

		float32 init_min_radius_ = 0.01f;
		float32 init_dilation_factor_ = 1.75f;

		UpdateMethod update_method_ = FIT;
		float32 mean_update_curvature_weight_ = 0.01f;
		bool auto_split_outside_spheres_ = false;

		Scalar total_error_ = 0.0;
		Scalar min_error_ = 0.0;
		Scalar max_error_ = 0.0;

		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool slow_down_ = true;
		uint32 update_rate_ = 20;
	};

public:
	SphereFitting(const App& app)
		: ViewModule(app, "SphereFitting (" + std::string{mesh_traits<SURFACE>::name} + "," +
							  std::string{mesh_traits<POINTS>::name} + ")")
	{
	}
	~SphereFitting()
	{
	}

	void set_selected_surface(SURFACE& s)
	{
		selected_surface_ = &s;
		picked_sphere_ = PVertex();
	}

	void set_surface_vertex_position(SURFACE& s, const std::shared_ptr<SAttribute<Vec3>>& surface_vertex_position)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_vertex_position_ = surface_vertex_position;
	}

	void init_surface_data(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_ = &s;

		if (!p.surface_vertex_position_)
		{
			std::cout << "No surface vertex position attribute set" << std::endl;
			return;
		}

		if (p.initialized_)
		{
			std::cout << "Surface data already initialized" << std::endl;
			return;
		}

		// create BVH and KDTree for the surface

		MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
		uint32 nb_vertices = md.template nb_cells<SVertex>();
		uint32 nb_faces = md.template nb_cells<SFace>();

		auto bvh_vertex_index = add_attribute<uint32, SVertex>(s, "__bvh_vertex_index");

		p.surface_kdt_vertices_.clear();
		p.surface_kdt_vertices_.reserve(nb_vertices);
		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(s, [&](SVertex v) -> bool {
			p.surface_kdt_vertices_.push_back(v);
			value<uint32>(s, bvh_vertex_index, v) = idx++;
			vertex_position_vector.push_back(value<Vec3>(s, p.surface_vertex_position_, v));
			return true;
		});

		p.surface_bvh_faces_.clear();
		p.surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(s, [&](SFace f) -> bool {
			p.surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(s, f, [&](SVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(s, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		p.surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);
		p.surface_kdt_ = new acc::KDTree<3, uint32>(vertex_position_vector);

		remove_attribute<SVertex>(s, bvh_vertex_index);

		// compute normals

		p.surface_face_normal_ = get_or_add_attribute<Vec3, SFace>(s, "normal");
		geometry::compute_normal<SFace>(s, p.surface_vertex_position_.get(), p.surface_face_normal_.get());

		p.surface_vertex_normal_ = get_or_add_attribute<Vec3, SVertex>(s, "normal");
		geometry::compute_normal<SVertex>(s, p.surface_vertex_position_.get(), p.surface_vertex_normal_.get());

		// compute shrinking balls for the surface vertices

		p.medial_axis_position_ = get_or_add_attribute<Vec3, SVertex>(s, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, SVertex>(s, "medial_axis_radius");
		p.medial_axis_secondary_vertex_ = get_or_add_attribute<SVertex, SVertex>(s, "medial_axis_secondary_vertex_");

		parallel_foreach_cell(s, [&](SVertex v) -> bool {
			uint32 v_index = index_of(s, v);
			auto [c, r, q] = geometry::shrinking_ball_center(
				s, (*p.surface_vertex_position_)[v_index], (*p.surface_vertex_normal_)[v_index],
				p.surface_vertex_position_.get(), p.surface_bvh_, p.surface_bvh_faces_, p.surface_kdt_,
				p.surface_kdt_vertices_);
			(*p.medial_axis_position_)[v_index] = c;
			(*p.medial_axis_radius_)[v_index] = r;
			(*p.medial_axis_secondary_vertex_)[v_index] = q;
			return true;
		});

		// filter the medial samples

		p.medial_axis_selected_ = get_or_add_attribute<bool, SVertex>(s, "medial_axis_selected");
		filter_medial_samples(p);

		// create the spheres mesh

		p.spheres_ = points_provider_->add_mesh(surface_provider_->mesh_name(s) + "_spheres");
		p.spheres_position_ = add_attribute<Vec3, PVertex>(*p.spheres_, "position");
		p.spheres_radius_ = add_attribute<Scalar, PVertex>(*p.spheres_, "radius");
		p.spheres_color_ = add_attribute<Vec4, PVertex>(*p.spheres_, "color");
		p.spheres_error_ = add_attribute<Scalar, PVertex>(*p.spheres_, "error");
		p.spheres_deviation_ = add_attribute<Scalar, PVertex>(*p.spheres_, "deviation");
		p.selected_spheres_error_ = p.spheres_error_.get();

		// cluster info

		p.spheres_cluster_ =
			add_attribute<std::vector<SVertex>, PVertex>(*p.spheres_, "cluster"); // surface vertices in the cluster
		p.surface_vertex_cluster_ =
			get_or_add_attribute<PVertex, SVertex>(s, "cluster"); // cluster of the surface vertex
		p.spheres_neighbor_clusters_ =
			get_or_add_attribute<std::set<PVertex>, PVertex>(s, "neighbor_clusters"); // neighbor clusters

		// create the skeleton mesh

		p.skeleton_ = non_manifold_provider_->add_mesh(surface_provider_->mesh_name(s) + "_skeleton");
		p.skeleton_position_ = add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");

		p.initialized_ = true;
	}

	void filter_medial_samples(SurfaceParameters& p)
	{
		std::lock_guard<std::mutex> lock(p.mutex_);

		parallel_foreach_cell(*p.surface_, [&](SVertex v) -> bool {
			uint32 v_index = index_of(*p.surface_, v);
			const Vec3& c = (*p.medial_axis_position_)[v_index];
			SVertex sv = (*p.medial_axis_secondary_vertex_)[v_index];
			if (sv.is_valid())
			{
				const Vec3& c1 = (*p.surface_vertex_position_)[v_index];
				const Vec3& c2 = value<Vec3>(*p.surface_, p.surface_vertex_position_, sv);
				const Scalar r = (*p.medial_axis_radius_)[v_index];

				(*p.medial_axis_selected_)[v_index] =
					r > p.filter_radius_threshold_ && geometry::angle(c1 - c, c2 - c) > p.filter_angle_threshold_;
			}
			else
				(*p.medial_axis_selected_)[v_index] = false;

			return true;
		});

		if (p.initialized_ && !p.running_)
			update_render_data(p);
	}

	void set_selected_spheres_error(SurfaceParameters& p, PAttribute<Scalar>* attribute)
	{
		p.selected_spheres_error_ = attribute;
		update_spheres_color(p);
		if (!p.running_)
			update_render_data(p);
	}

	uint32 nb_covered_vertices(SurfaceParameters& p, SVertex v)
	{
		uint32 v_index = index_of(*p.surface_, v);
		const Vec3& vp = (*p.medial_axis_position_)[v_index];
		Scalar vr = (*p.medial_axis_radius_)[v_index];
		uint32 count = 0;
		foreach_cell(*p.surface_, [&](SVertex w) -> bool {
			if ((value<Vec3>(*p.surface_, p.surface_vertex_position_, w) - vp).norm() < p.init_dilation_factor_ * vr)
				++count;
			return true;
		});
		return count;
	}

	void init_spheres(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];

		std::lock_guard<std::mutex> lock(p.mutex_);

		clear(*p.spheres_);

		// auto nb_covered = add_attribute<uint32, SVertex>(s, "__nb_covered");
		// parallel_foreach_cell(s, [&](SVertex v) -> bool {
		// 	value<uint32>(s, nb_covered, v) = nb_covered_vertices(s, v);
		// 	return true;
		// });

		std::vector<SVertex> sorted_vertices;
		foreach_cell(*p.surface_, [&](SVertex v) -> bool {
			sorted_vertices.push_back(v);
			return true;
		});
		std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&](SVertex a, SVertex b) {
			// sort candidate spheres by decreasing radius
			return value<Scalar>(*p.surface_, p.medial_axis_radius_, a) >
				   value<Scalar>(*p.surface_, p.medial_axis_radius_, b);
			// sort candidate spheres by decreasing number of covered vertices
			// return value<uint32>(*p.surface_, nb_covered, a) > value<uint32>(*p.surface_, nb_covered, b);
		});

		// remove_attribute<SVertex>(*p.surface_, nb_covered);

		auto covered = add_attribute<bool, SVertex>(*p.surface_, "__covered");
		covered->fill(false);

		for (SVertex v : sorted_vertices)
		{
			uint32 v_index = index_of(*p.surface_, v);

			// do not add spheres with radius smaller than init_min_radius_
			if ((*p.medial_axis_radius_)[v_index] < p.init_min_radius_)
				break;

			if ((*covered)[v_index])
				continue;

			if (!(*p.medial_axis_selected_)[v_index])
				continue;

			const Vec3& vp = (*p.medial_axis_position_)[v_index];
			Scalar vr = (*p.medial_axis_radius_)[v_index];

			PVertex sphere = add_vertex(*p.spheres_);
			value<Vec3>(*p.spheres_, p.spheres_position_, sphere) = vp;
			value<Scalar>(*p.spheres_, p.spheres_radius_, sphere) = vr;

			std::stack<SVertex> stack;
			stack.push(v);
			while (!stack.empty())
			{
				SVertex w = stack.top();
				stack.pop();
				value<bool>(*p.surface_, covered, w) = true;
				foreach_adjacent_vertex_through_edge(*p.surface_, w, [&](SVertex u) -> bool {
					uint32 u_index = index_of(*p.surface_, u);
					if (!(*covered)[u_index] &&
						((*p.surface_vertex_position_)[u_index] - vp).norm() < p.init_dilation_factor_ * vr)
						stack.push(u);
					return true;
				});
			}
			stack.push((*p.medial_axis_secondary_vertex_)[v_index]);
			while (!stack.empty())
			{
				SVertex w = stack.top();
				stack.pop();
				value<bool>(*p.surface_, covered, w) = true;
				foreach_adjacent_vertex_through_edge(*p.surface_, w, [&](SVertex u) -> bool {
					uint32 u_index = index_of(*p.surface_, u);
					if (!(*covered)[u_index] &&
						((*p.surface_vertex_position_)[u_index] - vp).norm() < p.init_dilation_factor_ * vr)
						stack.push(u);
					return true;
				});
			}
		}

		remove_attribute<SVertex>(*p.surface_, covered);

		compute_clusters(p);
		compute_spheres_error(p);

		if (!p.running_)
			update_render_data(p);
	}

	void compute_clusters(SurfaceParameters& p)
	{
		// clean cluster affectation
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_neighbor_clusters_)[v_index].clear();
			return true;
		});
		p.surface_vertex_cluster_->fill(PVertex());

		// auto start = std::chrono::high_resolution_clock::now();

		// build spheres BVH
		MeshData<POINTS>& md = points_provider_->mesh_data(*p.spheres_);
		uint32 nb_spheres = md.template nb_cells<PVertex>();
		std::vector<Vec3> sphere_centers;
		sphere_centers.reserve(nb_spheres);
		std::vector<Scalar> sphere_radii;
		sphere_radii.reserve(nb_spheres);
		std::vector<PVertex> spheres_bvh_vertices;
		spheres_bvh_vertices.reserve(nb_spheres);
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			spheres_bvh_vertices.push_back(v);
			sphere_centers.push_back((*p.spheres_position_)[v_index]);
			sphere_radii.push_back((*p.spheres_radius_)[v_index]);
			return true;
		});
		acc::BVHTreeSpheres<uint32, Vec3> spheres_bvh(sphere_centers, sphere_radii);

		// assign each surface vertex to its closest sphere cluster
		parallel_foreach_cell(*p.surface_, [&](SVertex v) -> bool {
			uint32 v_index = index_of(*p.surface_, v);
			if (!(*p.medial_axis_selected_)[v_index])
				return true;

			const Vec3& vp = (*p.surface_vertex_position_)[v_index];
			PVertex closest_sphere;

			std::pair<uint32, Vec3> bvh_res;
			spheres_bvh.closest_point(vp, &bvh_res);
			closest_sphere = spheres_bvh_vertices[bvh_res.first];

			(*p.surface_vertex_cluster_)[v_index] = closest_sphere;

			std::lock_guard<std::mutex> lock(spheres_mutex_[bvh_res.first % spheres_mutex_.size()]);
			value<std::vector<SVertex>>(*p.spheres_, p.spheres_cluster_, closest_sphere).push_back(v);

			return true;
		});

		// auto end = std::chrono::high_resolution_clock::now();

		// std::cout << "Cluster computation time: " << std::chrono::duration<Scalar>(end - start).count() << "s"
		// 		  << std::endl;

		// remove small clusters
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			std::vector<SVertex>& cluster = value<std::vector<SVertex>>(*p.spheres_, p.spheres_cluster_, v);
			if (cluster.size() < 4)
			{
				// std::cout << "small cluster - removed sphere" << std::endl;
				for (SVertex sv : cluster)
					value<PVertex>(*p.surface_, p.surface_vertex_cluster_, sv) = PVertex();
				remove_vertex(*p.spheres_, v);
			}
			return true;
		});

		// compute neighbor clusters
		foreach_cell(*p.surface_, [&](SEdge e) -> bool {
			std::vector<SVertex> vertices = incident_vertices(*p.surface_, e);

			PVertex v1_cluster = value<PVertex>(*p.surface_, p.surface_vertex_cluster_, vertices[0]);
			PVertex v2_cluster = value<PVertex>(*p.surface_, p.surface_vertex_cluster_, vertices[1]);
			if (v1_cluster.is_valid() && v2_cluster.is_valid() && v1_cluster != v2_cluster)
			{
				value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, v1_cluster).insert(v2_cluster);
				value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, v2_cluster).insert(v1_cluster);
			}

			return true;
		});
	}

	void compute_spheres_error(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) {
			uint32 v_index = index_of(*p.spheres_, v);

			Scalar error = 0.0;
			const Vec3& center = (*p.spheres_position_)[v_index];
			Scalar radius = (*p.spheres_radius_)[v_index];
			const std::vector<SVertex>& cluster = (*p.spheres_cluster_)[v_index];
			for (SVertex sv : cluster)
			{
				Scalar d = (value<Vec3>(*p.surface_, p.surface_vertex_position_, sv) - center).norm() - radius;
				error += d * d;
			}
			(*p.spheres_error_)[v_index] = error;

			Scalar mean_error = error / cluster.size();
			Scalar variance = 0.0;
			for (SVertex sv : cluster)
			{
				Scalar d = (value<Vec3>(*p.surface_, p.surface_vertex_position_, sv) - center).norm() - radius;
				Scalar err = d * d;
				variance += (err - mean_error) * (err - mean_error);
			}
			variance /= cluster.size();
			Scalar deviation = std::sqrt(variance);
			(*p.spheres_deviation_)[v_index] = deviation;

			return true;
		});

		p.min_error_ = std::numeric_limits<Scalar>::max();
		p.max_error_ = std::numeric_limits<Scalar>::min();
		p.total_error_ = 0.0;
		for (Scalar e : *p.selected_spheres_error_)
		{
			p.min_error_ = std::min(p.min_error_, e);
			p.max_error_ = std::max(p.max_error_, e);
			p.total_error_ += e;
		}
	}

	void update_spheres_color(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_color_)[v_index] = color_map((*p.selected_spheres_error_)[v_index], p.min_error_, p.max_error_);
			return true;
		});
	}

	Vec4 color_map(Scalar x, Scalar min, Scalar max)
	{
		x = (x - min) / (max - min);
		x = std::clamp(x, 0.0, 1.0);

		Scalar x2 = 2.0 * x;
		switch (int(std::floor(std::max(0.0, x2 + 1.0))))
		{
		case 0:
			return Vec4(0.0, 0.0, 1.0, 0.5);
		case 1:
			return Vec4(x2, x2, 1.0, 0.5);
		case 2:
			return Vec4(1.0, 2.0 - x2, 2.0 - x2, 0.5);
		}
		return Vec4(1.0, 0.0, 0.0, 0.5);
	}

	// search for the 25 closest vertices to the exterior center "center" and keep the first (p1) and the one among the
	// others which maximizes the angle [p1-pos, p2-pos] (p2). Then:
	// - update the given sphere with center and radius of the shrinking ball adjacent to p1
	// - adds a new sphere with center and radius of the shrinking ball adjacent to p2
	void split_outside_sphere(SurfaceParameters& p, const Vec3& center, PVertex sphere)
	{
		std::cout << "outside center" << std::endl;

		std::vector<std::pair<uint32, Scalar>> k_res;
		p.surface_kdt_->find_nns(center, 25, &k_res);

		// look for the closest vertex with maximal angle with the first closest vertex
		Scalar max_angle = 0.0;
		SVertex max_angle_vertex;
		for (uint32 i = 1; i < k_res.size(); ++i)
		{
			Scalar a = geometry::angle(p.surface_kdt_->vertex(k_res[0].first) - center,
									   p.surface_kdt_->vertex(k_res[i].first) - center);
			if (a > max_angle)
			{
				max_angle = a;
				max_angle_vertex = p.surface_kdt_vertices_[k_res[i].first];
			}
		}

		const Vec3& pos0 =
			value<Vec3>(*p.surface_, p.surface_vertex_position_, p.surface_kdt_vertices_[k_res[0].first]);
		const Vec3& pos0n = value<Vec3>(*p.surface_, p.surface_vertex_normal_, p.surface_kdt_vertices_[k_res[0].first]);

		auto [center0, radius0, secondary0] =
			geometry::shrinking_ball_center(*p.surface_, pos0, pos0n, p.surface_vertex_position_.get(), p.surface_bvh_,
											p.surface_bvh_faces_, p.surface_kdt_, p.surface_kdt_vertices_);

		value<Vec3>(*p.spheres_, p.spheres_position_, sphere) = center0;
		value<Scalar>(*p.spheres_, p.spheres_radius_, sphere) = radius0;

		const Vec3& pos1 = value<Vec3>(*p.surface_, p.surface_vertex_position_, max_angle_vertex);
		const Vec3& pos1n = value<Vec3>(*p.surface_, p.surface_vertex_normal_, max_angle_vertex);

		auto [center1, radius1, secondary1] =
			geometry::shrinking_ball_center(*p.surface_, pos1, pos1n, p.surface_vertex_position_.get(), p.surface_bvh_,
											p.surface_bvh_faces_, p.surface_kdt_, p.surface_kdt_vertices_);

		PVertex new_sphere = add_vertex(*p.spheres_);
		value<Vec3>(*p.spheres_, p.spheres_position_, new_sphere) = center1;
		value<Scalar>(*p.spheres_, p.spheres_radius_, new_sphere) = radius1;
	}

	void update_sphere_mean(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<SVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = Vec3(0.0, 0.0, 0.0);
		Scalar n = 0.0;

		for (SVertex v : cluster)
		{
			uint32 v_index = index_of(*p.surface_, v);
			Scalar r = (*p.medial_axis_radius_)[v_index];
			// avoid division by zero (or too small value)
			r = std::max(r, 0.001);
			// w varies between 1 and r when mean_update_curvature_weight_ varies between 0 and 1
			Scalar w = 1.0 + (1.0 - p.mean_update_curvature_weight_) * (r - 1.0);
			c += (w / r) * (*p.medial_axis_position_)[v_index];
			n += (w / r);
		}
		c /= n;

		std::pair<uint32, Vec3> bvh_res;
		p.surface_bvh_->closest_point(c, &bvh_res);
		Vec3 closest_point_position = bvh_res.second;
		Scalar closest_point_dist = (closest_point_position - c).norm();
		Vec3 closest_point_dir = (closest_point_position - c) / closest_point_dist;

		const Vec3& closest_face_normal =
			value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
		// TODO: exterior detection is not reliable
		if (closest_point_dir.dot(closest_face_normal) <= 0.0)
		{
			if (p.auto_split_outside_spheres_)
			{
				split_outside_sphere(p, c, sphere);
				return;
			}
			else
				closest_point_dir = -closest_point_dir;
		}

		auto [center, radius, secondary] = geometry::shrinking_ball_center(
			*p.surface_, closest_point_position, closest_point_dir, p.surface_vertex_position_.get(), p.surface_bvh_,
			p.surface_bvh_faces_, p.surface_kdt_, p.surface_kdt_vertices_);

		(*p.spheres_position_)[sphere_index] = center;
		(*p.spheres_radius_)[sphere_index] = radius;
	}

	void update_sphere_fit(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<SVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Eigen::MatrixXd A(cluster.size(), 4);
		Eigen::VectorXd b(cluster.size());
		uint32 idx = 0;
		for (SVertex v : cluster)
		{
			const Vec3& pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, v);
			A.row(idx) = Eigen::Vector4d(-2.0 * pos[0], -2.0 * pos[1], -2.0 * pos[2], 1.0);
			b(idx) = -(pos[0] * pos[0]) - (pos[1] * pos[1]) - (pos[2] * pos[2]);
			++idx;
		};
		Eigen::LDLT<Eigen::MatrixXd> solver(A.transpose() * A);
		Eigen::MatrixXd s1 = solver.solve(A.transpose() * b);
		s1(3) = std::sqrt(s1(0) * s1(0) + s1(1) * s1(1) + s1(2) * s1(2) - s1(3));

		Vec3 s1c = Vec3(s1(0), s1(1), s1(2));
		Scalar s1r = s1(3);

		Eigen::MatrixXd J(cluster.size(), 4);
		Eigen::VectorXd r(cluster.size());
		Eigen::VectorXd s2(4);
		s2 << s1(0), s1(1), s1(2), s1(3);
		for (uint32 i = 0; i < 5; ++i) // TODO: check number of iterations
		{
			idx = 0;
			for (SVertex v : cluster)
			{
				const Vec3& pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, v);
				Vec3 d = pos - Vec3(s2(0), s2(1), s2(2));
				Scalar l = d.norm();
				J.row(idx) = Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0);
				r(idx) = -(l - s2(3));
				++idx;
			}
			Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
			s2 += solver.solve(J.transpose() * r);
		}

		Vec3 s2c = Vec3(s2(0), s2(1), s2(2));
		Scalar s2r = s2(3);

		std::pair<uint32, Vec3> bvh_res;
		p.surface_bvh_->closest_point(s2c, &bvh_res);
		Vec3 closest_point_position = bvh_res.second;
		Vec3 closest_point_dir = (closest_point_position - s2c).normalized();

		Vec3 closest_face_normal =
			value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
		// TODO: exterior detection is not reliable
		if (closest_point_dir.dot(closest_face_normal) <= 0.0)
		{
			if (p.auto_split_outside_spheres_)
			{
				split_outside_sphere(p, s2c, sphere);
				return;
			}
			else
				closest_point_dir = -closest_point_dir;
		}

		auto [center, radius, secondary] = geometry::shrinking_ball_center(
			*p.surface_, closest_point_position, closest_point_dir, p.surface_vertex_position_.get(), p.surface_bvh_,
			p.surface_bvh_faces_, p.surface_kdt_, p.surface_kdt_vertices_);

		(*p.spheres_position_)[sphere_index] = center;
		(*p.spheres_radius_)[sphere_index] = radius;
	}

	void update_spheres(SurfaceParameters& p)
	{
		// auto start = std::chrono::high_resolution_clock::now();

		std::lock_guard<std::mutex> lock(p.mutex_);

		switch (p.update_method_)
		{
		case MEAN: {
			parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				update_sphere_mean(p, v);
				return true;
			});
			break;
		}
		case FIT: {
			parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				update_sphere_fit(p, v);
				return true;
			});
			break;
		}
		}

		compute_clusters(p);
		compute_spheres_error(p);

		// auto end = std::chrono::high_resolution_clock::now();
		// std::cout << "Update spheres: " << std::chrono::duration<Scalar>(end - start).count() << "s" << std::endl;

		if (!p.running_)
			update_render_data(p);
	}

	void split_sphere(SurfaceParameters& p, PVertex v)
	{
		uint32 v_index = index_of(*p.spheres_, v);

		std::lock_guard<std::mutex> lock(p.mutex_);

		// find the vertex of the cluster with maximal distance to the sphere

		const Vec3& pos = (*p.spheres_position_)[v_index];
		const Scalar r = (*p.spheres_radius_)[v_index];

		SVertex farthest_vertex;
		Scalar farthest_vertex_dist = 0.0;
		const std::vector<SVertex>& cluster = (*p.spheres_cluster_)[v_index];
		for (SVertex sv : cluster)
		{
			Scalar dist = (value<Vec3>(*p.surface_, p.surface_vertex_position_, sv) - pos).norm() - r;
			if (dist > farthest_vertex_dist)
			{
				farthest_vertex_dist = dist;
				farthest_vertex = sv;
			}
		}
		uint32 farthest_vertex_index = index_of(*p.surface_, farthest_vertex);

		// insert a new sphere for this vertex

		PVertex new_sphere = add_vertex(*p.spheres_);
		uint32 new_sphere_index = index_of(*p.spheres_, new_sphere);
		(*p.spheres_position_)[new_sphere_index] = (*p.medial_axis_position_)[farthest_vertex_index];
		(*p.spheres_radius_)[new_sphere_index] = (*p.medial_axis_radius_)[farthest_vertex_index];

		compute_clusters(p);
		compute_spheres_error(p);

		if (!p.running_)
			update_render_data(p);
	}

	void remove_sphere(SurfaceParameters& p, PVertex v)
	{
		std::lock_guard<std::mutex> lock(p.mutex_);

		const std::vector<SVertex>& cluster = value<std::vector<SVertex>>(*p.spheres_, p.spheres_cluster_, v);
		for (SVertex sv : cluster)
			value<PVertex>(*p.surface_, p.surface_vertex_cluster_, sv) = PVertex();
		remove_vertex(*p.spheres_, v);

		compute_clusters(p);
		compute_spheres_error(p);

		if (!p.running_)
			update_render_data(p);
	}

	PVertex max_error_sphere(SurfaceParameters& p)
	{
		// find sphere with maximal error according to the selected error attribute

		PVertex max_sphere;
		Scalar max = 0.0;
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			Scalar error = value<Scalar>(*p.spheres_, p.selected_spheres_error_, v);
			if (error > max)
			{
				max = error;
				max_sphere = v;
			}
			return true;
		});

		return max_sphere;
	}

	void compute_skeleton(SurfaceParameters& p)
	{
		clear(*p.skeleton_);

		auto spheres_skeleton_vertex_map =
			add_attribute<NMVertex, PVertex>(*p.spheres_, "__spheres_skeleton_vertex_map");

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv = add_vertex(*p.skeleton_);
			value<Vec3>(*p.skeleton_, p.skeleton_position_, nmv) = (*p.spheres_position_)[pv_index];
			(*spheres_skeleton_vertex_map)[pv_index] = nmv;
			return true;
		});

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv1 = (*spheres_skeleton_vertex_map)[pv_index];
			const std::set<PVertex>& neighbors = (*p.spheres_neighbor_clusters_)[pv_index];
			for (PVertex neighbor : neighbors)
			{
				NMVertex nmv2 = value<NMVertex>(*p.spheres_, spheres_skeleton_vertex_map, neighbor);
				std::vector<NMVertex> av = adjacent_vertices_through_edge(*p.skeleton_, nmv1);
				if (std::find(av.begin(), av.end(), nmv2) == av.end())
					add_edge(*p.skeleton_, nmv1, nmv2);
			}
			return true;
		});

		remove_attribute<PVertex>(*p.spheres_, spheres_skeleton_vertex_map);
	}

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));

		non_manifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() {
			SurfaceParameters& p = surface_parameters_[selected_surface_];
			update_render_data(p);
		});
	}

	void update_render_data(SurfaceParameters& p)
	{
		if (p.running_)
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

			update_spheres_color(p);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());

			compute_skeleton(p);
		}
		else
		{
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

			update_spheres_color(p);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());

			compute_skeleton(p);
		}

		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
	}

	void start_spheres_update(SurfaceParameters& p)
	{
		p.running_ = true;

		launch_thread([&]() {
			while (true)
			{
				update_spheres(p);
				if (p.slow_down_)
					std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
				else
					std::this_thread::yield();
				if (p.stopping_)
				{
					p.running_ = false;
					p.stopping_ = false;
					break;
				}
			}
		});

		app_.start_timer(100, [&]() -> bool { return !p.running_; });
	}

	void stop_spheres_update(SurfaceParameters& p)
	{
		p.stopping_ = true;
	}

	void key_press_event(View* view, int32 key_code) override
	{
		SurfaceParameters& p = surface_parameters_[selected_surface_];

		if (key_code == GLFW_KEY_I || key_code == GLFW_KEY_S || key_code == GLFW_KEY_D)
		{
			int32 x = view->mouse_x();
			int32 y = view->mouse_y();

			rendering::GLVec3d near = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near.x(), near.y(), near.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			Vec3 picked_sphere_center;
			foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				if (!picked_sphere_.is_valid())
				{
					picked_sphere_ = v;
					picked_sphere_center = value<Vec3>(*p.spheres_, p.spheres_position_, picked_sphere_);
					return true;
				}
				const Vec3& sp = value<Vec3>(*p.spheres_, p.spheres_position_, v);
				if (geometry::squared_distance_line_point(A, B, sp) <
					geometry::squared_distance_line_point(A, B, picked_sphere_center))
				{
					picked_sphere_ = v;
					picked_sphere_center = sp;
				}
				return true;
			});

			if (key_code == GLFW_KEY_S && picked_sphere_.is_valid())
				split_sphere(p, picked_sphere_);
			else if (key_code == GLFW_KEY_D && picked_sphere_.is_valid())
				remove_sphere(p, picked_sphere_);
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface",
							[&](SURFACE& s) { set_selected_surface(s); });

		if (selected_surface_)
		{
			SurfaceParameters& p = surface_parameters_[selected_surface_];

			imgui_combo_attribute<SVertex, Vec3>(
				*selected_surface_, p.surface_vertex_position_, "Position",
				[&](const std::shared_ptr<SAttribute<Vec3>>& attribute) { p.surface_vertex_position_ = attribute; });

			if (p.surface_vertex_position_ && !p.initialized_)
			{
				if (ImGui::Button("Init surface data"))
					init_surface_data(*selected_surface_);
			}

			if (p.initialized_)
			{
				ImGui::SliderFloat("Init min radius", &p.init_min_radius_, 0.001, 0.02);
				ImGui::SliderFloat("Init dilation factor", &p.init_dilation_factor_, 1.0, 4.0);
				if (ImGui::Button("Init spheres"))
					init_spheres(*selected_surface_);

				if (ImGui::SliderFloat("Samples min radius", &p.filter_radius_threshold_, 0.0001f, 1.0f, "%.4f",
									   ImGuiSliderFlags_Logarithmic))
					filter_medial_samples(p);
				if (ImGui::SliderFloat("Samples min angle", &p.filter_angle_threshold_, 0.0001f, M_PI, "%.4f"))
					filter_medial_samples(p);

				ImGui::RadioButton("Fit", (int*)&p.update_method_, FIT);
				ImGui::SameLine();
				ImGui::RadioButton("Mean", (int*)&p.update_method_, MEAN);

				ImGui::SliderFloat("Mean update curvature weight", &p.mean_update_curvature_weight_, 0.0f, 1.0f, "%.3f",
								   ImGuiSliderFlags_Logarithmic);

				ImGui::Checkbox("Auto split outside spheres", &p.auto_split_outside_spheres_);

				if (ImGui::Button("Update spheres"))
					update_spheres(p);

				ImGui::Checkbox("Slow down", &p.slow_down_);
				if (p.slow_down_)
					ImGui::SliderInt("Update rate", (int*)&p.update_rate_, 1, 1000);
				if (!p.running_)
				{
					if (ImGui::Button("Start spheres update"))
						start_spheres_update(p);
				}
				else
				{
					if (ImGui::Button("Stop spheres update"))
						stop_spheres_update(p);
				}

				ImGui::Separator();

				imgui_combo_attribute<PVertex, Scalar>(*p.spheres_, p.selected_spheres_error_, "Error measure",
													   [&](const std::shared_ptr<PAttribute<Scalar>>& attribute) {
														   set_selected_spheres_error(p, attribute.get());
													   });

				if (ImGui::Button("Split max error sphere"))
				{
					PVertex v = max_error_sphere(p);
					split_sphere(p, v);
				}

				ImGui::Separator();

				ImGui::Text("Total error: %f", p.total_error_);
				ImGui::Text("Min error: %f", p.min_error_);
				ImGui::Text("Max error: %f", p.max_error_);

				ImGui::Separator();

				ImGui::Text("Pick the sphere under the mouse with I and split it with S");
				if (picked_sphere_.is_valid())
				{
					ImGui::Text("Picked sphere:");
					const Vec3& sp = value<Vec3>(*p.spheres_, p.spheres_position_, picked_sphere_);
					ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
					ImGui::Text("Radius: %f", value<Scalar>(*p.spheres_, p.spheres_radius_, picked_sphere_));
				}
			}
		}
	}

private:
	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<POINTS>* points_provider_ = nullptr;
	MeshProvider<NONMANIFOLD>* non_manifold_provider_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	SURFACE* selected_surface_ = nullptr;
	PVertex picked_sphere_;

	std::array<std::mutex, 43> spheres_mutex_;

	std::shared_ptr<boost::synapse::connection> timer_connection_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SPHERE_FITTING_H_
