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

#ifndef CGOGN_MODULE_TUBULAR_MESH_H_
#define CGOGN_MODULE_TUBULAR_MESH_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/algos/distance.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/algos/filtering.h>
#include <cgogn/geometry/algos/hex_quality.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>

#include <cgogn/modeling/algos/graph_resampling.h>
#include <cgogn/modeling/algos/graph_to_hex.h>
#include <cgogn/modeling/algos/subdivision.h>

#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>

namespace cgogn
{

namespace ui
{

template <typename GRAPH, typename SURFACE, typename VOLUME>
class TubularMesh : public ViewModule
{
	template <typename T>
	using GraphAttribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	template <typename T>
	using VolumeAttribute = typename mesh_traits<VOLUME>::template Attribute<T>;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	using VolumeVertex = typename mesh_traits<VOLUME>::Vertex;
	using VolumeVertex2 = typename mesh_traits<VOLUME>::Vertex2;
	using VolumeEdge = typename mesh_traits<VOLUME>::Edge;
	using VolumeFace = typename mesh_traits<VOLUME>::Face;
	using VolumeVolume = typename mesh_traits<VOLUME>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;
	using Mat3 = geometry::Mat3;

public:
	TubularMesh(const App& app)
		: ViewModule(app, "TubularMesh"), graph_(nullptr), graph_vertex_position_(nullptr),
		  graph_vertex_radius_(nullptr), surface_(nullptr), surface_vertex_position_(nullptr), surface_bvh_(nullptr),
		  volume_vertex_position_(nullptr), volume_edge_target_length_(nullptr), volume_(nullptr)
	{
	}

	~TubularMesh()
	{
		if (transversal_faces_marker_)
			delete transversal_faces_marker_;
	}

protected:
	void init() override
	{
		graph_provider_ = static_cast<ui::MeshProvider<GRAPH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
		volume_provider_ = static_cast<ui::MeshProvider<VOLUME>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<VOLUME>::name} + ")"));
	}

public:
	void init_graph_radius_from_edge_length()
	{
		Scalar l = geometry::mean_edge_length(*graph_, graph_vertex_position_.get());
		graph_vertex_radius_->fill(l / 4.0);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
	}

	void extend_graph_extremities()
	{
		using SelectedFace = std::tuple<SurfaceFace, Vec3, Scalar>;
		CellCache<Graph> cache(*graph_);
		cache.template build<GraphVertex>();
		foreach_cell(cache, [&](GraphVertex v) -> bool {
			if (degree(*graph_, v) == 1)
			{
				std::vector<GraphVertex> av = adjacent_vertices_through_edge(*graph_, v);
				const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
				const Vec3& q = value<Vec3>(*graph_, graph_vertex_position_, av[0]);
				Vec3 dir = p - q;

				acc::Ray<Vec3> r{p, dir, 0, acc::inf};
				acc::BVHTree<uint32, Vec3>::Hit h;
				if (surface_bvh_->intersect(r, &h))
				{
					SurfaceFace f = surface_faces_[h.idx];
					std::vector<SurfaceVertex> vertices = incident_vertices(*surface_, f);
					Vec3 pos = h.bcoords[0] * value<Vec3>(*surface_, surface_vertex_position_, vertices[0]) +
							   h.bcoords[1] * value<Vec3>(*surface_, surface_vertex_position_, vertices[1]) +
							   h.bcoords[2] * value<Vec3>(*surface_, surface_vertex_position_, vertices[2]);
					GraphVertex nv = add_vertex(*graph_);
					connect_vertices(*graph_, v, nv);
					value<Vec3>(*graph_, graph_vertex_position_, nv) = p + 0.6 * (pos - p);
				}
			}
			return true;
		});

		graph_provider_->emit_connectivity_changed(graph_);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
	}

	void recenter_graph_from_surface()
	{
		parallel_foreach_cell(*graph_, [&](Graph::Vertex v) -> bool {
			for(uint32 i = 0; i < 10; ++i){
				const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
				Vec3 cp = surface_bvh_->closest_point(p);
				Vec3 cp_p = p - cp;
				value<Vec3>(*graph_, graph_vertex_position_, v) += 0.05 * cp_p;
			}

			const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
			Vec3 cp = surface_bvh_->closest_point(p);
			Vec3 cp_p = p - cp;
			value<Scalar>(*graph_, graph_vertex_radius_, v) = cp_p.norm();
			
			return true;
		});
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
		cgogn::io::export_CGR(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), "export.cgr");
	}

	void init_graph_radius_from_surface()
	{
		parallel_foreach_cell(*graph_, [&](Graph::Vertex v) -> bool {
			const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
			Vec3 cp = surface_bvh_->closest_point(p);
			value<Scalar>(*graph_, graph_vertex_radius_, v) = (cp - p).norm();
			return true;
		});
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
		cgogn::io::export_CGR(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), "export.cgr");
	}

	GRAPH* resample_graph(Scalar density)
	{
		static uint32 count = 0;
		GRAPH* resampled_graph = graph_provider_->add_mesh("resampled_" + std::to_string(count++));
		auto resampled_graph_vertex_position = add_attribute<Vec3, GraphVertex>(*resampled_graph, "position");
		auto resampled_graph_vertex_radius = add_attribute<Scalar, GraphVertex>(*resampled_graph, "radius");

		modeling::resample_graph(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), *resampled_graph,
								 resampled_graph_vertex_position.get(), resampled_graph_vertex_radius.get(), density);

		graph_provider_->emit_connectivity_changed(graph_);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());

		graph_provider_->emit_connectivity_changed(resampled_graph);
		graph_provider_->emit_attribute_changed(resampled_graph, resampled_graph_vertex_position.get());
		graph_provider_->emit_attribute_changed(resampled_graph, resampled_graph_vertex_radius.get());

		return resampled_graph;
	}

	VOLUME* build_hex_mesh()
	{
		// Scalar min_radius = std::numeric_limits<Scalar>::max();
		// for (Scalar r : *graph_vertex_radius_)
		// 	if (r < min_radius)
		// 		min_radius = r;

		// auto radius_copy = add_attribute<Scalar, GraphVertex>(*graph_, "radius_copy");
		// radius_copy->copy(graph_vertex_radius_.get());
		// graph_vertex_radius_->fill(min_radius);

		for (Scalar& r : *graph_vertex_radius_)
			r = r / 1.1;

		contact_surface_ = surface_provider_->add_mesh("contact");
		volume_ = volume_provider_->add_mesh("hex");

		hex_building_attributes_ = modeling::graph_to_hex(*graph_, *contact_surface_, *volume_);

		// if (!transversal_faces_marker_)
		// {
		// 	transversal_faces_marker_ = new CellMarker<VOLUME, VolumeFace>(*volume_);
		// 	modeling::mark_tranversal_faces(*volume_, *contact_surface_, std::get<1>(hex_building_attributes_),
		// 									*transversal_faces_marker_);
		// }

		// graph_vertex_radius_->swap(radius_copy.get());
		// remove_attribute<GraphVertex>(*graph_, radius_copy);

		if (check_integrity(*volume_))
			std::cout << "Volume mesh OK!" << std::endl;

		surface_provider_->emit_connectivity_changed(contact_surface_);
		volume_provider_->emit_connectivity_changed(volume_);

		volume_vertex_position_ = get_attribute<Vec3, VolumeVertex>(*volume_, "position");
		volume_edge_target_length_ = add_attribute<Scalar, VolumeEdge>(*volume_, "target_length");
		volume_provider_->set_mesh_bb_vertex_position(volume_, volume_vertex_position_);

		return volume_;
	}

	void project_on_surface()
	{
		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");
		auto volume_skin_vertex_normal = add_attribute<Vec3, SurfaceVertex>(volume_skin, "normal");
		auto volume_skin_vertex_volume_vertex = add_attribute<VolumeVertex, SurfaceVertex>(volume_skin, "hex_vertex");

		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get(), volume_skin_vertex_volume_vertex.get());
		geometry::apply_ear_triangulation(volume_skin, volume_skin_vertex_position.get());

		// auto normal_filtered = add_attribute<Vec3, SurfaceVertex>(volume_skin, "normal_filtered");

		geometry::compute_normal(volume_skin, volume_skin_vertex_position.get(), volume_skin_vertex_normal.get());
		// geometry::filter_average<Vec3>(volume_skin, volume_skin_vertex_normal.get(), normal_filtered.get());
		// geometry::filter_average<Vec3>(volume_skin, normal_filtered.get(), volume_skin_vertex_normal.get());
		// volume_skin_vertex_normal->swap(normal_filtered.get());

		// for (Vec3& n : *volume_skin_vertex_normal)
		// 	n.normalize();

		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			const Vec3& p = value<Vec3>(volume_skin, volume_skin_vertex_position, v);
			const Vec3& n = value<Vec3>(volume_skin, volume_skin_vertex_normal, v);
			Vec3 pos;

			Scalar local_size = 0.0;
			uint32 nb_neigh = 0;
			foreach_adjacent_vertex_through_edge(volume_skin, v, [&](SurfaceVertex av) -> bool {
				local_size += (value<Vec3>(volume_skin, volume_skin_vertex_position, av) - p).norm();
				++nb_neigh;
				return true;
			});
			local_size /= nb_neigh;

			acc::Ray<Vec3> r{p, n, 0, 1.5 * local_size};
			acc::BVHTree<uint32, Vec3>::Hit h;
			if (surface_bvh_->intersect(r, &h))
			{
				SurfaceFace f = surface_faces_[h.idx];
				std::vector<SurfaceVertex> vertices = incident_vertices(*surface_, f);
				pos = h.bcoords[0] * value<Vec3>(*surface_, surface_vertex_position_, vertices[0]) +
					  h.bcoords[1] * value<Vec3>(*surface_, surface_vertex_position_, vertices[1]) +
					  h.bcoords[2] * value<Vec3>(*surface_, surface_vertex_position_, vertices[2]);
			}
			else
				pos = surface_bvh_->closest_point(p);
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v)) = pos;
			return true;
		});

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void add_volume_padding()
	{
		modeling::padding(*volume_);

		volume_provider_->emit_connectivity_changed(volume_);
		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());

		refresh_edge_target_length_ = true;
	}

	// void subdivide_volume_length_wise()
	// {
	// 	modeling::subdivide_length_wise(*volume_, std::get<2>(hex_building_attributes_), *transversal_faces_marker_,
	// 									*graph_, std::get<0>(hex_building_attributes_));

	// 	volume_provider_->emit_connectivity_changed(volume_);
	// 	volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());

	// 	refresh_edge_target_length_ = true;
	// }

	// void subdivide_volume_width_wise()
	// {
	// 	modeling::subdivide_width_wise(*volume_, std::get<2>(hex_building_attributes_), *transversal_faces_marker_,
	// 								   *graph_, std::get<0>(hex_building_attributes_));

	// 	graph_provider_->emit_connectivity_changed(graph_);
	// 	graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
	// 	graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());

	// 	volume_provider_->emit_connectivity_changed(volume_);
	// 	volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());

	// 	refresh_edge_target_length_ = true;
	// }

	void subdivide_volume()
	{
		modeling::primal_cut_all_volumes(
			*volume_,
			[&](VolumeVertex v) {
				std::vector<VolumeVertex> av = adjacent_vertices_through_edge(*volume_, v);
				cgogn::value<Vec3>(*volume_, volume_vertex_position_, v) =
					0.5 * (cgogn::value<Vec3>(*volume_, volume_vertex_position_, av[0]) +
						   cgogn::value<Vec3>(*volume_, volume_vertex_position_, av[1]));
			},
			[&](VolumeVertex v) {
				Vec3 center;
				center.setZero();
				uint32 count = 0;
				foreach_adjacent_vertex_through_edge(*volume_, v, [&](VolumeVertex av) -> bool {
					center += cgogn::value<Vec3>(*volume_, volume_vertex_position_, av);
					++count;
					return true;
				});
				center /= Scalar(count);
				cgogn::value<Vec3>(*volume_, volume_vertex_position_, v) = center;
			},
			[&](VolumeVertex v) {
				Vec3 center;
				center.setZero();
				uint32 count = 0;
				foreach_adjacent_vertex_through_edge(*volume_, v, [&](VolumeVertex av) -> bool {
					center += cgogn::value<Vec3>(*volume_, volume_vertex_position_, av);
					++count;
					return true;
				});
				center /= Scalar(count);
				cgogn::value<Vec3>(*volume_, volume_vertex_position_, v) = center;
			});

		volume_provider_->emit_connectivity_changed(volume_);
		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());

		refresh_edge_target_length_ = true;
	}

	void smooth_volume_surface()
	{
		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");
		auto volume_skin_vertex_laplacian = add_attribute<Vec3, SurfaceVertex>(volume_skin, "laplacian");
		auto volume_skin_vertex_volume_vertex = add_attribute<VolumeVertex, SurfaceVertex>(volume_skin, "hex_vertex");

		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get(), volume_skin_vertex_volume_vertex.get());
		geometry::apply_ear_triangulation(volume_skin, volume_skin_vertex_position.get());

		geometry::compute_laplacian(volume_skin, volume_skin_vertex_position.get(), volume_skin_vertex_laplacian.get());
		parallel_foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			value<Vec3>(volume_skin, volume_skin_vertex_position, v) +=
				0.15 * value<Vec3>(volume_skin, volume_skin_vertex_laplacian, v);
			return true;
		});

		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			Vec3 cp = surface_bvh_->closest_point(value<Vec3>(volume_skin, volume_skin_vertex_position, v));
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v)) = cp;
			return true;
		});

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void regularize_surface_vertices(Scalar fit_to_data = 5.0)
	{
		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");
		auto volume_skin_vertex_laplacian = add_attribute<Vec3, SurfaceVertex>(volume_skin, "laplacian");
		auto volume_skin_vertex_volume_vertex = add_attribute<VolumeVertex, SurfaceVertex>(volume_skin, "hex_vertex");

		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get(), volume_skin_vertex_volume_vertex.get());
		// geometry::apply_ear_triangulation(volume_skin, volume_skin_vertex_position.get());

		auto vertex_index = add_attribute<uint32, SurfaceVertex>(volume_skin, "__vertex_index");

		uint32 nb_vertices = 0;
		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			value<uint32>(volume_skin, vertex_index, v) = nb_vertices++;
			return true;
		});

		// Scalar fit_to_data = geometry::mean_edge_length(volume_skin, volume_skin_vertex_position.get()) * 5.0;

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * nb_vertices, nb_vertices);
		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_vertices * 10);
		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			uint32 vidx = value<uint32>(volume_skin, vertex_index, v);
			// const Vec3& p = value<Vec3>(volume_skin, volume_skin_vertex_position, v);
			auto vertices = adjacent_vertices_through_edge(volume_skin, v);
			// Scalar d = 0;
			auto d = vertices.size();
			for (SurfaceVertex av : vertices)
			{
				uint32 avidx = value<uint32>(volume_skin, vertex_index, av);
				// Scalar dist = (value<Vec3>(volume_skin, volume_skin_vertex_position, av) - p).norm();
				// Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), dist));
				// d += dist;
				Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
			}
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(d)));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), fit_to_data));
			return true;
		});
		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(2 * nb_vertices, 3);

		parallel_foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			uint32 vidx = value<uint32>(volume_skin, vertex_index, v);
			b(vidx, 0) = 0;
			b(vidx, 1) = 0;
			b(vidx, 2) = 0;
			const Vec3& pos = value<Vec3>(volume_skin, volume_skin_vertex_position, v);
			b(nb_vertices + vidx, 0) = fit_to_data * pos[0];
			b(nb_vertices + vidx, 1) = fit_to_data * pos[1];
			b(nb_vertices + vidx, 2) = fit_to_data * pos[2];
			x(vidx, 0) = pos[0];
			x(vidx, 1) = pos[1];
			x(vidx, 2) = pos[2];
			return true;
		});

		x = solver.solveWithGuess(b, x);

		parallel_foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			uint32 vidx = value<uint32>(volume_skin, vertex_index, v);
			Vec3& pos = value<Vec3>(volume_skin, volume_skin_vertex_position, v);
			pos[0] = x(vidx, 0);
			pos[1] = x(vidx, 1);
			pos[2] = x(vidx, 2);
			return true;
		});

		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			Vec3 cp = surface_bvh_->closest_point(value<Vec3>(volume_skin, volume_skin_vertex_position, v));
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v)) = cp;
			return true;
		});

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void relocate_interior_vertices()
	{
		auto vertex_index = add_attribute<uint32, VolumeVertex>(*volume_, "__vertex_index");
		auto vertex_laplacian = add_attribute<Vec3, VolumeVertex>(*volume_, "laplacian");

		uint32 nb_vertices = 0;
		foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			value<uint32>(*volume_, vertex_index, v) = nb_vertices++;
			return true;
		});

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(nb_vertices, nb_vertices);
		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_vertices * 10);
		// foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
		// 	uint32 vidx = value<uint32>(*volume_, vertex_index, v);
		// 	auto vertices = adjacent_vertices_through_edge(*volume_, v);
		// 	auto d = vertices.size();
		// 	for (VolumeVertex av : vertices)
		// 	{
		// 		uint32 avidx = value<uint32>(*volume_, vertex_index, av);
		// 		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
		// 	}
		// 	Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(d)));
		// 	return true;
		// });
		foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			auto vertices = incident_vertices(*volume_, e);
			uint32 vidx1 = value<uint32>(*volume_, vertex_index, vertices[0]);
			uint32 vidx2 = value<uint32>(*volume_, vertex_index, vertices[1]);
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), 1));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), 1));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -1));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -1));
			return true;
		});
		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		// // set constrained vertices
		// foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
		// 	if (is_incident_to_boundary(*volume_, v))
		// 	{
		// 		int idx = int(value<uint32>(*volume_, vertex_index, v));
		// 		A.prune([&](int i, int, Scalar) { return i != idx; });
		// 		A.coeffRef(idx, idx) = 1.0;
		// 	}
		// 	return true;
		// });
		// A.makeCompressed();

		// Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);
		// Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(LAPL);

		// Eigen::MatrixXd x(nb_vertices, 3);
		// Eigen::MatrixXd b(nb_vertices, 3);

		// parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
		// 	uint32 vidx = value<uint32>(*volume_, vertex_index, v);
		// 	const Vec3& pos = value<Vec3>(*volume_, volume_vertex_position_, v);
		// 	if (is_incident_to_boundary(*volume_, v))
		// 	{
		// 		b.coeffRef(vidx, 0) = pos[0];
		// 		b.coeffRef(vidx, 1) = pos[1];
		// 		b.coeffRef(vidx, 2) = pos[2];
		// 	}
		// 	else
		// 	{
		// 		b.coeffRef(vidx, 0) = 0;
		// 		b.coeffRef(vidx, 1) = 0;
		// 		b.coeffRef(vidx, 2) = 0;
		// 	}
		// 	x(vidx, 0) = pos[0];
		// 	x(vidx, 1) = pos[1];
		// 	x(vidx, 2) = pos[2];
		// 	return true;
		// });

		// x = solver.solveWithGuess(b, x);

		Eigen::MatrixXd vpos(nb_vertices, 3);
		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			const Vec3& pv = value<Vec3>(*volume_, volume_vertex_position_, v);
			uint32 vidx = value<uint32>(*volume_, vertex_index, v);
			vpos(vidx, 0) = pv[0];
			vpos(vidx, 1) = pv[1];
			vpos(vidx, 2) = pv[2];
			return true;
		});

		Eigen::MatrixXd lapl(nb_vertices, 3);
		lapl = A * vpos;

		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			if (!is_incident_to_boundary(*volume_, v))
			{
				uint32 vidx = value<uint32>(*volume_, vertex_index, v);
				Vec3 l;
				l[0] = lapl(vidx, 0);
				l[1] = lapl(vidx, 1);
				l[2] = lapl(vidx, 2);
				value<Vec3>(*volume_, volume_vertex_position_, v) += 0.1 * l;
			}
			return true;
		});

		// parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
		// 	uint32 vidx = value<uint32>(*volume_, vertex_index, v);
		// 	Vec3& pos = value<Vec3>(*volume_, volume_vertex_position_, v);
		// 	pos[0] = x(vidx, 0);
		// 	pos[1] = x(vidx, 1);
		// 	pos[2] = x(vidx, 2);
		// 	return true;
		// });

		remove_attribute<VolumeVertex>(*volume_, vertex_index);

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void compute_target_edge_length()
	{
		foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			// auto vertices = incident_vertices(*volume_, e);

			std::vector<VolumeEdge> parallel_edges;
			parallel_edges.reserve(16);
			// std::vector<VolumeEdge> perpendicular_edges;
			// perpendicular_edges.reserve(16);

			Dart ed = e.dart;
			parallel_edges.push_back(VolumeEdge(ed)); // the edge itself
			do
			{
				Dart vd = phi<211>(*volume_, ed);
				parallel_edges.push_back(VolumeEdge(vd));
				if (!is_boundary(*volume_, ed))
				{
					vd = phi<211>(*volume_, vd);
					parallel_edges.push_back(VolumeEdge(vd));
				}
				else
					parallel_edges.push_back(VolumeEdge(ed)); // edge is on the boundary -> count twice
				// perpendicular_edges.push_back(VolumeEdge(phi1(*volume_, ed)));
				// perpendicular_edges.push_back(VolumeEdge(phi_1(*volume_, ed)));
				ed = phi<32>(*volume_, ed);
			} while (ed != e.dart);

			Scalar parallel_edges_mean_length = 0.0;
			for (VolumeEdge pe : parallel_edges)
				parallel_edges_mean_length += geometry::length(*volume_, pe, volume_vertex_position_.get());
			parallel_edges_mean_length /= parallel_edges.size();

			// Scalar perpendicular_edges_mean_length = 0.0;
			// for (VolumeEdge pe : perpendicular_edges)
			// 	perpendicular_edges_mean_length += geometry::length(*volume_, pe, volume_vertex_position_.get());
			// perpendicular_edges_mean_length /= perpendicular_edges.size();

			// Scalar target_length = edge_length;
			// if (is_incident_to_boundary(*volume_, vertices[0]) && is_incident_to_boundary(*volume_, vertices[1]))
			// 	target_length = 0.5 * edge_length + 0.5 * parallel_edges_mean_length;
			// else
			// 	target_length = parallel_edges_mean_length;

			value<Scalar>(*volume_, volume_edge_target_length_, e) = parallel_edges_mean_length;
			// 0.5 * (parallel_edges_mean_length + perpendicular_edges_mean_length);

			return true;
		});
		refresh_edge_target_length_ = false;
	}

	void optimize_volume_vertices(Scalar fit_to_data = 50.0)
	{
		if (refresh_edge_target_length_)
			compute_target_edge_length();

		auto vertex_index = add_attribute<uint32, VolumeVertex>(*volume_, "__vertex_index");

		uint32 nb_vertices = 0;
		uint32 nb_boundary_vertices = 0;
		foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			value<uint32>(*volume_, vertex_index, v) = nb_vertices++;
			if (is_incident_to_boundary(*volume_, v))
				++nb_boundary_vertices;
			return true;
		});

		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");
		auto volume_skin_vertex_normal = add_attribute<Vec3, SurfaceVertex>(volume_skin, "normal");
		auto volume_skin_vertex_volume_vertex = add_attribute<VolumeVertex, SurfaceVertex>(volume_skin, "hex_vertex");
		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get(), volume_skin_vertex_volume_vertex.get());
		geometry::apply_ear_triangulation(volume_skin, volume_skin_vertex_position.get());
		geometry::compute_normal(volume_skin, volume_skin_vertex_position.get(), volume_skin_vertex_normal.get());
		auto normal_filtered = add_attribute<Vec3, SurfaceVertex>(volume_skin, "normal_filtered");
		geometry::filter_average<Vec3>(volume_skin, volume_skin_vertex_normal.get(), normal_filtered.get());
		geometry::filter_average<Vec3>(volume_skin, normal_filtered.get(), volume_skin_vertex_normal.get());

		uint32 nb_oriented_edges = nb_cells<VolumeEdge>(*volume_) * 2;

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(nb_oriented_edges + nb_boundary_vertices, nb_vertices);

		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_oriented_edges * 2);

		uint32 oriented_edge_idx = 0;
		foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			auto vertices = incident_vertices(*volume_, e);
			uint32 vidx1 = value<uint32>(*volume_, vertex_index, vertices[0]);
			uint32 vidx2 = value<uint32>(*volume_, vertex_index, vertices[1]);

			Scalar target_length = value<Scalar>(*volume_, volume_edge_target_length_, e);

			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(oriented_edge_idx), int(vidx1), -1 / target_length));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(oriented_edge_idx), int(vidx2), 1 / target_length));

			++oriented_edge_idx;

			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(oriented_edge_idx), int(vidx1), 1 / target_length));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(oriented_edge_idx), int(vidx2), -1 / target_length));

			++oriented_edge_idx;

			return true;
		});

		// set constrained vertices
		uint32 boundary_vertex_idx = 0;
		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_, vertex_index,
										value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v));
			Acoeffs.push_back(
				Eigen::Triplet<Scalar>(int(oriented_edge_idx + boundary_vertex_idx), int(vidx), fit_to_data));
			++boundary_vertex_idx;
			return true;
		});

		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(nb_oriented_edges + nb_boundary_vertices, 3);

		oriented_edge_idx = 0;
		foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			auto vertices = incident_vertices(*volume_, e);

			// uint32 vidx1 = value<uint32>(*volume_, vertex_index, vertices[0]);
			const Vec3& pos1 = value<Vec3>(*volume_, volume_vertex_position_, vertices[0]);
			// uint32 vidx2 = value<uint32>(*volume_, vertex_index, vertices[1]);
			const Vec3& pos2 = value<Vec3>(*volume_, volume_vertex_position_, vertices[1]);

			Vec3 edge1 = (pos2 - pos1).normalized();
			Vec3 edge2 = (pos1 - pos2).normalized();

			Vec3 target_n1(0, 0, 0);
			Vec3 target_n2(0, 0, 0);

			// if (!(is_incident_to_boundary(*volume_, vertices[0]) && is_incident_to_boundary(*volume_, vertices[1])))
			{
				Dart d = e.dart;
				do
				{
					if (!is_boundary(*volume_, d))
					{
						Vec3 n =
							geometry::normal(
								pos1,
								value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi<211>(*volume_, d))),
								value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi_1(*volume_, d))))
								.normalized();
						if (edge1.dot(n) > 0)
							target_n1 += n;
					}
					else
						target_n1 += edge1;
					d = phi<32>(*volume_, d);
				} while (d != e.dart);
				target_n1.normalize();

				d = phi2(*volume_, e.dart);
				do
				{
					if (!is_boundary(*volume_, d))
					{
						Vec3 n =
							geometry::normal(
								pos2,
								value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi<211>(*volume_, d))),
								value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi_1(*volume_, d))))
								.normalized();
						if (edge2.dot(n) > 0)
							target_n2 += n;
					}
					else
						target_n2 += edge2;
					d = phi<32>(*volume_, d);
				} while (d != phi2(*volume_, e.dart));
				target_n2.normalize();
			}

			b.coeffRef(oriented_edge_idx, 0) = target_n1[0];
			b.coeffRef(oriented_edge_idx, 1) = target_n1[1];
			b.coeffRef(oriented_edge_idx, 2) = target_n1[2];

			++oriented_edge_idx;

			b.coeffRef(oriented_edge_idx, 0) = target_n2[0];
			b.coeffRef(oriented_edge_idx, 1) = target_n2[1];
			b.coeffRef(oriented_edge_idx, 2) = target_n2[2];

			++oriented_edge_idx;

			return true;
		});

		boundary_vertex_idx = 0;
		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			const Vec3& p = value<Vec3>(volume_skin, volume_skin_vertex_position, v);
			const Vec3& n = value<Vec3>(volume_skin, volume_skin_vertex_normal, v);
			Vec3 pos;

			Scalar local_size = 0.0;
			uint32 nb_neigh = 0;
			foreach_adjacent_vertex_through_edge(volume_skin, v, [&](SurfaceVertex av) -> bool {
				local_size += (value<Vec3>(volume_skin, volume_skin_vertex_position, av) - p).norm();
				++nb_neigh;
				return true;
			});
			local_size /= nb_neigh;

			acc::Ray<Vec3> r{p, n, 0, 1.5 * local_size};
			acc::BVHTree<uint32, Vec3>::Hit h;
			if (surface_bvh_->intersect(r, &h))
			{
				SurfaceFace f = surface_faces_[h.idx];
				std::vector<SurfaceVertex> vertices = incident_vertices(*surface_, f);
				pos = h.bcoords[0] * value<Vec3>(*surface_, surface_vertex_position_, vertices[0]) +
					  h.bcoords[1] * value<Vec3>(*surface_, surface_vertex_position_, vertices[1]) +
					  h.bcoords[2] * value<Vec3>(*surface_, surface_vertex_position_, vertices[2]);
			}
			else
				pos = surface_bvh_->closest_point(p);
			b.coeffRef(oriented_edge_idx + boundary_vertex_idx, 0) = fit_to_data * pos[0];
			b.coeffRef(oriented_edge_idx + boundary_vertex_idx, 1) = fit_to_data * pos[1];
			b.coeffRef(oriented_edge_idx + boundary_vertex_idx, 2) = fit_to_data * pos[2];
			++boundary_vertex_idx;
			return true;
		});

		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_, vertex_index, v);
			const Vec3& pos = value<Vec3>(*volume_, volume_vertex_position_, v);
			x(vidx, 0) = pos[0];
			x(vidx, 1) = pos[1];
			x(vidx, 2) = pos[2];
			return true;
		});

		x = solver.solveWithGuess(b, x);

		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_, vertex_index, v);
			Vec3& pos = value<Vec3>(*volume_, volume_vertex_position_, v);
			pos[0] = x(vidx, 0);
			pos[1] = x(vidx, 1);
			pos[2] = x(vidx, 2);
			return true;
		});

		remove_attribute<VolumeVertex>(*volume_, vertex_index);

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void compute_volumes_quality()
	{
		auto corner_frame = add_attribute<Mat3, VolumeVertex2>(*volume_, "__corner_frame");
		auto hex_frame = add_attribute<Mat3, VolumeVolume>(*volume_, "__hex_frame");

		auto scaled_jacobian = get_attribute<Scalar, VolumeVolume>(*volume_, "scaled_jacobian");
		if (!scaled_jacobian)
			scaled_jacobian = add_attribute<Scalar, VolumeVolume>(*volume_, "scaled_jacobian");

		auto jacobian = get_attribute<Scalar, VolumeVolume>(*volume_, "jacobian");
		if (!jacobian)
			jacobian = add_attribute<Scalar, VolumeVolume>(*volume_, "jacobian");

		auto max_froebnius = get_attribute<Scalar, VolumeVolume>(*volume_, "max_froebnius");
		if (!max_froebnius)
			max_froebnius = add_attribute<Scalar, VolumeVolume>(*volume_, "max_froebnius");

		auto mean_froebnius = get_attribute<Scalar, VolumeVolume>(*volume_, "mean_froebnius");
		if (!mean_froebnius)
			mean_froebnius = add_attribute<Scalar, VolumeVolume>(*volume_, "mean_froebnius");

		geometry::compute_hex_frame(*volume_, volume_vertex_position_.get(), corner_frame.get(), hex_frame.get());
		geometry::compute_scaled_jacobian(*volume_, corner_frame.get(), hex_frame.get(), scaled_jacobian.get());
		geometry::compute_jacobian(*volume_, corner_frame.get(), hex_frame.get(), jacobian.get());
		geometry::compute_maximum_aspect_frobenius(*volume_, corner_frame.get(), max_froebnius.get());
		geometry::compute_mean_aspect_frobenius(*volume_, corner_frame.get(), mean_froebnius.get());

		Scalar mean_scaled_jacobian = 0.0;
		for (Scalar s : *scaled_jacobian)
			mean_scaled_jacobian += s;
		mean_scaled_jacobian /= scaled_jacobian->size();
		std::cout << "mean scaled jacobian = " << mean_scaled_jacobian << std::endl;

		volume_provider_->emit_attribute_changed(volume_, scaled_jacobian.get());
		volume_provider_->emit_attribute_changed(volume_, jacobian.get());
		volume_provider_->emit_attribute_changed(volume_, max_froebnius.get());
		volume_provider_->emit_attribute_changed(volume_, mean_froebnius.get());

		remove_attribute<VolumeVolume>(*volume_, hex_frame.get());
		remove_attribute<VolumeVertex2>(*volume_, corner_frame.get());
	}

	void export_subdivided_skin()
	{
		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");

		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get());
		// modeling::catmull_clark_approx(volume_skin, volume_skin_vertex_position.get(), 2);
		geometry::apply_ear_triangulation(volume_skin, volume_skin_vertex_position.get());
		surface_provider_->save_surface_to_file(volume_skin, volume_skin_vertex_position.get(), "off", "surface");
	}

	void set_current_graph(GRAPH* g)
	{
		graph_ = g;
		graph_vertex_position_ = nullptr;
		graph_vertex_radius_ = nullptr;
	}

	void set_current_surface(SURFACE* s)
	{
		surface_ = s;
		surface_vertex_position_ = nullptr;
	}

	void set_current_volume(VOLUME* v)
	{
		volume_ = v;
		volume_vertex_position_ = get_attribute<Vec3, VolumeVertex>(*volume_, "position");
	}

	void set_current_graph_vertex_position(const std::shared_ptr<GraphAttribute<Vec3>>& attribute)
	{
		if (graph_)
			graph_vertex_position_ = attribute;
	}

	void set_current_graph_vertex_radius(const std::shared_ptr<GraphAttribute<Scalar>>& attribute)
	{
		if (graph_)
			graph_vertex_radius_ = attribute;
	}

	void set_current_surface_vertex_position(const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute)
	{
		if (surface_)
		{
			surface_vertex_position_ = attribute;

			if (surface_bvh_)
				delete surface_bvh_;

			uint32 nb_vertices = surface_provider_->mesh_data(surface_)->template nb_cells<SurfaceVertex>();
			uint32 nb_faces = surface_provider_->mesh_data(surface_)->template nb_cells<SurfaceFace>();

			auto vertex_index = add_attribute<uint32, SurfaceVertex>(*surface_, "__vertex_index");

			std::vector<Vec3> vertex_position;
			vertex_position.reserve(nb_vertices);
			uint32 idx = 0;
			foreach_cell(*surface_, [&](SurfaceVertex v) -> bool {
				value<uint32>(*surface_, vertex_index, v) = idx++;
				vertex_position.push_back(value<Vec3>(*surface_, surface_vertex_position_, v));
				return true;
			});

			surface_faces_.clear();
			surface_faces_.reserve(nb_faces);
			std::vector<uint32> face_vertex_indices;
			face_vertex_indices.reserve(nb_faces * 3);
			foreach_cell(*surface_, [&](SurfaceFace f) -> bool {
				surface_faces_.push_back(f);
				foreach_incident_vertex(*surface_, f, [&](SurfaceVertex v) -> bool {
					face_vertex_indices.push_back(value<uint32>(*surface_, vertex_index, v));
					return true;
				});
				return true;
			});

			surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position);
		}
	}

protected:
	void interface() override
	{
		ImGui::TextUnformatted("Graph");
		imgui_mesh_selector(graph_provider_, graph_, "Graph", [&](GRAPH* g) { set_current_graph(g); });

		if (graph_)
		{
			imgui_combo_attribute<GraphVertex, Vec3>(*graph_, graph_vertex_position_, "Position##graph",
													 [&](const std::shared_ptr<GraphAttribute<Vec3>>& attribute) {
														 set_current_graph_vertex_position(attribute);
													 });
			imgui_combo_attribute<GraphVertex, Scalar>(*graph_, graph_vertex_radius_, "Radius##graph",
													   [&](const std::shared_ptr<GraphAttribute<Scalar>>& attribute) {
														   set_current_graph_vertex_radius(attribute);
													   });

			if (ImGui::Button("Create radius attribute"))
				graph_vertex_radius_ = add_attribute<Scalar, GraphVertex>(*graph_, "radius");
		}

		ImGui::Separator();
		ImGui::TextUnformatted("Surface");
		imgui_mesh_selector(surface_provider_, surface_, "Surface", [&](SURFACE* s) { set_current_surface(s); });

		if (surface_)
		{
			imgui_combo_attribute<SurfaceVertex, Vec3>(*surface_, surface_vertex_position_, "Position##surface",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   set_current_surface_vertex_position(attribute);
													   });
		}

		ImGui::Separator();
		ImGui::TextUnformatted("Operations");
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_)
		{
			if (ImGui::Button("Init radius from edge length"))
				init_graph_radius_from_edge_length();
		}
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_ && surface_ && surface_vertex_position_)
		{
			if (ImGui::Button("Recenter graph from surface"))
				recenter_graph_from_surface();
			if (ImGui::Button("Init radius from surface"))
				init_graph_radius_from_surface();
			if (ImGui::Button("Extend graph extremities"))
				extend_graph_extremities();
		}
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_)
		{
			static float graph_resample_density = 0.5f;
			ImGui::SliderFloat("Graph resampling density", &graph_resample_density, 0.0, 1.0);
			if (ImGui::Button("Resample graph"))
				resample_graph(graph_resample_density);
			if (ImGui::Button("Build hex mesh"))
				build_hex_mesh();
		}
		if (volume_)
		{
			if (ImGui::Button("Export subdivided skin"))
				export_subdivided_skin();
			if (ImGui::Button("Add volume padding"))
				add_volume_padding();
			// if (ImGui::Button("Subdivide length wise"))
			// 	subdivide_volume_length_wise();
			// if (ImGui::Button("Subdivide width wise"))
			// 	subdivide_volume_width_wise();
			if (ImGui::Button("Subdivide volume"))
				subdivide_volume();
			if (ImGui::Button("Project on surface"))
				project_on_surface();
			static float regularize_fit_to_data = 5.0f;
			ImGui::SliderFloat("Regularize surface - Fit to data", &regularize_fit_to_data, 0.0, 20.0);
			if (ImGui::Button("Regularize surface vertices"))
				regularize_surface_vertices(regularize_fit_to_data);
			if (ImGui::Button("Relocate interior vertices"))
				relocate_interior_vertices();
			static float optimize_fit_to_surface = 20.0f;
			ImGui::SliderFloat("Optimize volume - Fit to surface", &optimize_fit_to_surface, 0.0, 200.0);
			ImGui::Checkbox("Refresh edge target length", &refresh_edge_target_length_);
			if (ImGui::Button("Optimize volume vertices"))
				optimize_volume_vertices(optimize_fit_to_surface);
			if (ImGui::Button("Compute volumes quality"))
				compute_volumes_quality();
		}
	}

private:
	MeshProvider<GRAPH>* graph_provider_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<VOLUME>* volume_provider_;

	GRAPH* graph_;
	std::shared_ptr<GraphAttribute<Vec3>> graph_vertex_position_;
	std::shared_ptr<GraphAttribute<Scalar>> graph_vertex_radius_;

	SURFACE* surface_;
	std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_position_;
	acc::BVHTree<uint32, Vec3>* surface_bvh_;
	std::vector<SurfaceFace> surface_faces_;

	SURFACE* contact_surface_;
	SURFACE* contact_surface_2;

	VOLUME* volume_;
	std::shared_ptr<VolumeAttribute<Vec3>> volume_vertex_position_;
	std::shared_ptr<VolumeAttribute<Scalar>> volume_edge_target_length_;
	bool refresh_edge_target_length_ = true;
	CellMarker<VOLUME, VolumeFace>* transversal_faces_marker_ = nullptr;

	std::tuple<modeling::GAttributes, modeling::M2Attributes, modeling::M3Attributes> hex_building_attributes_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_TUBULAR_MESH_H_
