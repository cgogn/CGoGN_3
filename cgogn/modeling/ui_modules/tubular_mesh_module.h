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

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/distance.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/algos/filtering.h>
#include <cgogn/geometry/algos/hex_quality.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/registration.h>

#include <cgogn/modeling/algos/graph_resampling.h>
#include <cgogn/modeling/algos/graph_to_hex.h>
#include <cgogn/modeling/algos/incidenceGraph_to_hex.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/algos/volume_utils.h>

#include <cgogn/io/utils.h>

#include <Eigen/Sparse>
#include <boost/synapse/connect.hpp>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

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
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

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
	TubularMesh(const App& app) : ViewModule(app, "TubularMesh")
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

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() { animate_volume(); });
	}

	struct BVH_Hit
	{
		bool hit = false;
		SurfaceFace face;
		Vec3 bcoords;
		Scalar dist;
		Vec3 pos;
	};

	BVH_Hit intersect_bvh(const acc::Ray<Vec3>& r)
	{
		acc::BVHTree<uint32, Vec3>::Hit h;
		if (surface_bvh_->intersect(r, &h))
		{
			SurfaceFace f = surface_faces_[h.idx];
			std::vector<SurfaceVertex> vertices = incident_vertices(*surface_, f);
			Vec3 p = h.bcoords[0] * value<Vec3>(*surface_, surface_vertex_position_, vertices[0]) +
					 h.bcoords[1] * value<Vec3>(*surface_, surface_vertex_position_, vertices[1]) +
					 h.bcoords[2] * value<Vec3>(*surface_, surface_vertex_position_, vertices[2]);
			return {true, f, {h.bcoords[0], h.bcoords[1], h.bcoords[2]}, h.t, p};
		}
		else
			return BVH_Hit();
	}

public:
	void init_graph_radius_from_edge_length()
	{
		Scalar l = geometry::mean_edge_length(*graph_, graph_vertex_position_.get());
		graph_vertex_radius_->fill(l / 4.0);
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_radius_.get());
	}

	void extend_graph_extremities()
	{
		using SelectedFace = std::tuple<SurfaceFace, Vec3, Scalar>;
		CellCache<GRAPH> cache(*graph_);
		cache.template build<GraphVertex>();
		foreach_cell(cache, [&](GraphVertex v) -> bool {
			if (degree(*graph_, v) == 1)
			{
				std::vector<GraphVertex> av = adjacent_vertices_through_edge(*graph_, v);
				const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
				const Vec3& q = value<Vec3>(*graph_, graph_vertex_position_, av[0]);
				Vec3 dir = p - q;

				BVH_Hit h = intersect_bvh({p, dir, 0, acc::inf});
				if (h.hit)
				{
					GraphVertex nv = add_vertex(*graph_);
					connect_vertices(*graph_, v, nv);
					value<Vec3>(*graph_, graph_vertex_position_, nv) = p + 0.6 * (h.pos - p);
				}
			}
			return true;
		});

		graph_provider_->emit_connectivity_changed(*graph_);
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());
	}

	void recenter_graph_from_surface()
	{
		foreach_cell(*graph_, [&](GraphVertex v) -> bool {
			uint32 d = degree(*graph_, v);
			if (d > 2)
			{
				Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
				Vec3 displ;
				Vec3 prev_displ;
				std::pair<uint32, Scalar> k_res;
				surface_kdt_->find_nn(p, &k_res);
				Vec3& nnp = value<Vec3>(*surface_, surface_vertex_position_, surface_vertices_[k_res.first]);
				displ = 0.01 * (p - nnp);
				p += displ;
				do
				{
					prev_displ = displ;
					surface_kdt_->find_nn(p, &k_res);
					Vec3& nnp = value<Vec3>(*surface_, surface_vertex_position_, surface_vertices_[k_res.first]);
					displ = 0.01 * (p - nnp);
					p += displ;
				} while (prev_displ.dot(displ) > 0);

				Vec3 cp = surface_bvh_->closest_point(p);
				value<Scalar>(*graph_, graph_vertex_radius_, v) = (cp - p).norm();
			}
			else if (d == 2)
			{
				Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
				std::vector<GraphVertex> av = adjacent_vertices_through_edge(*graph_, v);

				Vec3 dir = (value<Vec3>(*graph_, graph_vertex_position_, av[0]) - p).normalized() -
						   (p - value<Vec3>(*graph_, graph_vertex_position_, av[1])).normalized();
				dir.normalize();

				Vec3 displ;
				Vec3 prev_displ;
				std::pair<uint32, Scalar> k_res;
				surface_kdt_->find_nn(p, &k_res);
				Vec3& nnp = value<Vec3>(*surface_, surface_vertex_position_, surface_vertices_[k_res.first]);
				displ = dir * dir.dot(0.01 * (p - nnp));
				p += displ;
				do
				{
					prev_displ = displ;
					surface_kdt_->find_nn(p, &k_res);
					Vec3& nnp = value<Vec3>(*surface_, surface_vertex_position_, surface_vertices_[k_res.first]);
					displ = dir * dir.dot(0.01 * (p - nnp));
					p += displ;
				} while (prev_displ.dot(displ) > 0);

				Vec3 cp = surface_bvh_->closest_point(p);
				value<Scalar>(*graph_, graph_vertex_radius_, v) = (cp - p).norm();
			}
			else if (d == 1)
			{
				// TODO
			}

			return true;
		});

		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_radius_.get());
		// cgogn::io::export_CGR(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), "export.cgr");
	}

	void init_graph_radius_from_surface()
	{
		parallel_foreach_cell(*graph_, [&](GraphVertex v) -> bool {
			const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
			Vec3 cp = surface_bvh_->closest_point(p);
			value<Scalar>(*graph_, graph_vertex_radius_, v) = (cp - p).norm();
			return true;
		});
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_radius_.get());
		// cgogn::io::export_CGR(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), "export.cgr");
	}

	GRAPH* resample_graph(Scalar density)
	{
		if constexpr (std::is_same_v<GRAPH, cgogn::Graph>)
		{
			static uint32 count = 0;
			GRAPH* resampled_graph = graph_provider_->add_mesh("resampled_" + std::to_string(count++));
			auto resampled_graph_vertex_position = add_attribute<Vec3, GraphVertex>(*resampled_graph, "position");
			auto resampled_graph_vertex_radius = add_attribute<Scalar, GraphVertex>(*resampled_graph, "radius");

			modeling::resample_graph(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(),
									 *resampled_graph, resampled_graph_vertex_position.get(),
									 resampled_graph_vertex_radius.get(), density);

			graph_provider_->emit_connectivity_changed(*graph_);
			graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());
			graph_provider_->emit_attribute_changed(*graph_, graph_vertex_radius_.get());

			graph_provider_->emit_connectivity_changed(*resampled_graph);
			graph_provider_->emit_attribute_changed(*resampled_graph, resampled_graph_vertex_position.get());
			graph_provider_->emit_attribute_changed(*resampled_graph, resampled_graph_vertex_radius.get());

			return resampled_graph;
		}
		else
			return nullptr;
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

		// for (Scalar& r : *graph_vertex_radius_)
		// 	r = r / 1.1;

		contact_surface_ = surface_provider_->add_mesh("contact");
		volume_ = volume_provider_->add_mesh("hex");

		if constexpr (std::is_same_v<GRAPH, Graph>)
			hex_building_attributes_ = modeling::graph_to_hex(*graph_, *contact_surface_, *volume_);

		if constexpr (std::is_same_v<GRAPH, IncidenceGraph>)
			hex_building_attributes_ig_ = modeling::incidenceGraph_to_hex(*graph_, *contact_surface_, *volume_);

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

		surface_provider_->emit_connectivity_changed(*contact_surface_);
		volume_provider_->emit_connectivity_changed(*volume_);

		volume_vertex_position_ = get_attribute<Vec3, VolumeVertex>(*volume_, "position");
		volume_vertex_index_ = add_attribute<uint32, VolumeVertex>(*volume_, "vertex_index");
		volume_edge_index_ = add_attribute<uint32, VolumeEdge>(*volume_, "edge_index");
		volume_edge_target_length_ = add_attribute<Scalar, VolumeEdge>(*volume_, "target_length");

		volume_provider_->set_mesh_bb_vertex_position(*volume_, volume_vertex_position_);

		graph_provider_->emit_connectivity_changed(*graph_);
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());

		refresh_edge_target_length_ = true;
		refresh_volume_cells_indexing_ = true;
		refresh_volume_skin_ = true;
		refresh_solver_ = true;

		return volume_;
	}

	void add_volume_padding()
	{
		modeling::padding(*volume_);

		volume_provider_->emit_connectivity_changed(*volume_);
		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

		refresh_edge_target_length_ = true;
		refresh_volume_cells_indexing_ = true;
		refresh_volume_skin_ = true;
		refresh_solver_ = true;
	}

	void subdivide_slice()
	{
		if (selected_volume_faces_set_->size() == 1)
		{
			VolumeEdge e = modeling::find_fiber_dir(*volume_, *(selected_volume_faces_set_->begin()));
			CellCache<VOLUME> slice = modeling::get_slice(*volume_, e);
			modeling::cut_slice(*volume_, volume_vertex_position_.get(), slice);

			volume_provider_->emit_connectivity_changed(*volume_);
			volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

			refresh_edge_target_length_ = true;
			refresh_volume_cells_indexing_ = true;
			refresh_volume_skin_ = true;
			refresh_solver_ = true;
		}
	}

	void subdivide_all_slices()
	{
		if (selected_volume_faces_set_->size() == 1)
		{
			VolumeEdge e = modeling::find_fiber_dir(*volume_, *(selected_volume_faces_set_->begin()));
			CellCache<VOLUME> slices = modeling::surface_fiber_spread(*volume_, e);

			foreach_cell(slices, [&](VolumeEdge e) -> bool {
				CellCache<VOLUME> slice = modeling::get_slice(*volume_, e);
				modeling::cut_slice(*volume_, volume_vertex_position_.get(), slice);
				return true;
			});

			volume_provider_->emit_connectivity_changed(*volume_);
			volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

			refresh_edge_target_length_ = true;
			refresh_volume_cells_indexing_ = true;
			refresh_volume_skin_ = true;
			refresh_solver_ = true;
		}
	}

	void fiber_aligned_subdivision_from_input()
	{
		if (selected_volume_faces_set_->size() == 1)
		{
			CellMarker<VOLUME, VolumeEdge> edge_fibers(*volume_);
			VolumeEdge e = modeling::find_fiber_dir(*volume_, *(selected_volume_faces_set_->begin()));
			modeling::mark_mesh_fibers(*volume_, e, edge_fibers);

			modeling::fiber_aligned_subdivision(*volume_, edge_fibers);
			volume_provider_->emit_connectivity_changed(*volume_);
			volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

			refresh_edge_target_length_ = true;
			refresh_volume_cells_indexing_ = true;
			refresh_volume_skin_ = true;
			refresh_solver_ = true;
		}
	}

	// void subdivide_skin()
	// {
	// 	CellMarker<VOLUME, VolumeFace> visited_faces(*volume_);
	// 	CellCache<VOLUME> skin_cells(*volume_);

	// 	skin_cells.template build<VolumeVolume>([&](VolumeVolume w) {
	// 		bool adjacent_boundary = false;
	// 		foreach_incident_face(*volume_, w, [&](VolumeFace wf) -> bool {
	// 			adjacent_boundary = is_incident_to_boundary(*volume_, wf);
	// 			return !adjacent_boundary;
	// 		});
	// 		return adjacent_boundary;
	// 	});
	// 	foreach_cell(skin_cells, [&](VolumeVolume w) {
	// 		foreach_incident_face(*volume_, w, [&](VolumeFace wf) -> bool {
	// 			if (!visited_faces.is_marked(wf))
	// 			{
	// 				skin_cells.add(wf);
	// 				visited_faces.mark(wf);
	// 			}
	// 			return true;
	// 		});
	// 		return true;
	// 	});
	// 	modeling::primal_cut_all_volumes(
	// 		skin_cells,
	// 		[&](VolumeVertex v) {
	// 			std::vector<VolumeVertex> av = adjacent_vertices_through_edge(*volume_, v);
	// 			cgogn::value<Vec3>(*volume_, volume_vertex_position_, v) =
	// 				0.5 * (cgogn::value<Vec3>(*volume_, volume_vertex_position_, av[0]) +
	// 					   cgogn::value<Vec3>(*volume_, volume_vertex_position_, av[1]));
	// 		},
	// 		[&](VolumeVertex v) {
	// 			Vec3 center;
	// 			center.setZero();
	// 			uint32 count = 0;
	// 			foreach_adjacent_vertex_through_edge(*volume_, v, [&](VolumeVertex av) -> bool {
	// 				center += cgogn::value<Vec3>(*volume_, volume_vertex_position_, av);
	// 				++count;
	// 				return true;
	// 			});
	// 			center /= Scalar(count);
	// 			cgogn::value<Vec3>(*volume_, volume_vertex_position_, v) = center;
	// 		},
	// 		[&](VolumeVertex v) {
	// 			Vec3 center;
	// 			center.setZero();
	// 			uint32 count = 0;
	// 			foreach_adjacent_vertex_through_edge(*volume_, v, [&](VolumeVertex av) -> bool {
	// 				center += cgogn::value<Vec3>(*volume_, volume_vertex_position_, av);
	// 				++count;
	// 				return true;
	// 			});
	// 			center /= Scalar(count);
	// 			cgogn::value<Vec3>(*volume_, volume_vertex_position_, v) = center;
	// 		});

	// 	volume_provider_->emit_connectivity_changed(*volume_);
	// 	volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

	// 	return;
	// }

	void subdivide_volume()
	{
		CellCache<CMap3> cache(*volume_);
		cache.build<VolumeVolume>();
		cache.build<VolumeFace>();

		modeling::primal_cut_all_volumes(
			cache,
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

		volume_provider_->emit_connectivity_changed(*volume_);
		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

		refresh_edge_target_length_ = true;
		refresh_volume_cells_indexing_ = true;
		refresh_volume_skin_ = true;
		refresh_solver_ = true;
	}

	void refresh_volume_cells_indexing()
	{
		// std::cout << "refresh_volume_cells_indexing" << std::endl;

		uint32 vertex_idx = 0;
		foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			value<uint32>(*volume_, volume_vertex_index_, v) = vertex_idx++;
			return true;
		});

		uint32 edge_idx = 0;
		foreach_cell(*volume_, [&](VolumeEdge v) -> bool {
			value<uint32>(*volume_, volume_edge_index_, v) = edge_idx++;
			edge_idx++; // for the second orientation of the edge
			return true;
		});

		refresh_volume_cells_indexing_ = false;
	}

	void refresh_edge_target_length()
	{
		// std::cout << "refresh_edge_target_length" << std::endl;

		length_mean_ = 0.0;
		uint32 nbe = 0;

		foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			// auto vertices = incident_vertices(*volume_, e);

			std::vector<VolumeEdge> parallel_edges;
			parallel_edges.reserve(16);
			// std::vector<VolumeEdge> perpendicular_edges;
			// perpendicular_edges.reserve(16);

			Dart ed = e.dart;
			parallel_edges.push_back(VolumeEdge(ed)); // the edge itself
			// Dart c = phi<1, 2, 3>(*volume_, ed);
			// if (!is_boundary(*volume_, c))
			// 	parallel_edges.push_back(VolumeEdge(phi<2, 1>(*volume_, c)));
			do
			{
				Dart vd = phi<2, 1, 1>(*volume_, ed);
				parallel_edges.push_back(VolumeEdge(vd));
				// c = phi<1, 2, 3>(*volume_, vd);
				// if (!is_boundary(*volume_, c))
				// 	parallel_edges.push_back(VolumeEdge(phi<2, 1>(*volume_, c)));
				if (!is_boundary(*volume_, ed))
				{
					vd = phi<2, 1, 1>(*volume_, vd);
					parallel_edges.push_back(VolumeEdge(vd));
					// c = phi<1, 2, 3>(*volume_, vd);
					// if (!is_boundary(*volume_, c))
					// 	parallel_edges.push_back(VolumeEdge(phi<2, 1>(*volume_, c)));
				}
				else
					parallel_edges.push_back(VolumeEdge(ed)); // edge is on the boundary -> count twice
				// perpendicular_edges.push_back(VolumeEdge(phi1(*volume_, ed)));
				// perpendicular_edges.push_back(VolumeEdge(phi_1(*volume_, ed)));
				ed = phi<3, 2>(*volume_, ed);
			} while (ed != e.dart);

			// Dart ed2 = phi2(*volume_, e.dart);
			// c = phi<1, 2, 3>(*volume_, ed2);
			// if (!is_boundary(*volume_, c))
			// 	parallel_edges.push_back(VolumeEdge(phi<2, 1>(*volume_, c)));
			// do
			// {
			// 	Dart vd = phi<2, 1, 1>(*volume_, ed2);
			// 	c = phi<1, 2, 3>(*volume_, vd);
			// 	if (!is_boundary(*volume_, c))
			// 		parallel_edges.push_back(VolumeEdge(phi<2, 1>(*volume_, c)));
			// 	if (!is_boundary(*volume_, ed))
			// 	{
			// 		vd = phi<2, 1, 1>(*volume_, vd);
			// 		c = phi<1, 2, 3>(*volume_, vd);
			// 		if (!is_boundary(*volume_, c))
			// 			parallel_edges.push_back(VolumeEdge(phi<2, 1>(*volume_, c)));
			// 	}
			// 	ed2 = phi<3, 2>(*volume_, ed2);
			// } while (ed2 != phi2(*volume_, e.dart));

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
			// value<Scalar>(*volume_, volume_edge_target_length_, e) =
			// 	(2.0 * parallel_edges_mean_length + perpendicular_edges_mean_length) / 3.0;

			length_mean_ += parallel_edges_mean_length;
			++nbe;

			return true;
		});

		length_mean_ /= nbe;

		refresh_edge_target_length_ = false;
		refresh_solver_matrix_values_only_ = true;
	}

	void refresh_volume_skin()
	{
		// std::cout << "refresh_volume_skin" << std::endl;

		if (!volume_skin_)
			volume_skin_ = surface_provider_->add_mesh("volume_skin");

		surface_provider_->clear_mesh(*volume_skin_);

		volume_skin_vertex_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "position");
		volume_skin_vertex_normal_ = get_or_add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "normal");
		volume_skin_vertex_index_ = get_or_add_attribute<uint32, SurfaceVertex>(*volume_skin_, "vertex_index");
		volume_skin_vertex_volume_vertex_ =
			get_or_add_attribute<VolumeVertex, SurfaceVertex>(*volume_skin_, "hex_vertex");
		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), *volume_skin_,
										 volume_skin_vertex_position_.get(), volume_skin_vertex_volume_vertex_.get());

		uint32 nb_vertices = 0;
		foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			value<uint32>(*volume_skin_, volume_skin_vertex_index_, v) = nb_vertices++;
			return true;
		});

		surface_provider_->emit_connectivity_changed(*volume_skin_);

		refresh_volume_skin_ = false;
	}

	void refresh_solver_matrix_values(Scalar fit_to_data)
	{
		// std::cout << "refresh_solver_matrix_values" << std::endl;

		if (refresh_volume_skin_)
			refresh_volume_skin();
		if (refresh_volume_cells_indexing_)
			refresh_volume_cells_indexing();
		if (refresh_edge_target_length_)
			refresh_edge_target_length();

		MeshData<VOLUME>& mdv = volume_provider_->mesh_data(*volume_);
		MeshData<SURFACE>& mds = surface_provider_->mesh_data(*volume_skin_);
		uint32 nb_volume_vertices = mdv.template nb_cells<VolumeVertex>();
		uint32 nb_boundary_vertices = mds.template nb_cells<SurfaceVertex>();
		uint32 nb_volume_edges = mdv.template nb_cells<VolumeEdge>();
		uint32 nb_oriented_edges = nb_volume_edges * 2;

		std::vector<Eigen::Triplet<Scalar>> triplets;
		triplets.reserve(nb_oriented_edges * 2 + nb_boundary_vertices);

		foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			uint32 eidx = value<uint32>(*volume_, volume_edge_index_, e);

			Scalar target_length = value<Scalar>(*volume_, volume_edge_target_length_, e) / length_mean_;

			auto vertices = incident_vertices(*volume_, e);
			uint32 vidx1 = value<uint32>(*volume_, volume_vertex_index_, vertices[0]);
			uint32 vidx2 = value<uint32>(*volume_, volume_vertex_index_, vertices[1]);

			triplets.push_back(Eigen::Triplet<Scalar>(int(eidx), int(vidx1), -1 / target_length));
			triplets.push_back(Eigen::Triplet<Scalar>(int(eidx), int(vidx2), 1 / target_length));

			triplets.push_back(Eigen::Triplet<Scalar>(int(eidx + 1), int(vidx1), 1 / target_length));
			triplets.push_back(Eigen::Triplet<Scalar>(int(eidx + 1), int(vidx2), -1 / target_length));

			return true;
		});

		// set constrained vertices
		foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			uint32 boundary_vertex_idx = value<uint32>(*volume_skin_, volume_skin_vertex_index_, v);
			uint32 volume_vertex_idx =
				value<uint32>(*volume_, volume_vertex_index_,
							  value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v));
			triplets.push_back(Eigen::Triplet<Scalar>(int(nb_oriented_edges + boundary_vertex_idx),
													  int(volume_vertex_idx), fit_to_data));
			return true;
		});

		solver_matrix_.setZero();
		solver_matrix_.resize(nb_oriented_edges + nb_boundary_vertices, nb_volume_vertices);
		solver_matrix_.setFromTriplets(triplets.begin(), triplets.end());

		refresh_solver_matrix_values_only_ = false;
	}

	bool is_inside(const Vec3& p)
	{
		std::pair<uint32, Vec3> cp;
		surface_bvh_->closest_point(p, &cp);
		Vec3 dir = (cp.second - p).normalized();
		Vec3 n = geometry::normal(*surface_, surface_faces_[cp.first], surface_vertex_position_.get());
		return dir.dot(n) >= 0.0;
	}

	void project_on_surface()
	{
		if (refresh_volume_skin_)
			refresh_volume_skin();

		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			const Vec3& p = value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v);
			Vec3 n{0, 0, 0};
			foreach_incident_face(*volume_skin_, v, [&](SurfaceFace f) -> bool {
				Vec3 nf = geometry::normal(*volume_skin_, f, volume_skin_vertex_position_.get());
				Vec3 cf = geometry::centroid<Vec3>(*volume_skin_, f, volume_skin_vertex_position_.get());
				bool inside = is_inside(cf);
				if (!inside)
					nf *= -1;
				BVH_Hit h = intersect_bvh({cf, nf, 0, acc::inf});
				if (h.hit)
					n += inside ? h.pos - cf : cf - h.pos;
				return true;
			});
			n.normalize();

			if (!is_inside(p))
				n *= -1;

			BVH_Hit h = intersect_bvh({p, n, 0, acc::inf});
			Vec3 pos;
			if (h.hit)
				pos = h.pos;
			else
				pos = surface_bvh_->closest_point(p);

			value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v) = pos;
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v)) = pos;
			return true;
		});

		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
		surface_provider_->emit_attribute_changed(*volume_skin_, volume_skin_vertex_position_.get());

		refresh_edge_target_length_ = true;
	}

	void regularize_surface_vertices(Scalar fit_to_data)
	{
		if (refresh_volume_skin_)
			refresh_volume_skin();

		MeshData<SURFACE>& mds = surface_provider_->mesh_data(*volume_skin_);
		uint32 nb_vertices = mds.template nb_cells<SurfaceVertex>();

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * nb_vertices, nb_vertices);
		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_vertices * 10);
		Eigen::MatrixXd b(2 * nb_vertices, 3);

		foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_skin_, volume_skin_vertex_index_, v);
			uint32 nbv = 0;
			foreach_adjacent_vertex_through_edge(*volume_skin_, v, [&](SurfaceVertex av) -> bool {
				uint32 avidx = value<uint32>(*volume_skin_, volume_skin_vertex_index_, av);
				Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
				++nbv;
				return true;
			});
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(nbv)));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), fit_to_data));
			b(vidx, 0) = 0;
			b(vidx, 1) = 0;
			b(vidx, 2) = 0;
			Vec3 pos = surface_bvh_->closest_point(value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v));
			b(nb_vertices + vidx, 0) = fit_to_data * pos[0];
			b(nb_vertices + vidx, 1) = fit_to_data * pos[1];
			b(nb_vertices + vidx, 2) = fit_to_data * pos[2];
			return true;
		});
		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At = A.transpose();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(At * A);

		Eigen::MatrixXd vpos(nb_vertices, 3);
		vpos = solver.solve(At * b);

		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_skin_, volume_skin_vertex_index_, v);
			Vec3& pos = value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v);
			pos[0] = vpos(vidx, 0);
			pos[1] = vpos(vidx, 1);
			pos[2] = vpos(vidx, 2);
			return true;
		});

		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			const Vec3& pos = value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v);
			// Vec3 cp = surface_bvh_->closest_point(pos);
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v)) = pos;
			return true;
		});

		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
		surface_provider_->emit_attribute_changed(*volume_skin_, volume_skin_vertex_position_.get());

		refresh_edge_target_length_ = true;
	}

	void relocate_interior_vertices()
	{
		if (refresh_volume_cells_indexing_)
			refresh_volume_cells_indexing();

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A =
			geometry::topo_laplacian_matrix(*volume_, volume_vertex_index_.get());
		A = A.diagonal().asDiagonal().inverse() * A;

		Eigen::MatrixXd vpos(A.cols(), 3);
		Eigen::MatrixXd lapl(A.rows(), 3);

		for (uint32 i = 0; i < 10; ++i)
		{
			parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
				const Vec3& pv = value<Vec3>(*volume_, volume_vertex_position_, v);
				uint32 vidx = value<uint32>(*volume_, volume_vertex_index_, v);
				vpos(vidx, 0) = pv[0];
				vpos(vidx, 1) = pv[1];
				vpos(vidx, 2) = pv[2];
				return true;
			});

			lapl = A * vpos;

			parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
				if (!is_incident_to_boundary(*volume_, v))
				{
					uint32 vidx = value<uint32>(*volume_, volume_vertex_index_, v);
					Vec3 l;
					l[0] = lapl(vidx, 0);
					l[1] = lapl(vidx, 1);
					l[2] = lapl(vidx, 2);
					value<Vec3>(*volume_, volume_vertex_position_, v) -= 0.1 * l;
				}
				return true;
			});
		}

		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());

		refresh_edge_target_length_ = true;
	}

	void optimize_volume_vertices(Scalar fit_to_data, bool use_skin_as_boundary = false)
	{
		MeshData<SURFACE>& mds = surface_provider_->mesh_data(*surface_);

		if (refresh_edge_target_length_)
			refresh_edge_target_length();
		if (refresh_volume_skin_)
			refresh_volume_skin();
		if (refresh_volume_cells_indexing_)
			refresh_volume_cells_indexing();
		if (refresh_solver_)
		{
			refresh_solver_matrix_values(fit_to_data);
			if (solver_)
				delete solver_;
			solver_ = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>>();
			Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A = solver_matrix_.transpose() * solver_matrix_;
			solver_->analyzePattern(A);
			solver_->factorize(A);
			refresh_solver_ = false;
			refresh_solver_matrix_values_only_ = false;
		}
		else if (refresh_solver_matrix_values_only_)
		{
			refresh_solver_matrix_values(fit_to_data);
			Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A = solver_matrix_.transpose() * solver_matrix_;
			solver_->factorize(A);
			refresh_solver_matrix_values_only_ = false;
		}

		MeshData<VOLUME>& mdv = volume_provider_->mesh_data(*volume_);
		uint32 nb_volume_edges = mdv.template nb_cells<VolumeEdge>();
		uint32 nb_oriented_edges = nb_volume_edges * 2;

		Eigen::MatrixXd x(solver_matrix_.cols(), 3);
		Eigen::MatrixXd b(solver_matrix_.rows(), 3);

		parallel_foreach_cell(*volume_, [&](VolumeEdge e) -> bool {
			uint32 eidx = value<uint32>(*volume_, volume_edge_index_, e);

			auto vertices = incident_vertices(*volume_, e);
			const Vec3& pos1 = value<Vec3>(*volume_, volume_vertex_position_, vertices[0]);
			const Vec3& pos2 = value<Vec3>(*volume_, volume_vertex_position_, vertices[1]);

			Vec3 edge1 = (pos2 - pos1).normalized();
			Vec3 edge2 = (pos1 - pos2).normalized();

			Vec3 target_n1 = edge1; // (0, 0, 0);
			Vec3 target_n2 = edge2; // (0, 0, 0);

			Dart d = e.dart;
			do
			{
				if (!is_boundary(*volume_, d))
				{
					const Vec3& p2 =
						value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi<2, 1, 1>(*volume_, d)));
					const Vec3& p3 = value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi_1(*volume_, d)));
					Vec3 n = geometry::normal(pos1, p2, p3).normalized();
					if (edge1.dot(n) > 0)
						target_n1 += n;
				}
				d = phi<3, 2>(*volume_, d);
			} while (d != e.dart);
			target_n1.normalize();

			d = phi2(*volume_, e.dart);
			do
			{
				if (!is_boundary(*volume_, d))
				{
					const Vec3& p2 =
						value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi<2, 1, 1>(*volume_, d)));
					const Vec3& p3 = value<Vec3>(*volume_, volume_vertex_position_, VolumeVertex(phi_1(*volume_, d)));
					Vec3 n = geometry::normal(pos2, p2, p3).normalized();
					if (edge2.dot(n) > 0)
						target_n2 += n;
				}
				d = phi<3, 2>(*volume_, d);
			} while (d != phi2(*volume_, e.dart));
			target_n2.normalize();

			target_n1 *= length_mean_;
			target_n2 *= length_mean_;

			b.coeffRef(eidx, 0) = target_n1[0];
			b.coeffRef(eidx, 1) = target_n1[1];
			b.coeffRef(eidx, 2) = target_n1[2];

			b.coeffRef(eidx + 1, 0) = target_n2[0];
			b.coeffRef(eidx + 1, 1) = target_n2[1];
			b.coeffRef(eidx + 1, 2) = target_n2[2];

			return true;
		});

		if (use_skin_as_boundary)
		{
			parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
				const Vec3& pos = value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v);
				uint32 boundary_vertex_idx = value<uint32>(*volume_skin_, volume_skin_vertex_index_, v);
				b.coeffRef(nb_oriented_edges + boundary_vertex_idx, 0) = fit_to_data * pos[0];
				b.coeffRef(nb_oriented_edges + boundary_vertex_idx, 1) = fit_to_data * pos[1];
				b.coeffRef(nb_oriented_edges + boundary_vertex_idx, 2) = fit_to_data * pos[2];
				return true;
			});
		}
		else
		{
			parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
				const Vec3& p = value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v);
				Vec3 n{0, 0, 0};
				foreach_incident_face(*volume_skin_, v, [&](SurfaceFace f) -> bool {
					Vec3 nf = geometry::normal(*volume_skin_, f, volume_skin_vertex_position_.get());
					Vec3 cf = geometry::centroid<Vec3>(*volume_skin_, f, volume_skin_vertex_position_.get());
					bool inside = is_inside(cf);
					if (!inside)
						nf *= -1;
					BVH_Hit h = intersect_bvh({cf, nf, 0, acc::inf});
					if (h.hit)
						n += inside ? h.pos - cf : cf - h.pos;
					return true;
				});
				n.normalize();

				if (!is_inside(p))
					n *= -1;

				BVH_Hit h = intersect_bvh({p, n, 0, acc::inf});
				Vec3 pos;
				if (h.hit)
					pos = h.pos;
				else
					pos = surface_bvh_->closest_point(p);

				uint32 boundary_vertex_idx = value<uint32>(*volume_skin_, volume_skin_vertex_index_, v);
				b.coeffRef(nb_oriented_edges + boundary_vertex_idx, 0) = fit_to_data * pos[0];
				b.coeffRef(nb_oriented_edges + boundary_vertex_idx, 1) = fit_to_data * pos[1];
				b.coeffRef(nb_oriented_edges + boundary_vertex_idx, 2) = fit_to_data * pos[2];
				return true;
			});
		}

		x = solver_->solve(solver_matrix_.transpose() * b);

		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_, volume_vertex_index_, v);
			Vec3& pos = value<Vec3>(*volume_, volume_vertex_position_, v);
			pos[0] = x(vidx, 0);
			pos[1] = x(vidx, 1);
			pos[2] = x(vidx, 2);
			return true;
		});

		// update volume_skin vertex position
		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v) =
				value<Vec3>(*volume_, volume_vertex_position_,
							value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v));
			return true;
		});

		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
		surface_provider_->emit_attribute_changed(*volume_skin_, volume_skin_vertex_position_.get());
	}

	void rigid_register_volume_mesh()
	{
		if (refresh_volume_skin_)
			refresh_volume_skin();

		geometry::rigid_register_mesh(*volume_, volume_vertex_position_.get(), *surface_,
									  surface_vertex_position_.get());

		// update volume_skin vertex position
		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v) =
				value<Vec3>(*volume_, volume_vertex_position_,
							value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v));
			return true;
		});

		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
		surface_provider_->emit_attribute_changed(*volume_skin_, volume_skin_vertex_position_.get());
	}

	void non_rigid_register_volume_mesh(Scalar registration_fit, Scalar optimization_fit, bool init_surface_steady_pos,
										geometry::ProximityPolicy proximity)
	{
		if (refresh_volume_skin_)
			refresh_volume_skin();

		for (uint32 i = 0; i < 5; ++i)
		{
			geometry::non_rigid_register_mesh(*volume_skin_, volume_skin_vertex_position_, *surface_,
											  surface_vertex_position_.get(), registration_fit, false,
											  init_surface_steady_pos, proximity);
			optimize_volume_vertices(optimization_fit, true);
		}
	}

	void snapshot_volume_vertex_position()
	{
		static uint32 count = 1;
		auto pos = add_attribute<Vec3, VolumeVertex>(*volume_, "position_" + std::to_string(count));
		pos->copy(volume_vertex_position_.get());
		animate_volume_vertex_positions_.push_back(pos);
		++count;
	}

	void save_volume_animation()
	{
		std::ofstream out_file;
		out_file.open(volume_provider_->mesh_name(*volume_) + "_animation");
		uint32 nb_vertices = volume_provider_->mesh_data(*volume_).template nb_cells<VolumeVertex>();
		out_file << nb_vertices << " " << animate_volume_vertex_positions_.size() << "\n";

		for (uint32 i = volume_->attribute_containers_[VolumeVertex::ORBIT].first_index(),
					end = volume_->attribute_containers_[VolumeVertex::ORBIT].last_index();
			 i != end; i = volume_->attribute_containers_[VolumeVertex::ORBIT].next_index(i))
		{
			for (uint32 j = 0; j < animate_volume_vertex_positions_.size(); ++j)
			{
				const Vec3& p = (*animate_volume_vertex_positions_[j])[i];
				out_file << p[0] << " " << p[1] << " " << p[2] << " ";
			}
			out_file << "\n";
		}
		out_file.close();
	}

	void load_volume_animation()
	{
		for (uint32 i = 0; i < animate_volume_vertex_positions_.size(); ++i)
			remove_attribute<VolumeVertex>(*volume_, animate_volume_vertex_positions_[i]);
		animate_volume_vertex_positions_.clear();

		std::ifstream fp(volume_provider_->mesh_name(*volume_) + "_animation", std::ios::in);
		std::string line;
		line.reserve(512u);

		const uint32 nb_vertices = io::read_uint(fp, line);
		const uint32 nb_animation_position = io::read_uint(fp, line);

		for (uint32 i = 0u; i < nb_animation_position; ++i)
		{
			auto pos = add_attribute<Vec3, VolumeVertex>(*volume_, "position_" + std::to_string(i));
			animate_volume_vertex_positions_.push_back(pos);
		}

		uint32 vertex_id = volume_->attribute_containers_[VolumeVertex::ORBIT].first_index();
		for (uint32 i = 0u; i < nb_vertices; ++i)
		{
			for (uint32 j = 0u; j < nb_animation_position; ++j)
			{
				float64 x = io::read_double(fp, line);
				float64 y = io::read_double(fp, line);
				float64 z = io::read_double(fp, line);
				(*animate_volume_vertex_positions_[j])[vertex_id] = {x, y, z};
			}
			vertex_id = volume_->attribute_containers_[VolumeVertex::ORBIT].next_index(vertex_id);
		}

		volume_vertex_position_->copy(animate_volume_vertex_positions_.back().get());
		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
	}

	void compute_volumes_quality()
	{
		auto corner_frame = add_attribute<Mat3, VolumeVertex2>(*volume_, "__corner_frame");
		auto hex_frame = add_attribute<Mat3, VolumeVolume>(*volume_, "__hex_frame");

		auto scaled_jacobian = get_or_add_attribute<Scalar, VolumeVolume>(*volume_, "scaled_jacobian");
		auto jacobian = get_or_add_attribute<Scalar, VolumeVolume>(*volume_, "jacobian");
		auto max_froebnius = get_or_add_attribute<Scalar, VolumeVolume>(*volume_, "max_froebnius");
		auto mean_froebnius = get_or_add_attribute<Scalar, VolumeVolume>(*volume_, "mean_froebnius");

		geometry::compute_hex_frame(*volume_, volume_vertex_position_.get(), corner_frame.get(), hex_frame.get());
		geometry::compute_scaled_jacobian(*volume_, corner_frame.get(), hex_frame.get(), scaled_jacobian.get());
		geometry::compute_jacobian(*volume_, corner_frame.get(), hex_frame.get(), jacobian.get());
		geometry::compute_maximum_aspect_frobenius(*volume_, corner_frame.get(), max_froebnius.get());
		geometry::compute_mean_aspect_frobenius(*volume_, corner_frame.get(), mean_froebnius.get());

		Scalar min_scaled_jacobian = std::numeric_limits<Scalar>::max();
		Scalar mean_scaled_jacobian = 0.0;
		for (Scalar s : *scaled_jacobian)
		{
			mean_scaled_jacobian += s;
			if (s < min_scaled_jacobian)
				min_scaled_jacobian = s;
		}
		mean_scaled_jacobian /= scaled_jacobian->size();
		std::cout << "mean scaled jacobian = " << mean_scaled_jacobian
				  << " / min scaled jacobian = " << min_scaled_jacobian << std::endl;

		volume_provider_->emit_attribute_changed(*volume_, scaled_jacobian.get());
		volume_provider_->emit_attribute_changed(*volume_, jacobian.get());
		volume_provider_->emit_attribute_changed(*volume_, max_froebnius.get());
		volume_provider_->emit_attribute_changed(*volume_, mean_froebnius.get());

		remove_attribute<VolumeVolume>(*volume_, hex_frame.get());
		remove_attribute<VolumeVertex2>(*volume_, corner_frame.get());
	}

	void export_subdivided_skin()
	{
		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");
		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get());
		// modeling::catmull_clark_approx(volume_skin, volume_skin_vertex_position_.get(), 2);
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
		volume_vertex_index_ = get_or_add_attribute<uint32, VolumeVertex>(*volume_, "vertex_index");
		volume_edge_index_ = get_or_add_attribute<uint32, VolumeEdge>(*volume_, "edge_index");
		volume_edge_target_length_ = get_or_add_attribute<Scalar, VolumeEdge>(*volume_, "target_length");
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

			uint32 nb_vertices = surface_provider_->mesh_data(*surface_).template nb_cells<SurfaceVertex>();
			uint32 nb_faces = surface_provider_->mesh_data(*surface_).template nb_cells<SurfaceFace>();

			auto vertex_index = get_or_add_attribute<uint32, SurfaceVertex>(*surface_, "__bvh_vertex_index");

			std::vector<Vec3> vertex_position;
			vertex_position.reserve(nb_vertices);
			surface_vertices_.clear();
			surface_vertices_.reserve(nb_vertices);
			uint32 idx = 0;
			foreach_cell(*surface_, [&](SurfaceVertex v) -> bool {
				value<uint32>(*surface_, vertex_index, v) = idx++;
				surface_vertices_.push_back(v);
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

			if (surface_kdt_)
				delete surface_kdt_;

			surface_kdt_ = new acc::KDTree<3, uint32>(vertex_position);
		}
	}

protected:
	void start_animate_volume()
	{
		if (animate_volume_vertex_positions_.size() > 1)
		{
			animate_volume_ = true;
			animate_volume_start_time_ = App::frame_time_;
			app_.start_timer(40, [this]() -> bool { return !animate_volume_; });
		}
	}

	void stop_animate_volume()
	{
		animate_volume_ = false;
		volume_vertex_position_->copy(animate_volume_vertex_positions_.back().get());
		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
	}

	void animate_volume()
	{
		double cycle_duration = animate_volume_slot_duration_ * (animate_volume_vertex_positions_.size() - 1);
		double time_pos = std::fmod(App::frame_time_ - animate_volume_start_time_, cycle_duration);

		uint32 slot_index = static_cast<uint32>(std::floor(time_pos / animate_volume_slot_duration_));
		double slot_pos = std::fmod(time_pos, animate_volume_slot_duration_);

		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			value<Vec3>(*volume_, volume_vertex_position_, v) =
				(1.0 - slot_pos) * value<Vec3>(*volume_, animate_volume_vertex_positions_[slot_index], v) +
				slot_pos * value<Vec3>(*volume_, animate_volume_vertex_positions_[slot_index + 1], v);
			return true;
		});

		volume_provider_->emit_attribute_changed(*volume_, volume_vertex_position_.get());
	}

	void left_panel() override
	{
		ImGui::TextUnformatted("Graph");
		imgui_mesh_selector(graph_provider_, graph_, "Graph", [&](GRAPH& g) { set_current_graph(&g); });
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
		imgui_mesh_selector(surface_provider_, surface_, "Surface", [&](SURFACE& s) { set_current_surface(&s); });
		if (surface_)
		{
			imgui_combo_attribute<SurfaceVertex, Vec3>(*surface_, surface_vertex_position_, "Position##surface",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   set_current_surface_vertex_position(attribute);
													   });
		}

		ImGui::Separator();
		ImGui::TextUnformatted("Graph Operations");
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_)
		{
			if (ImGui::Button("Init radius from edge length"))
				init_graph_radius_from_edge_length();
		}
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_ && surface_ && surface_vertex_position_)
		{
			if (ImGui::Button("Init radius from surface"))
				init_graph_radius_from_surface();
			if (ImGui::Button("Recenter graph from surface"))
				recenter_graph_from_surface();
			if (ImGui::Button("Extend graph extremities"))
				extend_graph_extremities();
		}
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_)
		{
			static float graph_resample_density = 0.5f;
			ImGui::SliderFloat("Graph resampling density", &graph_resample_density, 0.0, 2.0);
			if (ImGui::Button("Resample graph"))
				resample_graph(graph_resample_density);
		}

		ImGui::Separator();
		if (graph_ && graph_vertex_position_)
		{
			if (!volume_)
				if (ImGui::Button("Build hex mesh"))
					build_hex_mesh();
		}
		if (volume_)
		{
			MeshData<VOLUME>& md = volume_provider_->mesh_data(*volume_);
			// float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			ImGui::TextUnformatted("HexMesh Connectivity Ops");

			// if (ImGui::Button("Export subdivided skin"))
			// 	export_subdivided_skin();
			if (ImGui::Button("Add volume padding"))
				add_volume_padding();
			// if (ImGui::Button("Subdivide length wise"))
			// 	subdivide_volume_length_wise();
			// if (ImGui::Button("Subdivide width wise"))
			// 	subdivide_volume_width_wise();
			// if (ImGui::Button("Find Fibers"))
			// 	fiber_aligned_subdivision_from_input();
			if (ImGui::Button("Subdivide volume"))
				subdivide_volume();
			// if (ImGui::Button("Subdivide skin"))
			// 	subdivide_skin();
			imgui_combo_cells_set(md, selected_volume_faces_set_, "Faces",
								  [&](CellsSet<VOLUME, VolumeFace>* cs) { selected_volume_faces_set_ = cs; });
			if (selected_volume_faces_set_)
			{
				if (ImGui::Button("Subdivide slice"))
					subdivide_slice();
				if (ImGui::Button("Subdivide all slices"))
					subdivide_all_slices();
			}

			ImGui::Separator();
			ImGui::TextUnformatted("HexMesh Geometry Ops");

			if (ImGui::Button("Relocate interior vertices"))
				relocate_interior_vertices();
			ImGui::Checkbox("Refresh edge target length", &refresh_edge_target_length_);
			refresh_solver_matrix_values_only_ = true;

			static float optimize_fit_to_surface = 1.0f;
			ImGui::SliderFloat("Optimize volume - Fit to surface", &optimize_fit_to_surface, 0.1, 10.0);
			if (ImGui::Button("Optimize volume vertices (inside skin)"))
				optimize_volume_vertices(optimize_fit_to_surface, true);

			if (surface_ && surface_vertex_position_)
			{
				if (ImGui::Button("Optimize volume vertices"))
					optimize_volume_vertices(optimize_fit_to_surface);

				ImGui::Separator();

				if (ImGui::Button("Project on surface"))
					project_on_surface();
				static float regularize_fit_to_data = 5.0f;
				ImGui::SliderFloat("Regularize surface - Fit to data", &regularize_fit_to_data, 0.0, 20.0);
				if (ImGui::Button("Regularize surface vertices"))
					regularize_surface_vertices(regularize_fit_to_data);

				ImGui::Separator();

				if (ImGui::Button("Rigid register volume mesh"))
					rigid_register_volume_mesh();
				static bool init_steady_pos = false;
				ImGui::Checkbox("Init steady pos", &init_steady_pos);
				static float registration_fit_to_target = 0.05f;
				ImGui::SliderFloat("Registration - Fit to target", &registration_fit_to_target, 0.01, 0.5);
				static int proximity = geometry::NEAREST_POINT;
				ImGui::RadioButton("Nearest", &proximity, geometry::NEAREST_POINT);
				ImGui::SameLine();
				ImGui::RadioButton("Normal ray", &proximity, geometry::NORMAL_RAY);
				if (ImGui::Button("Non-rigid register volume mesh"))
					non_rigid_register_volume_mesh(registration_fit_to_target, optimize_fit_to_surface, init_steady_pos,
												   geometry::ProximityPolicy(proximity));
				if (ImGui::Button("Snapshot volume position"))
					snapshot_volume_vertex_position();
				if (!animate_volume_)
				{
					if (ImGui::Button("Start animation"))
						start_animate_volume();
				}
				else
				{
					if (ImGui::Button("Stop animation"))
						stop_animate_volume();
				}
				if (ImGui::Button("Save volume animation"))
					save_volume_animation();
				if (ImGui::Button("Load volume animation"))
					load_volume_animation();
			}

			ImGui::Separator();

			if (ImGui::Button("Compute volumes quality"))
				compute_volumes_quality();
		}
	}

private:
	MeshProvider<GRAPH>* graph_provider_ = nullptr;
	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<VOLUME>* volume_provider_ = nullptr;

	GRAPH* graph_ = nullptr;
	std::shared_ptr<GraphAttribute<Vec3>> graph_vertex_position_ = nullptr;
	std::shared_ptr<GraphAttribute<Scalar>> graph_vertex_radius_ = nullptr;

	SURFACE* surface_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_position_ = nullptr;
	acc::BVHTree<uint32, Vec3>* surface_bvh_ = nullptr;
	std::vector<SurfaceFace> surface_faces_;
	acc::KDTree<3, uint32>* surface_kdt_ = nullptr;
	std::vector<SurfaceVertex> surface_vertices_;

	SURFACE* contact_surface_ = nullptr;

	VOLUME* volume_ = nullptr;
	std::shared_ptr<VolumeAttribute<Vec3>> volume_vertex_position_ = nullptr;
	std::shared_ptr<VolumeAttribute<uint32>> volume_vertex_index_ = nullptr;
	std::shared_ptr<VolumeAttribute<uint32>> volume_edge_index_ = nullptr;
	std::shared_ptr<VolumeAttribute<Scalar>> volume_edge_target_length_ = nullptr;
	Scalar length_mean_;
	bool refresh_edge_target_length_ = true;
	bool refresh_volume_cells_indexing_ = true;

	SURFACE* volume_skin_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_position_ = nullptr;
	std::shared_ptr<SurfaceAttribute<uint32>> volume_skin_vertex_index_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_normal_ = nullptr;
	std::shared_ptr<SurfaceAttribute<VolumeVertex>> volume_skin_vertex_volume_vertex_ = nullptr;
	bool refresh_volume_skin_ = true;

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> solver_matrix_;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>>* solver_ = nullptr;
	bool refresh_solver_matrix_values_only_ = true;
	bool refresh_solver_ = true;

	bool animate_volume_ = false;
	double animate_volume_start_time_;
	double animate_volume_slot_duration_ = 1.0;
	std::vector<std::shared_ptr<VolumeAttribute<Vec3>>> animate_volume_vertex_positions_;

	CellMarker<VOLUME, VolumeFace>* transversal_faces_marker_ = nullptr;
	CellsSet<VOLUME, VolumeFace>* selected_volume_faces_set_ = nullptr;

	std::tuple<modeling::GAttributes, modeling::M2Attributes, modeling::M3Attributes> hex_building_attributes_;
	std::tuple<modeling::IG_GAttributes, modeling::IG_M2Attributes, modeling::IG_M3Attributes>
		hex_building_attributes_ig_;

	std::shared_ptr<boost::synapse::connection> timer_connection_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_TUBULAR_MESH_H_
