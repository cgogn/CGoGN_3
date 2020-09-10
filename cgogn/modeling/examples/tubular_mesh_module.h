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
#include <cgogn/geometry/types/grid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/algos/distance.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/algos/filtering.h>
#include <cgogn/geometry/algos/hex_quality.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/picking.h>

#include <cgogn/modeling/algos/graph_resampling.h>
#include <cgogn/modeling/algos/graph_to_hex.h>

#include <Eigen/Sparse>

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

	using Grid = geometry::Grid<10, 10, 10, SURFACE>;

public:
	TubularMesh(const App& app)
		: ViewModule(app, "TubularMesh"), graph_(nullptr), graph_vertex_position_(nullptr),
		  graph_vertex_radius_(nullptr), surface_(nullptr), surface_vertex_position_(nullptr), surface_grid_(nullptr),
		  volume_(nullptr)
	{
	}

	~TubularMesh()
	{
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

	void init_graph_radius_from_edge_length()
	{
		Scalar l = geometry::mean_edge_length(*graph_, graph_vertex_position_.get());
		graph_vertex_radius_->fill(l / 4.0);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
	}

	void extend_graph_extremities()
	{
		using SelectedFace = std::tuple<SurfaceFace, Vec3, Scalar>;
		foreach_cell(*graph_, [&](GraphVertex v) -> bool {
			if (degree(*graph_, v) == 1)
			{
				std::vector<GraphVertex> av = adjacent_vertices_through_edge(*graph_, v);
				const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
				const Vec3& q = value<Vec3>(*graph_, graph_vertex_position_, av[0]);
				Vec3 dir = p - q;
				std::vector<SelectedFace> selectedfaces =
					geometry::internal::picking(*surface_, surface_vertex_position_.get(), p, p + dir);
				if (selectedfaces.size() > 0)
				{
					Vec3 pos = std::get<1>(selectedfaces[0]);
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

	void init_graph_radius_from_surface()
	{
		parallel_foreach_cell(*graph_, [&](Graph::Vertex v) -> bool {
			const Vec3& p = value<Vec3>(*graph_, graph_vertex_position_, v);
			Vec3 cp = geometry::closest_point_on_surface(*surface_, surface_vertex_position_.get(), p);
			value<Scalar>(*graph_, graph_vertex_radius_, v) = (cp - p).norm();
			return true;
		});
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
	}

	void resample_graph()
	{
		GRAPH* resampled_graph = graph_provider_->add_mesh("resampled");
		auto resampled_graph_vertex_position = add_attribute<Vec3, GraphVertex>(*resampled_graph, "position");
		auto resampled_graph_vertex_radius = add_attribute<Scalar, GraphVertex>(*resampled_graph, "radius");

		modeling::resample_graph(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), *resampled_graph,
								 resampled_graph_vertex_position.get(), resampled_graph_vertex_radius.get());

		graph_provider_->emit_connectivity_changed(graph_);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());

		graph_provider_->emit_connectivity_changed(resampled_graph);
		graph_provider_->emit_attribute_changed(resampled_graph, resampled_graph_vertex_position.get());
		graph_provider_->emit_attribute_changed(resampled_graph, resampled_graph_vertex_radius.get());
	}

	void build_hex_mesh()
	{
		Scalar min_radius = std::numeric_limits<Scalar>::max();
		for (Scalar r : *graph_vertex_radius_)
			if (r < min_radius)
				min_radius = r;

		auto radius_copy = add_attribute<Scalar, GraphVertex>(*graph_, "radius_copy");
		radius_copy->copy(graph_vertex_radius_.get());
		graph_vertex_radius_->fill(min_radius);

		contact_surface_ = surface_provider_->add_mesh("contact");
		volume_ = volume_provider_->add_mesh("hex");

		hex_building_attributes_ = modeling::graph_to_hex(*graph_, *contact_surface_, *volume_);

		graph_vertex_radius_->swap(radius_copy.get());
		remove_attribute<GraphVertex>(*graph_, radius_copy);

		surface_provider_->emit_connectivity_changed(contact_surface_);
		volume_provider_->emit_connectivity_changed(volume_);

		volume_vertex_position_ = get_attribute<Vec3, VolumeVertex>(*volume_, "position");
	}

	void test_ortho()
	{
		Scalar min_radius = std::numeric_limits<Scalar>::max();
		for (Scalar r : *graph_vertex_radius_)
			if (r < min_radius)
				min_radius = r;

		auto radius_copy = add_attribute<Scalar, GraphVertex>(*graph_, "radius_copy");
		radius_copy->copy(graph_vertex_radius_.get());
		graph_vertex_radius_->fill(min_radius);

		contact_surface_ = surface_provider_->add_mesh("contact");
		contact_surface_2 = surface_provider_->add_mesh("contact2");
		volume_ = volume_provider_->add_mesh("hex");
		std::cout << "indexed volume: " << is_indexed<CMap3::Volume>(*volume_) << std::endl;
		modeling::create_ortho_hex(*graph_, *contact_surface_, *contact_surface_2, *volume_);

		graph_vertex_radius_->swap(radius_copy.get());
		remove_attribute<GraphVertex>(*graph_, radius_copy);

		surface_provider_->emit_connectivity_changed(contact_surface_);
		surface_provider_->emit_connectivity_changed(contact_surface_2);
		volume_provider_->emit_connectivity_changed(volume_);

		// volume_vertex_position_ = get_attribute<Vec3, VolumeVertex>(*volume_, "position");
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

		using SelectedFace = std::tuple<SurfaceFace, Vec3, Scalar>;

		auto normal_filtered = add_attribute<Vec3, SurfaceVertex>(volume_skin, "normal_filtered");

		geometry::compute_normal(volume_skin, volume_skin_vertex_position.get(), volume_skin_vertex_normal.get());
		geometry::filter_average<Vec3>(volume_skin, volume_skin_vertex_normal.get(), normal_filtered.get());
		volume_skin_vertex_normal->swap(normal_filtered.get());

		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			Vec3& p = value<Vec3>(volume_skin, volume_skin_vertex_position, v);
			const Vec3& n = value<Vec3>(volume_skin, volume_skin_vertex_normal, v);
			std::vector<SelectedFace> selectedfaces =
				geometry::internal::picking(*surface_, surface_vertex_position_.get(), p, p + n);
			Vec3 pos = selectedfaces.size() > 0 ? std::get<1>(selectedfaces[0]) : p;
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v)) = pos;
			p = pos;
			return true;
		});

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void subdivide_volume()
	{
		CellMarker<VOLUME, VolumeFace> cm(*volume_);
		modeling::mark_tranversal_faces(*volume_, *contact_surface_, std::get<1>(hex_building_attributes_), cm);
		modeling::trisect_length_wise(*volume_, std::get<2>(hex_building_attributes_), cm, *graph_,
									  std::get<0>(hex_building_attributes_));
		modeling::subdivide_width_wise(*volume_, std::get<2>(hex_building_attributes_), cm, *graph_,
									   std::get<0>(hex_building_attributes_));
		modeling::subdivide_length_wise(*volume_, std::get<2>(hex_building_attributes_), cm, *graph_,
										std::get<0>(hex_building_attributes_));
		modeling::subdivide_width_wise(*volume_, std::get<2>(hex_building_attributes_), cm, *graph_,
									   std::get<0>(hex_building_attributes_));
		// modeling::subdivide_length_wise(*volume_, std::get<2>(hex_building_attributes_), cm, *graph_,
		// 								std::get<0>(hex_building_attributes_));
		// modeling::subdivide_width_wise(*volume_, std::get<2>(hex_building_attributes_), cm, *graph_,
		// 							   std::get<0>(hex_building_attributes_));

		graph_provider_->emit_connectivity_changed(graph_);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());

		volume_provider_->emit_connectivity_changed(volume_);
		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
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

		// auto position2 = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position2");
		// geometry::filter_average<Vec3>(volume_skin, volume_skin_vertex_position.get(), position2.get());
		// volume_skin_vertex_position->swap(position2.get());

		foreach_cell(volume_skin, [&](SurfaceVertex v) -> bool {
			Vec3 p = geometry::closest_point_on_surface(*surface_, surface_vertex_position_.get(),
														value<Vec3>(volume_skin, volume_skin_vertex_position, v));
			// value<Vec3>(volume_skin, volume_skin_vertex_position, v) = p;
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v)) = p;
			return true;
		});

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void regularize_surface_vertices()
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
		Scalar fit_to_data = 5.0;

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
			Vec3 p = geometry::closest_point_on_surface(*surface_, surface_vertex_position_.get(), *surface_grid_,
														value<Vec3>(volume_skin, volume_skin_vertex_position, v));
			// value<Vec3>(volume_skin, volume_skin_vertex_position, v) = p;
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(volume_skin, volume_skin_vertex_volume_vertex, v)) = p;
			return true;
		});

		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void optimize_volume_vertices()
	{
		auto vertex_index = add_attribute<uint32, VolumeVertex>(*volume_, "__vertex_index");

		uint32 nb_vertices = 0;
		foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			value<uint32>(*volume_, vertex_index, v) = nb_vertices++;
			return true;
		});
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(nb_vertices, nb_vertices);
		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_vertices * 10);
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

		// set constrained vertices
		foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			if (is_incident_to_boundary(*volume_, v))
			{
				int idx = int(value<uint32>(*volume_, vertex_index, v));
				A.prune([&](int i, int, Scalar) { return i != idx; });
				A.coeffRef(idx, idx) = 1.0;
			}
			return true;
		});
		A.makeCompressed();

		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);
		// Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(LAPL);

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(nb_vertices, 3);

		parallel_foreach_cell(*volume_, [&](VolumeVertex v) -> bool {
			uint32 vidx = value<uint32>(*volume_, vertex_index, v);
			const Vec3& pos = value<Vec3>(*volume_, volume_vertex_position_, v);
			if (is_incident_to_boundary(*volume_, v))
			{
				b.coeffRef(vidx, 0) = pos[0];
				b.coeffRef(vidx, 1) = pos[1];
				b.coeffRef(vidx, 2) = pos[2];
			}
			else
			{
				b.coeffRef(vidx, 0) = 0;
				b.coeffRef(vidx, 1) = 0;
				b.coeffRef(vidx, 2) = 0;
			}
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

		remove_attribute<VolumeVolume>(*volume_, hex_frame.get());
		remove_attribute<VolumeVertex2>(*volume_, corner_frame.get());
	}

	void export_subdivided_skin()
	{
		SURFACE volume_skin;
		auto volume_skin_vertex_position = add_attribute<Vec3, SurfaceVertex>(volume_skin, "position");

		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), volume_skin,
										 volume_skin_vertex_position.get());
		modeling::catmull_clark_approx(volume_skin, volume_skin_vertex_position.get(), 2);
		geometry::apply_ear_triangulation(volume_skin, volume_skin_vertex_position.get());
		modeling::export_surface_off(volume_skin, volume_skin_vertex_position.get(), "surface.off");
	}

	void interface() override
	{
		ImGui::TextUnformatted("Graph");
		imgui_mesh_selector(graph_provider_, graph_, "Graph", [&](GRAPH* g) {
			graph_ = g;
			graph_vertex_position_ = nullptr;
			graph_vertex_radius_ = nullptr;
		});

		if (graph_)
		{
			imgui_combo_attribute<GraphVertex, Vec3>(
				*graph_, graph_vertex_position_, "Position##graph",
				[&](const std::shared_ptr<GraphAttribute<Vec3>>& attribute) { graph_vertex_position_ = attribute; });
			imgui_combo_attribute<GraphVertex, Scalar>(
				*graph_, graph_vertex_radius_, "Radius##graph",
				[&](const std::shared_ptr<GraphAttribute<Scalar>>& attribute) { graph_vertex_radius_ = attribute; });

			if (ImGui::Button("Create radius attribute"))
				graph_vertex_radius_ = add_attribute<Scalar, GraphVertex>(*graph_, "radius");
		}

		ImGui::Separator();
		ImGui::TextUnformatted("Surface");
		imgui_mesh_selector(surface_provider_, surface_, "Surface", [&](SURFACE* s) {
			surface_ = s;
			surface_vertex_position_ = nullptr;
		});

		if (surface_)
		{
			imgui_combo_attribute<SurfaceVertex, Vec3>(*surface_, surface_vertex_position_, "Position##surface",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   surface_vertex_position_ = attribute;
														   if (surface_grid_)
															   delete surface_grid_;
														   surface_grid_ =
															   new Grid(*surface_, surface_vertex_position_);
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
			if (ImGui::Button("Init radius from surface"))
				init_graph_radius_from_surface();
			if (ImGui::Button("Extend graph extremities"))
				extend_graph_extremities();
		}
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_)
		{
			if (ImGui::Button("Resample graph"))
				resample_graph();
			if (ImGui::Button("Build hex mesh"))
				build_hex_mesh();
			if (ImGui::Button("Test Ortho"))
				test_ortho();
		}
		if (volume_)
		{
			if (ImGui::Button("Export subdivided skin"))
				export_subdivided_skin();
			if (ImGui::Button("Subdivide volume"))
				subdivide_volume();
			if (ImGui::Button("Project on surface"))
				project_on_surface();
			if (ImGui::Button("Smooth volume surface"))
				smooth_volume_surface();
			if (ImGui::Button("Regularize surface vertices"))
				regularize_surface_vertices();
			if (ImGui::Button("Optimize volume vertices"))
				optimize_volume_vertices();
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
	Grid* surface_grid_;

	SURFACE* contact_surface_;
	SURFACE* contact_surface_2;

	VOLUME* volume_;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_vertex_position_;

	std::tuple<modeling::GAttributes, modeling::M2Attributes, modeling::M3Attributes> hex_building_attributes_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_TUBULAR_MESH_H_
