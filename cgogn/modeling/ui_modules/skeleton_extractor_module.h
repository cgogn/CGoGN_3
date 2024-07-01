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

#ifndef CGOGN_MODULE_SKELETON_EXTRACTOR_H_
#define CGOGN_MODULE_SKELETON_EXTRACTOR_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/convert.h>

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>

#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>
#include <cgogn/modeling/algos/skeleton.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Mat3;
using geometry::Scalar;
using geometry::Vec3;

template <typename SURFACE, typename NONMANIFOLD, typename GRAPH>
class SkeletonExtractor : public Module
{
	template <typename T>
	using GraphAttribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	template <typename T>
	using NonManifoldAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NonManifoldFace = typename mesh_traits<NONMANIFOLD>::Face;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	struct SurfaceParameters
	{
		SURFACE* mesh_;
		std::shared_ptr<SurfaceAttribute<Vec3>> medial_axis_samples_position_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> medial_axis_samples_radius_ = nullptr;
		std::shared_ptr<SurfaceAttribute<SurfaceVertex>> medial_axis_samples_secondary_vertex_ = nullptr;

		float32 radius_threshold_ = 0.025f;
		float32 angle_threshold_ = M_PI / 2.0f;

		CellsSet<SURFACE, SurfaceVertex>* filtered_medial_axis_samples_set_ = nullptr;
		CellsSet<SURFACE, SurfaceVertex>* neutralized_surface_vertices_set_ = nullptr;
	};

public:
	SkeletonExtractor(const App& app)
		: Module(app, "SkeletonExtractor"), non_manifold_(nullptr), non_manifold_vertex_position_(nullptr),
		  selected_surface_(nullptr), selected_surface_vertex_position_(nullptr),
		  selected_surface_vertex_normal_(nullptr)
	{
	}

	~SkeletonExtractor()
	{
	}

private:
	void init_surface_mesh(SURFACE* s)
	{
		SurfaceParameters& p = surface_parameters_[s];
		p.mesh_ = s;
		p.medial_axis_samples_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(*s, "medial_axis_samples_position");
		p.medial_axis_samples_secondary_vertex_ =
			get_or_add_attribute<SurfaceVertex, SurfaceVertex>(*s, "medial_axis_samples_secondary_vertex");
		p.medial_axis_samples_radius_ = get_or_add_attribute<Scalar, SurfaceVertex>(*s, "medial_axis_samples_radius");

		MeshData<SURFACE>& md = surface_provider_->mesh_data(*s);
		p.filtered_medial_axis_samples_set_ =
			&md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
		p.neutralized_surface_vertices_set_ =
			&md.template get_or_add_cells_set<SurfaceVertex>("neutralized_surface_vertices");
	}

public:
	void sample_medial_axis(SURFACE& s, SurfaceAttribute<Vec3>* vertex_position, SurfaceAttribute<Vec3>* vertex_normal)
	{
		SurfaceParameters& p = surface_parameters_[&s];

		geometry::shrinking_ball_centers(s, vertex_position, vertex_normal, p.medial_axis_samples_position_.get(),
										 p.medial_axis_samples_radius_.get(),
										 p.medial_axis_samples_secondary_vertex_.get());

		p.filtered_medial_axis_samples_set_->select_if([&](SurfaceVertex v) { return true; });

		surface_provider_->emit_attribute_changed(s, p.medial_axis_samples_position_.get());
		surface_provider_->emit_attribute_changed(s, p.medial_axis_samples_radius_.get());

		surface_provider_->emit_cells_set_changed(s, p.filtered_medial_axis_samples_set_);
	}

	void filter_medial_axis_samples(SURFACE& s, SurfaceAttribute<Vec3>* vertex_position)
	{
		SurfaceParameters& p = surface_parameters_[&s];

		p.filtered_medial_axis_samples_set_->clear();
		p.neutralized_surface_vertices_set_->clear();
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			const Vec3& c = value<Vec3>(s, p.medial_axis_samples_position_, v);
			const SurfaceVertex sv = value<SurfaceVertex>(s, p.medial_axis_samples_secondary_vertex_, v);
			const Vec3& c1 = value<Vec3>(s, vertex_position, v);
			const Vec3& c2 = value<Vec3>(s, vertex_position, sv);
			const Scalar r = value<Scalar>(s, p.medial_axis_samples_radius_, v);
			if (r > p.radius_threshold_ && geometry::angle(c1 - c, c2 - c) > p.angle_threshold_)
				p.filtered_medial_axis_samples_set_->select(v);
			else
				p.neutralized_surface_vertices_set_->select(v);
			return true;
		});

		surface_provider_->emit_cells_set_changed(s, p.filtered_medial_axis_samples_set_);
	}

	void skeletonize(SURFACE& s, std::shared_ptr<SurfaceAttribute<Vec3>>& vertex_position, Scalar wL, Scalar wH,
					 Scalar wM, Scalar resampling_ratio)
	{
		modeling::mean_curvature_skeleton(s, vertex_position, wL, wH, wM, resampling_ratio);

		surface_provider_->emit_connectivity_changed(s);
		surface_provider_->emit_attribute_changed(s, vertex_position.get());
	}

	void create_non_manifold_from_surface(SURFACE& s, SurfaceAttribute<Vec3>* surface_vertex_position)
	{
		non_manifold_ = non_manifold_provider_->add_mesh("extracted_non_manifold");
		non_manifold_vertex_position_ = add_attribute<Vec3, NonManifoldVertex>(*non_manifold_, "position");

		non_manifold_from_surface(
			s, *non_manifold_,
			[&](NonManifoldVertex nmv, SurfaceVertex sv) {
				value<Vec3>(*non_manifold_, non_manifold_vertex_position_, nmv) =
					value<Vec3>(s, surface_vertex_position, sv);
			},
			[](NonManifoldEdge, SurfaceEdge) {}, [](NonManifoldFace, SurfaceFace) {});

		non_manifold_provider_->emit_connectivity_changed(*non_manifold_);
		non_manifold_provider_->emit_attribute_changed(*non_manifold_, non_manifold_vertex_position_.get());
	}

	void create_graph_from_non_manifold(NONMANIFOLD& nm, NonManifoldAttribute<Vec3>* non_manifold_vertex_position)
	{
		graph_ = graph_provider_->add_mesh("extracted_graph");
		graph_vertex_position_ = add_attribute<Vec3, GraphVertex>(*graph_, "position");

		graph_from_non_manifold(
			nm, *graph_,
			[&](GraphVertex gv, NonManifoldVertex nmv) {
				value<Vec3>(*graph_, graph_vertex_position_, gv) = value<Vec3>(nm, non_manifold_vertex_position, nmv);
			},
			[](GraphEdge, NonManifoldEdge) {});

		graph_provider_->emit_connectivity_changed(*graph_);
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());
	}

	void collapse_non_manifold()
	{
		using EdgeQueue = std::multimap<Scalar, NonManifoldEdge>;
		using EdgeQueueIt = typename EdgeQueue::const_iterator;
		using EdgeInfo = std::pair<bool, EdgeQueueIt>; // {valid, iterator}

		EdgeQueue queue;
		auto edge_queue_it = add_attribute<EdgeInfo, NonManifoldEdge>(*non_manifold_, "__non_manifold_edge_queue_it");
		foreach_cell(*non_manifold_, [&](NonManifoldEdge e) -> bool {
			value<EdgeInfo>(*non_manifold_, edge_queue_it, e) = {
				true, queue.emplace(geometry::length(*non_manifold_, e, non_manifold_vertex_position_.get()), e)};
			return true;
		});

		using PositionAccu = std::vector<Vec3>;
		auto vertex_position_accu =
			add_attribute<PositionAccu, NonManifoldVertex>(*non_manifold_, "__non_manifold_vertex_position_accu");
		foreach_cell(*non_manifold_, [&](NonManifoldVertex v) -> bool {
			value<PositionAccu>(*non_manifold_, vertex_position_accu,
								v) = {value<Vec3>(*non_manifold_, non_manifold_vertex_position_, v)};
			return true;
		});

		while (!queue.empty())
		{
			auto it = queue.begin();
			NonManifoldEdge e = (*it).second;
			queue.erase(it);
			value<EdgeInfo>(*non_manifold_, edge_queue_it, e).first = false;

			std::vector<NonManifoldFace> ifaces = incident_faces(*non_manifold_, e);
			if (ifaces.size() == 0)
				continue;

			// iv[0] will be removed and iv[1] will survive
			std::vector<NonManifoldVertex> iv = incident_vertices(*non_manifold_, e);
			PositionAccu& accu0 = value<PositionAccu>(*non_manifold_, vertex_position_accu, iv[0]);
			PositionAccu& accu1 = value<PositionAccu>(*non_manifold_, vertex_position_accu, iv[1]);
			accu1.insert(accu1.end(), accu0.begin(), accu0.end());

			auto [v, removed_edges] = collapse_edge(*non_manifold_, e);
			for (NonManifoldEdge re : removed_edges)
			{
				EdgeInfo einfo = value<EdgeInfo>(*non_manifold_, edge_queue_it, re);
				if (einfo.first)
					queue.erase(einfo.second);
			}

			foreach_incident_edge(*non_manifold_, v, [&](NonManifoldEdge ie) -> bool {
				EdgeInfo einfo = value<EdgeInfo>(*non_manifold_, edge_queue_it, ie);
				if (einfo.first)
					queue.erase(einfo.second);
				value<EdgeInfo>(*non_manifold_, edge_queue_it, ie) = {
					true, queue.emplace(geometry::length(*non_manifold_, ie, non_manifold_vertex_position_.get()), ie)};
				return true;
			});
		}

		foreach_cell(*non_manifold_, [&](NonManifoldVertex v) -> bool {
			Vec3 mean{0, 0, 0};
			uint32 count = 0;
			for (Vec3& p : value<PositionAccu>(*non_manifold_, vertex_position_accu, v))
			{
				mean += p;
				++count;
			}
			mean /= count;
			value<Vec3>(*non_manifold_, non_manifold_vertex_position_, v) = mean;
			return true;
		});

		remove_attribute<NonManifoldEdge>(*non_manifold_, edge_queue_it);
		remove_attribute<NonManifoldVertex>(*non_manifold_, vertex_position_accu);

		non_manifold_provider_->emit_connectivity_changed(*non_manifold_);
		non_manifold_provider_->emit_attribute_changed(*non_manifold_, non_manifold_vertex_position_.get());
	}

protected:
	void init() override
	{
		graph_provider_ = static_cast<ui::MeshProvider<GRAPH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));
		non_manifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		surface_provider_->foreach_mesh([this](SURFACE& m, const std::string&) { init_surface_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<SURFACE>::mesh_added>(
			surface_provider_, this, &SkeletonExtractor<SURFACE, NONMANIFOLD, GRAPH>::init_surface_mesh));
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface", [&](SURFACE& m) {
			selected_surface_ = &m;
			selected_surface_vertex_position_.reset();
			selected_surface_vertex_normal_.reset();
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_surface_)
		{
			imgui_combo_attribute<SurfaceVertex, Vec3>(*selected_surface_, selected_surface_vertex_position_,
													   "Position",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   selected_surface_vertex_position_ = attribute;
													   });

			imgui_combo_attribute<SurfaceVertex, Vec3>(*selected_surface_, selected_surface_vertex_normal_, "Normal",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   selected_surface_vertex_normal_ = attribute;
													   });

			if (selected_surface_vertex_position_)
			{
				static float wL = 0.6f, wH = 0.1f, wM = 0.1f, resampling_ratio = 0.9f;
				ImGui::SliderFloat("Smoothness", &wL, 0.01f, 1.0f);
				ImGui::SliderFloat("Velocity", &wH, 0.01f, 1.0f);
				ImGui::SliderFloat("Medial attraction", &wM, 0.01f, 1.0f);
				ImGui::SliderFloat("Resampling ratio", &resampling_ratio, 0.01f, 2.0f);
				if (ImGui::Button("Skeletonize"))
					skeletonize(*selected_surface_, selected_surface_vertex_position_, wL, wH, wM, resampling_ratio);
				if (ImGui::Button("Non-manifold from surface"))
					create_non_manifold_from_surface(*selected_surface_, selected_surface_vertex_position_.get());

				if (selected_surface_vertex_normal_)
				{
					SurfaceParameters& p = surface_parameters_[selected_surface_];

					ImGui::Separator();
					if (ImGui::Button("Sample medial axis"))
						sample_medial_axis(*selected_surface_, selected_surface_vertex_position_.get(),
										   selected_surface_vertex_normal_.get());
					if (ImGui::SliderFloat("Min radius (log)", &p.radius_threshold_, 0.0001f, 1.0f, "%.4f"))
						filter_medial_axis_samples(*selected_surface_, selected_surface_vertex_position_.get());
					if (ImGui::SliderFloat("Min angle (log)", &p.angle_threshold_, 0.0001f, M_PI, "%.4f"))
						filter_medial_axis_samples(*selected_surface_, selected_surface_vertex_position_.get());
					if (ImGui::Button("Filter medial axis samples"))
						filter_medial_axis_samples(*selected_surface_, selected_surface_vertex_position_.get());
				}
			}
			if (non_manifold_)
			{
				ImGui::Separator();
				if (ImGui::Button("Collapse non-manifold"))
					collapse_non_manifold();
				if (ImGui::Button("Graph from non-manifold"))
					create_graph_from_non_manifold(*non_manifold_, non_manifold_vertex_position_.get());
			}
		}
	}

private:
	MeshProvider<GRAPH>* graph_provider_ = nullptr;
	MeshProvider<NONMANIFOLD>* non_manifold_provider_ = nullptr;
	MeshProvider<SURFACE>* surface_provider_ = nullptr;

	GRAPH* graph_ = nullptr;
	std::shared_ptr<GraphAttribute<Vec3>> graph_vertex_position_ = nullptr;

	NONMANIFOLD* non_manifold_ = nullptr;
	std::shared_ptr<NonManifoldAttribute<Vec3>> non_manifold_vertex_position_ = nullptr;

	SURFACE* selected_surface_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_position_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_normal_ = nullptr;
	CellsSet<SURFACE, SurfaceVertex>* selected_surface_filtered_medial_axis_samples_set_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SKELETON_EXTRACTOR_H_
