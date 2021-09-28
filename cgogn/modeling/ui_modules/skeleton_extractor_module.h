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

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>

#include <cgogn/modeling/algos/medial_axis.h>
#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>
#include <cgogn/modeling/algos/skeleton.h>

#include <map>

namespace cgogn
{

namespace ui
{

template <typename GRAPH, typename SURFACE>
class SkeletonExtractor : public Module
{
	template <typename T>
	using GraphAttribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;
	using GraphFace = typename mesh_traits<GRAPH>::Face;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;
	using Mat3 = geometry::Mat3;

public:
	SkeletonExtractor(const App& app)
		: Module(app, "SkeletonExtractor"), graph_(nullptr), graph_vertex_position_(nullptr),
		  selected_surface_(nullptr), selected_surface_vertex_position_(nullptr),
		  selected_surface_vertex_normal_(nullptr)
	{
	}

	~SkeletonExtractor()
	{
	}

	void medial_axis(SURFACE& s, SurfaceAttribute<Vec3>* vertex_position, SurfaceAttribute<Vec3>* vertex_normal)
	{
		auto sbc = get_attribute<Vec3, SurfaceVertex>(s, "shrinking_ball_centers");
		if (!sbc)
			sbc = add_attribute<Vec3, SurfaceVertex>(s, "shrinking_ball_centers");
		modeling::shrinking_ball_centers(s, vertex_position, vertex_normal, sbc.get());
	}

	void skeletonize(SURFACE& s, std::shared_ptr<SurfaceAttribute<Vec3>>& vertex_position, Scalar wL, Scalar wH,
					 Scalar wM)
	{
		modeling::mean_curvature_skeleton(s, vertex_position, wL, wH, wM);

		surface_provider_->emit_connectivity_changed(s);
		surface_provider_->emit_attribute_changed(s, vertex_position.get());
	}

	void graph_from_surface(SURFACE& s, SurfaceAttribute<Vec3>* surface_vertex_position)
	{
		graph_ = graph_provider_->add_mesh("extracted_graph");
		graph_vertex_position_ = add_attribute<Vec3, GraphVertex>(*graph_, "position");

		auto surface_graph_vertex = add_attribute<GraphVertex, SurfaceVertex>(s, "__graph_vertex");
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			GraphVertex gv = add_vertex(*graph_);
			value<GraphVertex>(s, surface_graph_vertex, v) = gv;
			value<Vec3>(*graph_, graph_vertex_position_, gv) = value<Vec3>(s, surface_vertex_position, v);
			return true;
		});
		auto surface_graph_edge = add_attribute<GraphEdge, SurfaceEdge>(s, "__graph_edge");
		foreach_cell(s, [&](SurfaceEdge e) -> bool {
			std::vector<SurfaceVertex> iv = incident_vertices(s, e);
			value<GraphEdge>(s, surface_graph_edge, e) =
				add_edge(*graph_, value<GraphVertex>(s, surface_graph_vertex, iv[0]),
						 value<GraphVertex>(s, surface_graph_vertex, iv[1]));
			return true;
		});
		foreach_cell(s, [&](SurfaceFace f) -> bool {
			std::vector<SurfaceEdge> ie = incident_edges(s, f);
			std::vector<GraphEdge> edges;
			std::transform(ie.begin(), ie.end(), std::back_inserter(edges),
						   [&](SurfaceEdge e) { return value<GraphEdge>(s, surface_graph_edge, e); });
			add_face(*graph_, edges);
			return true;
		});
		remove_attribute<SurfaceVertex>(s, surface_graph_vertex);
		remove_attribute<SurfaceEdge>(s, surface_graph_edge);

		graph_provider_->emit_connectivity_changed(*graph_);
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());
	}

	void collapse_graph()
	{
		using EdgeQueue = std::multimap<Scalar, GraphEdge>;
		using EdgeQueueIt = typename EdgeQueue::const_iterator;
		using EdgeInfo = std::pair<bool, EdgeQueueIt>; // {valid, iterator}

		EdgeQueue queue;
		auto edge_queue_it = add_attribute<EdgeInfo, GraphEdge>(*graph_, "__graph_edge_queue_it");
		foreach_cell(*graph_, [&](GraphEdge e) -> bool {
			value<EdgeInfo>(*graph_, edge_queue_it,
							e) = {true, queue.emplace(geometry::length(*graph_, e, graph_vertex_position_.get()), e)};
			return true;
		});

		using PositionAccu = std::vector<Vec3>;
		auto vertex_position_accu = add_attribute<PositionAccu, GraphVertex>(*graph_, "__graph_vertex_position_accu");
		foreach_cell(*graph_, [&](GraphVertex v) -> bool {
			value<PositionAccu>(*graph_, vertex_position_accu, v) = {value<Vec3>(*graph_, graph_vertex_position_, v)};
			return true;
		});

		while (!queue.empty())
		{
			auto it = queue.begin();
			GraphEdge e = (*it).second;
			queue.erase(it);
			value<EdgeInfo>(*graph_, edge_queue_it, e).first = false;

			// if (graph_->attribute_containers_[GraphEdge::CELL_INDEX].nb_refs(e.index_) == 0)
			// 	continue;
			std::vector<GraphFace> ifaces = incident_faces(*graph_, e);
			if (ifaces.size() == 0)
				continue;

			// iv[0] will be removed and iv[1] will survive
			std::vector<GraphVertex> iv = incident_vertices(*graph_, e);
			PositionAccu& accu0 = value<PositionAccu>(*graph_, vertex_position_accu, iv[0]);
			PositionAccu& accu1 = value<PositionAccu>(*graph_, vertex_position_accu, iv[1]);
			accu1.insert(accu1.end(), accu0.begin(), accu0.end());

			auto [v, removed_edges] = collapse_edge(*graph_, e);
			for (GraphEdge re : removed_edges)
			{
				EdgeInfo einfo = value<EdgeInfo>(*graph_, edge_queue_it, re);
				if (einfo.first)
					queue.erase(einfo.second);
			}

			foreach_incident_edge(*graph_, v, [&](GraphEdge ie) -> bool {
				EdgeInfo einfo = value<EdgeInfo>(*graph_, edge_queue_it, ie);
				if (einfo.first)
					queue.erase(einfo.second);
				value<EdgeInfo>(*graph_, edge_queue_it, ie) = {
					true, queue.emplace(geometry::length(*graph_, ie, graph_vertex_position_.get()), ie)};
				return true;
			});
		}

		foreach_cell(*graph_, [&](GraphVertex v) -> bool {
			Vec3 mean{0, 0, 0};
			uint32 count = 0;
			for (Vec3& p : value<PositionAccu>(*graph_, vertex_position_accu, v))
			{
				mean += p;
				++count;
			}
			mean /= count;
			value<Vec3>(*graph_, graph_vertex_position_, v) = mean;
			return true;
		});

		remove_attribute<GraphEdge>(*graph_, edge_queue_it);
		remove_attribute<GraphVertex>(*graph_, vertex_position_accu);

		graph_provider_->emit_connectivity_changed(*graph_);
		graph_provider_->emit_attribute_changed(*graph_, graph_vertex_position_.get());
	}

protected:
	void init() override
	{
		graph_provider_ = static_cast<ui::MeshProvider<GRAPH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
	}

	void interface() override
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
				static float wL = 1.0, wH = 0.1, wM = 0.25;
				ImGui::SliderFloat("Smoothness", &wL, 0.01f, 1.0f);
				ImGui::SliderFloat("Velocity", &wH, 0.01f, 1.0f);
				ImGui::SliderFloat("Medial attraction", &wM, 0.01f, 1.0f);
				if (ImGui::Button("Skeletonize"))
					skeletonize(*selected_surface_, selected_surface_vertex_position_, wL, wH, wM);
				if (ImGui::Button("Graph from surface"))
					graph_from_surface(*selected_surface_, selected_surface_vertex_position_.get());

				if (selected_surface_vertex_normal_)
				{
					if (ImGui::Button("Medial axis"))
						medial_axis(*selected_surface_, selected_surface_vertex_position_.get(),
									selected_surface_vertex_normal_.get());
				}
			}
			if (graph_)
			{
				if (ImGui::Button("Collapse graph"))
					collapse_graph();
			}
		}
	}

private:
	MeshProvider<GRAPH>* graph_provider_ = nullptr;
	MeshProvider<SURFACE>* surface_provider_ = nullptr;

	GRAPH* graph_ = nullptr;
	std::shared_ptr<GraphAttribute<Vec3>> graph_vertex_position_ = nullptr;

	SURFACE* selected_surface_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_position_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_normal_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SKELETON_EXTRACTOR_H_
