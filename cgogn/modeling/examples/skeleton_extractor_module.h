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

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>

#include <cgogn/modeling/algos/medial_axis.h>
#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>
#include <cgogn/modeling/algos/skeleton.h>

#include <queue>

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

	void medial_axis(SURFACE& m, SurfaceAttribute<Vec3>* vertex_position, SurfaceAttribute<Vec3>* vertex_normal)
	{
		auto sbc = get_attribute<Vec3, SurfaceVertex>(m, "shrinking_ball_centers");
		if (!sbc)
			sbc = add_attribute<Vec3, SurfaceVertex>(m, "shrinking_ball_centers");
		modeling::shrinking_ball_centers(m, vertex_position, vertex_normal, sbc.get());
	}

	void skeletonize(SURFACE& m, SurfaceAttribute<Vec3>* vertex_position, SurfaceAttribute<Vec3>* vertex_normal)
	{
		modeling::pliant_remeshing(m, vertex_position);
		geometry::compute_normal(m, vertex_position, vertex_normal);

		auto sbc = get_attribute<Vec3, SurfaceVertex>(m, "shrinking_ball_centers");
		if (!sbc)
			sbc = add_attribute<Vec3, SurfaceVertex>(m, "shrinking_ball_centers");
		modeling::shrinking_ball_centers(m, vertex_position, vertex_normal, sbc.get());

		modeling::mean_curvature_skeleton(m, vertex_position, sbc.get());

		surface_provider_->emit_connectivity_changed(&m);
		surface_provider_->emit_attribute_changed(&m, vertex_position);
	}

	void graph_from_surface(SURFACE& s, SurfaceAttribute<Vec3>* surface_vertex_position, GRAPH& g,
							GraphAttribute<Vec3>* graph_vertex_position)
	{
		using WeightedEdge = std::tuple<Scalar, SurfaceEdge>;
		std::priority_queue<WeightedEdge> queue;
		foreach_cell(s, [&](SurfaceEdge e) -> bool {
			queue.emplace(geometry::length(s, e, surface_vertex_position), e);
			return true;
		});

		// GRAPH* g = graph_provider_->add_mesh("extracted_graph");
		// foreach_cell(s, [&](SurfaceVertex v) -> bool {
		// 	add_vertex();
		// 	return true;
		// });
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
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface", [&](SURFACE* m) {
			selected_surface_ = m;
			selected_surface_vertex_position_.reset();
			selected_surface_vertex_normal_.reset();
			surface_provider_->mesh_data(selected_surface_)->outlined_until_ = App::frame_time_ + 1.0;
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
				if (selected_surface_vertex_normal_)
				{
					if (ImGui::Button("Medial axis"))
						medial_axis(*selected_surface_, selected_surface_vertex_position_.get(),
									selected_surface_vertex_normal_.get());
					if (ImGui::Button("Skeletonize"))
						skeletonize(*selected_surface_, selected_surface_vertex_position_.get(),
									selected_surface_vertex_normal_.get());
				}
			}
		}
	}

private:
	MeshProvider<GRAPH>* graph_provider_;
	MeshProvider<SURFACE>* surface_provider_;

	GRAPH* graph_;
	std::shared_ptr<GraphAttribute<Vec3>> graph_vertex_position_;

	SURFACE* selected_surface_;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_position_;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_normal_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SKELETON_EXTRACTOR_H_
