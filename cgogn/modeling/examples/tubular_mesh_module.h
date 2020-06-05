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

#include <cgogn/geometry/algos/distance.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/picking.h>

#include <cgogn/modeling/algos/graph_resampling.h>
#include <cgogn/modeling/algos/graph_to_hex.h>

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
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;
	using VolumeVertex = typename mesh_traits<VOLUME>::Vertex;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	TubularMesh(const App& app)
		: ViewModule(app, "TubularMesh"), graph_(nullptr), graph_vertex_position_(nullptr),
		  graph_vertex_radius_(nullptr), surface_(nullptr), surface_vertex_position_(nullptr), volume_(nullptr)
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
		graph_vertex_radius_->fill(l);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
	}

	void init_graph_radius_from_surface()
	{
		modeling::compute_graph_radius_from_surface(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(),
													*surface_, surface_vertex_position_.get());
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());
	}

	void resample_graph()
	{
		resampled_graph_ = graph_provider_->add_mesh("resampled");
		resampled_graph_vertex_position_ = add_attribute<Vec3, GraphVertex>(*resampled_graph_, "position");
		resampled_graph_vertex_radius_ = add_attribute<Scalar, GraphVertex>(*resampled_graph_, "radius");

		modeling::resample_graph(*graph_, graph_vertex_position_.get(), graph_vertex_radius_.get(), *resampled_graph_,
								 resampled_graph_vertex_position_.get(), resampled_graph_vertex_radius_.get());

		graph_provider_->emit_connectivity_changed(graph_);
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(graph_, graph_vertex_radius_.get());

		graph_provider_->emit_connectivity_changed(resampled_graph_);
		graph_provider_->emit_attribute_changed(resampled_graph_, resampled_graph_vertex_position_.get());
		graph_provider_->emit_attribute_changed(resampled_graph_, resampled_graph_vertex_radius_.get());
	}

	void build_hex_mesh()
	{
		Scalar min_radius = std::numeric_limits<Scalar>::max();
		for (Scalar r : *resampled_graph_vertex_radius_)
			if (r < min_radius)
				min_radius = r;
		resampled_graph_vertex_radius_->fill(min_radius);

		contact_surface_ = surface_provider_->add_mesh("contact");
		volume_ = volume_provider_->add_mesh("hex");

		modeling::graph_to_hex(*resampled_graph_, *contact_surface_, *volume_);

		surface_provider_->emit_connectivity_changed(contact_surface_);
		volume_provider_->emit_connectivity_changed(volume_);

		volume_vertex_position_ = get_attribute<Vec3, VolumeVertex>(*volume_, "position");

		volume_skin_ = surface_provider_->add_mesh("skin");
		volume_skin_vertex_position_ = add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "position");
		volume_skin_vertex_normal_ = add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "normal");
		volume_skin_vertex_laplacian_ = add_attribute<Vec3, SurfaceVertex>(*volume_skin_, "laplacian");
		volume_skin_vertex_volume_vertex_ = add_attribute<VolumeVertex, SurfaceVertex>(*volume_skin_, "hex_vertex");

		modeling::extract_volume_surface(*volume_, volume_vertex_position_.get(), *volume_skin_,
										 volume_skin_vertex_position_.get(), volume_skin_vertex_volume_vertex_.get());
		geometry::apply_ear_triangulation(*volume_skin_, volume_skin_vertex_position_.get());

		surface_provider_->emit_connectivity_changed(volume_skin_);
	}

	void project_on_surface()
	{
		using SelectedFace = std::tuple<SurfaceFace, Vec3, Scalar>;

		geometry::compute_normal(*volume_skin_, volume_skin_vertex_position_.get(), volume_skin_vertex_normal_.get());

		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			Vec3& p = value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v);
			const Vec3& n = value<Vec3>(*volume_skin_, volume_skin_vertex_normal_, v);
			std::vector<SelectedFace> selectedfaces =
				geometry::internal::picking(*surface_, surface_vertex_position_.get(), p, p + n);
			Vec3 pos = selectedfaces.size() > 0 ? std::get<1>(selectedfaces[0]) : p;
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v)) = pos;
			p = pos;
			return true;
		});

		surface_provider_->emit_attribute_changed(volume_skin_, volume_skin_vertex_position_.get());
		surface_provider_->emit_attribute_changed(volume_skin_, volume_skin_vertex_normal_.get());
		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void smooth_volume_surface()
	{
		geometry::compute_laplacian(*volume_skin_, volume_skin_vertex_position_.get(),
									volume_skin_vertex_laplacian_.get());
		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v) +=
				0.1 * value<Vec3>(*volume_skin_, volume_skin_vertex_laplacian_, v);
			return true;
		});
		parallel_foreach_cell(*volume_skin_, [&](SurfaceVertex v) -> bool {
			Vec3 p = geometry::closest_point_on_surface(*surface_, surface_vertex_position_.get(),
														value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v));
			value<Vec3>(*volume_skin_, volume_skin_vertex_position_, v) = p;
			value<Vec3>(*volume_, volume_vertex_position_,
						value<VolumeVertex>(*volume_skin_, volume_skin_vertex_volume_vertex_, v)) = p;
			return true;
		});

		surface_provider_->emit_attribute_changed(volume_skin_, volume_skin_vertex_position_.get());
		surface_provider_->emit_attribute_changed(volume_skin_, volume_skin_vertex_laplacian_.get());
		volume_provider_->emit_attribute_changed(volume_, volume_vertex_position_.get());
	}

	void interface() override
	{
		ImGui::TextUnformatted("Graph");
		imgui_mesh_selector(graph_provider_, graph_, "Graph", [&](GRAPH* g) { graph_ = g; });

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
		imgui_mesh_selector(surface_provider_, surface_, "Surface", [&](SURFACE* s) { surface_ = s; });

		if (surface_)
		{
			imgui_combo_attribute<SurfaceVertex, Vec3>(*surface_, surface_vertex_position_, "Position##surface",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   surface_vertex_position_ = attribute;
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
		}
		if (graph_ && graph_vertex_position_ && graph_vertex_radius_)
		{
			if (ImGui::Button("Resample graph"))
				resample_graph();
		}
		if (resampled_graph_ && resampled_graph_vertex_position_ && resampled_graph_vertex_radius_)
		{
			if (ImGui::Button("Build hex mesh"))
				build_hex_mesh();
		}
		if (volume_ && volume_skin_)
		{
			if (ImGui::Button("Project on surface"))
				project_on_surface();
			if (ImGui::Button("Smooth volume surface"))
				smooth_volume_surface();
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

	GRAPH* resampled_graph_;
	std::shared_ptr<GraphAttribute<Vec3>> resampled_graph_vertex_position_;
	std::shared_ptr<GraphAttribute<Scalar>> resampled_graph_vertex_radius_;

	SURFACE* contact_surface_;

	VOLUME* volume_;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_vertex_position_;

	SURFACE* volume_skin_;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_position_;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_normal_;
	std::shared_ptr<SurfaceAttribute<Vec3>> volume_skin_vertex_laplacian_;
	std::shared_ptr<SurfaceAttribute<VolumeVertex>> volume_skin_vertex_volume_vertex_;
}; // namespace ui

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_TUBULAR_MESH_H_
