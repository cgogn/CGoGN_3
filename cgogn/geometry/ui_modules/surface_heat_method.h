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

#ifndef CGOGN_MODULE_SURFACE_HEAT_METHOD_H_
#define CGOGN_MODULE_SURFACE_HEAT_METHOD_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/distance.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_set>
#include <utility>

namespace cgogn
{

namespace ui
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

template <typename MESH>
class SurfaceHeatMethod : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SurfaceHeatMethod can only be used with meshes of dimension 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using Edge = typename mesh_traits<MESH>::Edge;

public:
	SurfaceHeatMethod(const App& app) : Module(app, "SurfaceHeatMethod (" + std::string{mesh_traits<MESH>::name} + ")")
	{
	}

	~SurfaceHeatMethod()
	{
	}

	void euclidean_distance(const MESH& m, const Attribute<Vec3>* vertex_position,
							const CellsSet<MESH, Vertex>* source_vertices, Attribute<Scalar>* vertex_euclidean_dist)
	{
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Scalar s = std::numeric_limits<Scalar>::max();
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			source_vertices->foreach_cell([&](Vertex vec) {
				Scalar d = (pos - value<Vec3>(m, vertex_position, vec)).norm();
				if (s > d)
					s = d;
			});
			value<Scalar>(m, vertex_euclidean_dist, v) = s;
			return true;
		});

		mesh_provider_->emit_attribute_changed(m, vertex_euclidean_dist);
	}

	void geodesic_distance(MESH& m, const Attribute<Vec3>* vertex_position,
						   const CellsSet<MESH, Vertex>* source_vertices, Attribute<Scalar>* vertex_geodesic_distance,
						   double t_multiplier)
	{
		geometry::compute_geodesic_distance(m, vertex_position, source_vertices, vertex_geodesic_distance,
											t_multiplier);
		mesh_provider_->emit_attribute_changed(m, vertex_geodesic_distance);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void left_panel() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			selected_vertex_position_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Vertex Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			imgui_combo_cells_set(md, selected_vertices_set_, "Source vertices",
								  [&](CellsSet<MESH, Vertex>* cs) { selected_vertices_set_ = cs; });

			if (selected_vertex_position_ && selected_vertices_set_)
			{
				imgui_combo_attribute<Vertex, Scalar>(*selected_mesh_, euclidean_distance_vertex_,
													  "Vertex Euclidean distance",
													  [&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
														  euclidean_distance_vertex_ = attribute;
													  });

				if (ImGui::Button("Compute Euclidean Distance"))
				{
					if (!euclidean_distance_vertex_)
						euclidean_distance_vertex_ =
							get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "euclidean distance");
					euclidean_distance(*selected_mesh_, selected_vertex_position_.get(), selected_vertices_set_,
									   euclidean_distance_vertex_.get());
				}

				imgui_combo_attribute<Vertex, Scalar>(*selected_mesh_, geodesic_distance_vertex_,
													  "Vertex Geodesic distance",
													  [&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
														  geodesic_distance_vertex_ = attribute;
													  });

				static double t_multiplier = 1.0;
				ImGui::InputDouble("Scalar t_multiplier", &t_multiplier, 0.01f, 100.0f, "%.3f");

				if (ImGui::Button("Compute Geodesic Distance"))
				{
					if (!geodesic_distance_vertex_)
						geodesic_distance_vertex_ =
							get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "geodesic distance");
					geodesic_distance(*selected_mesh_, selected_vertex_position_.get(), selected_vertices_set_,
									  geodesic_distance_vertex_.get(), t_multiplier);
				}
			}
		}
	}

private:
	MESH* selected_mesh_ = nullptr;

	// geometry::DistanceHeatSolver<MESH>* distanceHeatSolver = nullptr;

	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_ = nullptr;
	CellsSet<MESH, Vertex>* selected_vertices_set_ = nullptr;

	MeshProvider<MESH>* mesh_provider_ = nullptr;

	std::shared_ptr<Attribute<Scalar>> euclidean_distance_vertex_ = nullptr;
	std::shared_ptr<Attribute<Scalar>> geodesic_distance_vertex_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_HEAT_METHOD_H_
