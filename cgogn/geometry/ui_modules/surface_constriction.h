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

#ifndef CGOGN_MODULE_SURFACE_CONSTRICTION_H_
#define CGOGN_MODULE_SURFACE_CONSTRICTION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/curvature.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceConstriction : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2,
				  "SurfaceConstriction can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Volume = typename mesh_traits<MESH>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceConstriction(const App& app)
		: ViewModule(app, "SurfaceConstriction (" + std::string{mesh_traits<MESH>::name} + ")"),
		selected_mesh_(nullptr), selected_vertex_position_(nullptr), selected_vertex_normal_(nullptr),
		selected_vertex_kgauss_(nullptr)
	{
	}

	~SurfaceConstriction()
	{
	}

	void compute_vertex_area(const MESH& m, const Attribute<Vec3>* vertex_position, Attribute<Scalar>* vertex_area)
	{
		geometry::compute_area<Vertex>(m, vertex_position, vertex_area);
		mesh_provider_->emit_attribute_changed(m, vertex_area);
	}

	void compute_normal(const MESH& m, const Attribute<Vec3>* vertex_position, Attribute<Vec3>* vertex_normal)
	{
		geometry::compute_normal<Vertex>(m, vertex_position, vertex_normal);
		mesh_provider_->emit_attribute_changed(m, vertex_normal);
	}

	void compute_gaussian_curvature(const MESH& m, Scalar radius, const Attribute<Vec3>* vertex_position,
						   const Attribute<Vec3>* vertex_normal, const Attribute<Scalar>* edge_angle,
						   Attribute<Scalar>* vertex_kgauss)
	{
		geometry::compute_gaussian_curvature(m, radius, vertex_position, vertex_normal, edge_angle, vertex_kgauss);
		mesh_provider_->emit_attribute_changed(m, vertex_kgauss);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void draw(View* view) override
	{
	}

	void left_panel() override
	{
		bool need_update = false;

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_normal_, "Normal",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_normal_ = attribute; });

			if (selected_vertex_position_)
			{
				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

				if (ImGui::Button("Compute area"))
				{
					if (!selected_vertex_area_)
						selected_vertex_area_ = get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "area");
					compute_vertex_area(*selected_mesh_, selected_vertex_position_.get(), selected_vertex_area_.get());
				}

				if (ImGui::Button("Compute normal"))
				{
					if (!selected_vertex_normal_)
						selected_vertex_normal_ = get_or_add_attribute<Vec3, Vertex>(*selected_mesh_, "normal");
					compute_normal(*selected_mesh_, selected_vertex_position_.get(), selected_vertex_normal_.get());
				}
			}

			if (selected_vertex_position_ && selected_vertex_normal_)
			{
				if (ImGui::Button("Compute gaussian curvature"))
				{
					if (!selected_vertex_kgauss_)
						selected_vertex_kgauss_ = get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "kgauss");

					std::shared_ptr<Attribute<Scalar>> edge_angle =
						add_attribute<Scalar, Edge>(*selected_mesh_, "__edge_angle");
					geometry::compute_angle(*selected_mesh_, selected_vertex_position_.get(), edge_angle.get());	// TODO different centroid

					Scalar mean_edge_length =
						geometry::mean_edge_length(*selected_mesh_, selected_vertex_position_.get()); // TODO different neighborhoud

					compute_gaussian_curvature(*selected_mesh_, mean_edge_length * 2.5, selected_vertex_position_.get(),
									  selected_vertex_normal_.get(), edge_angle.get(), selected_vertex_kgauss_.get());

					remove_attribute<Edge>(*selected_mesh_, edge_angle);
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	MESH* selected_mesh_ = nullptr;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_normal_;
	std::shared_ptr<Attribute<Scalar>> selected_vertex_area_;
	std::shared_ptr<Attribute<Scalar>> selected_vertex_kgauss_;	// gaussian curvature

	MeshProvider<MESH>* mesh_provider_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_CONSTRICTION_H_
