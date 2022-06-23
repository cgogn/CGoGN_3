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

#ifndef CGOGN_MODULE_SURFACE_FILTERING_H_
#define CGOGN_MODULE_SURFACE_FILTERING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/filtering.h>
#include <cgogn/geometry/algos/laplacian.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceFiltering : public Module
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceFiltering can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceFiltering(const App& app)
		: Module(app, "SurfaceFiltering (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_vertex_attribute_(nullptr)
	{
	}
	~SurfaceFiltering()
	{
	}

	void filter_average(MESH& m, Attribute<Vec3>* vertex_attribute)
	{
		std::shared_ptr<Attribute<Vec3>> filtered_vertex_attribute =
			add_attribute<Vec3, Vertex>(m, "__filtered_attribute");
		geometry::filter_average<Vec3>(m, vertex_attribute, filtered_vertex_attribute.get());
		vertex_attribute->swap(filtered_vertex_attribute.get());
		remove_attribute<Vertex>(m, filtered_vertex_attribute);

		mesh_provider_->emit_attribute_changed(m, vertex_attribute);
	}

	void filter_bilateral(MESH& m, Attribute<Vec3>* vertex_position)
	{
		std::shared_ptr<Attribute<Vec3>> filtered_vertex_position =
			add_attribute<Vec3, Vertex>(m, "__filtered_position");
		geometry::filter_bilateral(m, vertex_position, filtered_vertex_position.get());
		vertex_position->swap(filtered_vertex_position.get());
		remove_attribute<Vertex>(m, filtered_vertex_position);

		mesh_provider_->emit_attribute_changed(m, vertex_position);
	}

	void regularize(MESH& m, const Attribute<Vec3>* vertex_position, Attribute<Vec3>* vertex_attribute,
					Scalar fit_to_data)
	{
		geometry::filter_regularize(m, vertex_position, vertex_attribute, fit_to_data);
		mesh_provider_->emit_attribute_changed(m, vertex_attribute);
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
			selected_vertex_attribute_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_attribute_, "Attribute to filter",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_attribute_ = attribute; });

			if (selected_vertex_attribute_)
			{
				if (ImGui::Button("Filter average"))
					filter_average(*selected_mesh_, selected_vertex_attribute_.get());
			}

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			if (selected_vertex_position_)
			{
				if (selected_vertex_attribute_)
				{
					static float fit_to_data = 10.0;
					ImGui::SliderFloat("Fit to data", &fit_to_data, 0.0, 100.0);
					if (ImGui::Button("Regularize"))
						regularize(*selected_mesh_, selected_vertex_position_.get(), selected_vertex_attribute_.get(),
								   fit_to_data);
				}

				if (ImGui::Button("Filter bilateral"))
					filter_bilateral(*selected_mesh_, selected_vertex_position_.get());
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_attribute_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_FILTERING_H_
