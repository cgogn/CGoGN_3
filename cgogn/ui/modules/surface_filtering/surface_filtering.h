/*******************************************************************************
* CGoGN                                                                        *
* Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/filtering.h>

namespace cgogn
{

namespace ui
{

class App;
template <typename MESH>
class MeshProvider;

template <typename MESH>
class SurfaceFiltering : public Module
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	SurfaceFiltering(const App& app) :
		Module(app, "SurfaceFiltering"),
		selected_mesh_(nullptr),
		selected_vertex_position_(nullptr)
	{}
	~SurfaceFiltering()
	{}

	void init()
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider"));
	}

	void filter_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		std::shared_ptr<Attribute<Vec3>> filtered_vertex_position = add_attribute<Vec3, Vertex>(m, "__filtered_position");
		for (auto it = filtered_vertex_position->begin(), end = filtered_vertex_position->end(); it != end; ++it)
			*it = (*vertex_position)[it.index()];
		
		geometry::filter_average<Vec3>(m, vertex_position, filtered_vertex_position.get());
		
		vertex_position->swap(filtered_vertex_position.get());
		remove_attribute<Vertex>(m, filtered_vertex_position);
	}

protected:

    void interface() override
	{
		ImGui::Begin("Filtering", nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Select mesh"))
		{
			mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
			{
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
				{
					selected_mesh_ = m;
					selected_vertex_position_ = nullptr;
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			std::string selected_vertex_position_name_ = selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [this] (Attribute<Vec3>* attribute)
				{
					bool is_selected = attribute == selected_vertex_position_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_position_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}

			if (selected_vertex_position_)
			{
				if (ImGui::Button("Filter"))
					filter_mesh(*selected_mesh_, selected_vertex_position_);
			}
		}
		
		ImGui::End();
	}

private:

	MESH* selected_mesh_;
	Attribute<Vec3>* selected_vertex_position_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_FILTERING_H_
