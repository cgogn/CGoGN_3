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

#include <cgogn/ui/modules/surface_differential_properties/surface_differential_properties.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/modules/cmap_provider/cmap_provider.h>

#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/algos/normal.h>

#include <imgui/imgui.h>

namespace cgogn
{

namespace ui
{

SurfaceDifferentialProperties::SurfaceDifferentialProperties(const App& app) :
	Module(app, "SurfaceDifferentialProperties"),
	selected_mesh_(nullptr)
{}

SurfaceDifferentialProperties::~SurfaceDifferentialProperties()
{}

void SurfaceDifferentialProperties::init()
{
	cmap_provider_ = static_cast<CMapProvider*>(app_.module("CMapProvider"));
}

void SurfaceDifferentialProperties::compute_normal(CMap2& m, const AttributePtr<Vec3>& vertex_position, AttributePtr<Vec3>& vertex_normal)
{
	if (!vertex_normal)
		vertex_normal = cgogn::add_attribute<Vec3, Vertex>(m, "normal");
	cgogn::geometry::compute_normal(m, vertex_position, vertex_normal);
}

void SurfaceDifferentialProperties::interface()
{
	ImGui::Begin("SurfaceDifferentialProperties", nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0, 0});

	if (ImGui::ListBoxHeader("Select mesh"))
	{
		cmap_provider_->foreach_cmap2([this] (Mesh* m, const std::string& name)
		{
			if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
			{
				selected_mesh_ = m;
				selected_vertex_position_.reset();
				selected_vertex_normal_.reset();
			}
		});
		ImGui::ListBoxFooter();
	}

	if (selected_mesh_)
	{
		std::string selected_vertex_position_name_ = selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
		if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
		{
			cgogn::foreach_attribute<Vec3, Vertex>(*selected_mesh_, [this] (const AttributePtr<Vec3>& attribute)
			{
				bool is_selected = attribute == selected_vertex_position_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_position_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			});
			ImGui::EndCombo();
		}
		std::string selected_vertex_normal_name_ = selected_vertex_normal_ ? selected_vertex_normal_->name() : "-- select --";
		if (ImGui::BeginCombo("Normal", selected_vertex_normal_name_.c_str()))
		{
			cgogn::foreach_attribute<Vec3, Vertex>(*selected_mesh_, [this] (const AttributePtr<Vec3>& attribute)
			{
				bool is_selected = attribute == selected_vertex_normal_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_normal_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			});
			ImGui::EndCombo();
		}

		if (selected_vertex_position_)
		{
			if (ImGui::Button("Compute normal"))
				compute_normal(*selected_mesh_, selected_vertex_position_, selected_vertex_normal_);
		}
	}

	ImGui::End();
}

} // namespace ui

} // namespace cgogn
