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
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/curvature.h>

#include <imgui/imgui.h>

namespace cgogn
{

namespace ui
{

SurfaceDifferentialProperties::SurfaceDifferentialProperties(const App& app) :
	Module(app, "SurfaceDifferentialProperties"),
	selected_mesh_(nullptr),
	selected_vertex_position_(nullptr),
	selected_vertex_normal_(nullptr),
	selected_vertex_kmax_(nullptr),
	selected_vertex_kmin_(nullptr),
	selected_vertex_Kmax_(nullptr),
	selected_vertex_Kmin_(nullptr),
	selected_vertex_Knormal_(nullptr)
{}

SurfaceDifferentialProperties::~SurfaceDifferentialProperties()
{}

void SurfaceDifferentialProperties::init()
{
	mesh_provider_ = static_cast<MeshProvider*>(app_.module("MeshProvider"));
}

void SurfaceDifferentialProperties::compute_normal(
	const Mesh& m,
	const Attribute<Vec3>* vertex_position,
	Attribute<Vec3>* vertex_normal
)
{
	geometry::compute_normal(m, vertex_position, vertex_normal);
}

void SurfaceDifferentialProperties::compute_curvature(
	const Mesh& m,
	Scalar radius,
	const Attribute<Vec3>* vertex_position,
	const Attribute<Vec3>* vertex_normal,
	const Attribute<Scalar>* edge_angle,
	Attribute<Scalar>* vertex_kmax,
	Attribute<Scalar>* vertex_kmin,
	Attribute<Vec3>* vertex_Kmax,
	Attribute<Vec3>* vertex_Kmin,
	Attribute<Vec3>* vertex_Knormal
)
{
	geometry::compute_curvature(
		m,
		radius,
		vertex_position,
		vertex_normal,
		edge_angle,
		vertex_kmax,
		vertex_kmin,
		vertex_Kmax,
		vertex_Kmin,
		vertex_Knormal
	);
}

void SurfaceDifferentialProperties::interface()
{
	ImGui::Begin("SurfaceDifferentialProperties", nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0, 0});

	if (ImGui::ListBoxHeader("Select mesh"))
	{
		mesh_provider_->foreach_mesh([this] (Mesh* m, const std::string& name)
		{
			if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
			{
				selected_mesh_ = m;
				selected_vertex_position_ = nullptr;
				selected_vertex_normal_ = nullptr;
			}
		});
		ImGui::ListBoxFooter();
	}

	if (selected_mesh_)
	{
		std::vector<Attribute<Vec3>*> vec3_attributes;
		foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (Attribute<Vec3>* attribute)
		{
			vec3_attributes.push_back(attribute);
		});
		std::vector<Attribute<Scalar>*> scalar_attributes;
		foreach_attribute<Scalar, Vertex>(*selected_mesh_, [&] (Attribute<Scalar>* attribute)
		{
			scalar_attributes.push_back(attribute);
		});

		std::string selected_vertex_position_name_ = selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
		if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
		{
			for (Attribute<Vec3>* attribute : vec3_attributes)
			{
				bool is_selected = attribute == selected_vertex_position_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_position_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}
		std::string selected_vertex_normal_name_ = selected_vertex_normal_ ? selected_vertex_normal_->name() : "-- select --";
		if (ImGui::BeginCombo("Normal", selected_vertex_normal_name_.c_str()))
		{
			for (Attribute<Vec3>* attribute : vec3_attributes)
			{
				bool is_selected = attribute == selected_vertex_normal_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_normal_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}
		std::string selected_vertex_kmax_name_ = selected_vertex_kmax_ ? selected_vertex_kmax_->name() : "-- select --";
		if (ImGui::BeginCombo("kmax", selected_vertex_kmax_name_.c_str()))
		{
			for (Attribute<Scalar>* attribute : scalar_attributes)
			{
				bool is_selected = attribute == selected_vertex_kmax_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_kmax_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}
		std::string selected_vertex_kmin_name_ = selected_vertex_kmin_ ? selected_vertex_kmin_->name() : "-- select --";
		if (ImGui::BeginCombo("kmin", selected_vertex_kmin_name_.c_str()))
		{
			for (Attribute<Scalar>* attribute : scalar_attributes)
			{
				bool is_selected = attribute == selected_vertex_kmin_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_kmin_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}
		std::string selected_vertex_Kmax_name_ = selected_vertex_Kmax_ ? selected_vertex_Kmax_->name() : "-- select --";
		if (ImGui::BeginCombo("Kmax", selected_vertex_Kmax_name_.c_str()))
		{
			for (Attribute<Vec3>* attribute : vec3_attributes)
			{
				bool is_selected = attribute == selected_vertex_Kmax_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_Kmax_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}
		std::string selected_vertex_Kmin_name_ = selected_vertex_Kmin_ ? selected_vertex_Kmin_->name() : "-- select --";
		if (ImGui::BeginCombo("Kmin", selected_vertex_Kmin_name_.c_str()))
		{
			for (Attribute<Vec3>* attribute : vec3_attributes)
			{
				bool is_selected = attribute == selected_vertex_Kmin_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_Kmin_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}
		std::string selected_vertex_Knormal_name_ = selected_vertex_Knormal_ ? selected_vertex_Knormal_->name() : "-- select --";
		if (ImGui::BeginCombo("Knormal", selected_vertex_Knormal_name_.c_str()))
		{
			for (Attribute<Vec3>* attribute : vec3_attributes)
			{
				bool is_selected = attribute == selected_vertex_Knormal_;
				if (ImGui::Selectable(attribute->name().c_str(), is_selected))
					selected_vertex_Knormal_ = attribute;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}

		if (selected_vertex_position_)
		{
			if (ImGui::Button("Compute normal"))
			{
				if (!selected_vertex_normal_)
				{
					std::shared_ptr<Attribute<Vec3>> vertex_normal = add_attribute<Vec3, Vertex>(*selected_mesh_, "normal");
					selected_vertex_normal_ = vertex_normal.get();
				}
				compute_normal(*selected_mesh_, selected_vertex_position_, selected_vertex_normal_);
			}
		}

		if (selected_vertex_position_ && selected_vertex_normal_)
		{
			if (ImGui::Button("Compute curvature"))
			{
				if (!selected_vertex_kmax_)
				{
					std::shared_ptr<Attribute<Scalar>> vertex_kmax = add_attribute<Scalar, Vertex>(*selected_mesh_, "kmax");
					selected_vertex_kmax_ = vertex_kmax.get();
				}
				if (!selected_vertex_kmin_)
				{
					std::shared_ptr<Attribute<Scalar>> vertex_kmin = add_attribute<Scalar, Vertex>(*selected_mesh_, "kmin");
					selected_vertex_kmin_ = vertex_kmin.get();
				}
				if (!selected_vertex_Kmax_)
				{
					std::shared_ptr<Attribute<Vec3>> vertex_Kmax = add_attribute<Vec3, Vertex>(*selected_mesh_, "Kmax");
					selected_vertex_Kmax_ = vertex_Kmax.get();
				}
				if (!selected_vertex_Kmin_)
				{
					std::shared_ptr<Attribute<Vec3>> vertex_Kmin = add_attribute<Vec3, Vertex>(*selected_mesh_, "Kmin");
					selected_vertex_Kmin_ = vertex_Kmin.get();
				}
				if (!selected_vertex_Knormal_)
				{
					std::shared_ptr<Attribute<Vec3>> vertex_Knormal = add_attribute<Vec3, Vertex>(*selected_mesh_, "Knormal");
					selected_vertex_Knormal_ = vertex_Knormal.get();
				}

				std::shared_ptr<Attribute<Scalar>> edge_angle = add_attribute<Scalar, Edge>(*selected_mesh_, "__edge_angle");
				geometry::compute_angle(*selected_mesh_, selected_vertex_position_, edge_angle.get());

				Scalar mean_edge_length = geometry::mean_edge_length(*selected_mesh_, selected_vertex_position_);

				compute_curvature(
					*selected_mesh_,
					mean_edge_length * 4.0,
					selected_vertex_position_,
					selected_vertex_normal_,
					edge_angle.get(),
					selected_vertex_kmax_,
					selected_vertex_kmin_,
					selected_vertex_Kmax_,
					selected_vertex_Kmin_,
					selected_vertex_Knormal_
				);

				remove_attribute<Edge>(*selected_mesh_, edge_angle);
			}
		}
	}

	ImGui::End();
}

} // namespace ui

} // namespace cgogn
