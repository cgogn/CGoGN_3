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

#ifndef CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_H_
#define CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/curvature.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceDifferentialProperties : public Module
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Edge = typename mesh_traits<MESH>::Edge;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	SurfaceDifferentialProperties(const App& app) :
		Module(app, "SurfaceDifferentialProperties (" + mesh_traits<MESH>::name + ")"),
		selected_mesh_(nullptr),
		selected_vertex_position_(nullptr),
		selected_vertex_normal_(nullptr),
		selected_vertex_kmax_(nullptr),
		selected_vertex_kmin_(nullptr),
		selected_vertex_Kmax_(nullptr),
		selected_vertex_Kmin_(nullptr),
		selected_vertex_Knormal_(nullptr)
	{}
	~SurfaceDifferentialProperties()
	{}
    
	void init()
	{
		mesh_provider_ = static_cast<MeshProvider<MESH>*>(app_.module("MeshProvider (" + mesh_traits<MESH>::name + ")"));
	}

	void compute_normal(
		const MESH& m,
		const Attribute<Vec3>* vertex_position,
		Attribute<Vec3>* vertex_normal
	)
	{
		geometry::compute_normal(m, vertex_position, vertex_normal);
		mesh_provider_->emit_attribute_changed(&m, vertex_normal);
	}

	void compute_curvature(
		const MESH& m,
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
			m, radius,
			vertex_position, vertex_normal,	edge_angle,
			vertex_kmax, vertex_kmin,
			vertex_Kmax, vertex_Kmin, vertex_Knormal
		);
		mesh_provider_->emit_attribute_changed(&m, vertex_kmax);
		mesh_provider_->emit_attribute_changed(&m, vertex_kmin);
		mesh_provider_->emit_attribute_changed(&m, vertex_Kmax);
		mesh_provider_->emit_attribute_changed(&m, vertex_Kmin);
		mesh_provider_->emit_attribute_changed(&m, vertex_Knormal);
	}

protected:

    void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Select mesh"))
		{
			mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
			{
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
				{
					selected_mesh_ = m;
					selected_vertex_position_.reset();
					selected_vertex_normal_.reset();
					selected_vertex_kmax_.reset();
					selected_vertex_kmin_.reset();
					selected_vertex_Kmax_.reset();
					selected_vertex_Kmin_.reset();
					selected_vertex_Knormal_.reset();
				}
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			std::string selected_vertex_position_name_ = selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
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
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					selected_vertex_position_.reset();
			}

			std::string selected_vertex_normal_name_ = selected_vertex_normal_ ? selected_vertex_normal_->name() : "-- select --";
			if (ImGui::BeginCombo("Normal", selected_vertex_normal_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_normal_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_normal_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_normal_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##normal"))
					selected_vertex_normal_.reset();
			}

			std::string selected_vertex_kmax_name_ = selected_vertex_kmax_ ? selected_vertex_kmax_->name() : "-- select --";
			if (ImGui::BeginCombo("kmax", selected_vertex_kmax_name_.c_str()))
			{
				foreach_attribute<Scalar, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Scalar>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_kmax_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_kmax_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_kmax_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##kmax"))
					selected_vertex_kmax_.reset();
			}

			std::string selected_vertex_kmin_name_ = selected_vertex_kmin_ ? selected_vertex_kmin_->name() : "-- select --";
			if (ImGui::BeginCombo("kmin", selected_vertex_kmin_name_.c_str()))
			{
				foreach_attribute<Scalar, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Scalar>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_kmin_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_kmin_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_kmin_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##kmin"))
					selected_vertex_kmin_.reset();
			}

			std::string selected_vertex_Kmax_name_ = selected_vertex_Kmax_ ? selected_vertex_Kmax_->name() : "-- select --";
			if (ImGui::BeginCombo("Kmax", selected_vertex_Kmax_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_Kmax_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_Kmax_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_Kmax_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##Kmax"))
					selected_vertex_Kmax_.reset();
			}

			std::string selected_vertex_Kmin_name_ = selected_vertex_Kmin_ ? selected_vertex_Kmin_->name() : "-- select --";
			if (ImGui::BeginCombo("Kmin", selected_vertex_Kmin_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_Kmin_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_Kmin_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_Kmin_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##Kmin"))
					selected_vertex_Kmin_.reset();
			}

			std::string selected_vertex_Knormal_name_ = selected_vertex_Knormal_ ? selected_vertex_Knormal_->name() : "-- select --";
			if (ImGui::BeginCombo("Knormal", selected_vertex_Knormal_name_.c_str()))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == selected_vertex_Knormal_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						selected_vertex_Knormal_ = attribute;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_vertex_Knormal_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##Knormal"))
					selected_vertex_Knormal_.reset();
			}

			if (selected_vertex_position_)
			{
				if (ImGui::Button("Compute normal"))
				{
					if (!selected_vertex_normal_)
						selected_vertex_normal_ = add_attribute<Vec3, Vertex>(*selected_mesh_, "normal");
					compute_normal(*selected_mesh_, selected_vertex_position_.get(), selected_vertex_normal_.get());
				}
			}

			if (selected_vertex_position_ && selected_vertex_normal_)
			{
				if (ImGui::Button("Compute curvature"))
				{
					if (!selected_vertex_kmax_)
						selected_vertex_kmax_ = add_attribute<Scalar, Vertex>(*selected_mesh_, "kmax");
					if (!selected_vertex_kmin_)
						selected_vertex_kmin_ = add_attribute<Scalar, Vertex>(*selected_mesh_, "kmin");
					if (!selected_vertex_Kmax_)
						selected_vertex_Kmax_ = add_attribute<Vec3, Vertex>(*selected_mesh_, "Kmax");
					if (!selected_vertex_Kmin_)
						selected_vertex_Kmin_ = add_attribute<Vec3, Vertex>(*selected_mesh_, "Kmin");
					if (!selected_vertex_Knormal_)
						selected_vertex_Knormal_ = add_attribute<Vec3, Vertex>(*selected_mesh_, "Knormal");

					std::shared_ptr<Attribute<Scalar>> edge_angle = add_attribute<Scalar, Edge>(*selected_mesh_, "__edge_angle");
					geometry::compute_angle(*selected_mesh_, selected_vertex_position_.get(), edge_angle.get());

					Scalar mean_edge_length = geometry::mean_edge_length(*selected_mesh_, selected_vertex_position_.get());

					compute_curvature(
						*selected_mesh_,
						mean_edge_length * 4.0,
						selected_vertex_position_.get(),
						selected_vertex_normal_.get(),
						edge_angle.get(),
						selected_vertex_kmax_.get(),
						selected_vertex_kmin_.get(),
						selected_vertex_Kmax_.get(),
						selected_vertex_Kmin_.get(),
						selected_vertex_Knormal_.get()
					);

					remove_attribute<Edge>(*selected_mesh_, edge_angle);
				}
			}
		}

		ImGui::End();
	}

private:

	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_normal_;
	std::shared_ptr<Attribute<Scalar>> selected_vertex_kmax_;
	std::shared_ptr<Attribute<Scalar>> selected_vertex_kmin_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_Kmax_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_Kmin_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_Knormal_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_DIFFERENTIAL_PROPERTIES_H_
