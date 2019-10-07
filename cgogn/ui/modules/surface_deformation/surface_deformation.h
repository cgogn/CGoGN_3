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

#ifndef CGOGN_MODULE_SURFACE_DEFORMATION_H_
#define CGOGN_MODULE_SURFACE_DEFORMATION_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <GLFW/glfw3.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceDeformation : public ViewModule
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters() :
			vertex_position_(nullptr),
			selected_free_vertices_set_(nullptr),
			selected_handle_vertices_set_(nullptr)
		{}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		CellsSet<MESH, Vertex>* selected_free_vertices_set_;
		CellsSet<MESH, Vertex>* selected_handle_vertices_set_;
	};

public:

	SurfaceDeformation(const App& app) :
		ViewModule(app, "SurfaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"),
		selected_mesh_(nullptr),
		dragging_(false)
	{}
	~SurfaceDeformation()
	{}

private:

	void init_mesh(MESH* m)
	{
		parameters_[m];
	}

public:

	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;
	}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
				mesh_provider_, this, &SurfaceDeformation<MESH>::init_mesh
			)
		);
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				if (p.vertex_position_ && p.selected_handle_vertices_set_)
					dragging_ = true;
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			dragging_ = false;
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (dragging_)
		{
			Parameters& p = parameters_[selected_mesh_];
			
			drag_z_ = 0.0;
			p.selected_handle_vertices_set_->foreach_cell([&] (Vertex v)
			{
				drag_z_ += value<Vec3>(*selected_mesh_, p.vertex_position_, v)[2];
			});
			drag_z_ /= p.selected_handle_vertices_set_->size();
			
			rendering::GLVec3d start = view->unproject(view->previous_mouse_x(), view->previous_mouse_y(), drag_z_);
			rendering::GLVec3d end = view->unproject(x, y, drag_z_);
			
			Vec3 t = end - start;

			p.selected_handle_vertices_set_->foreach_cell([&] (Vertex v)
			{
				value<Vec3>(*selected_mesh_, p.vertex_position_, v) += t;
			});

			mesh_provider_->emit_attribute_changed(selected_mesh_, p.vertex_position_.get());
		}
	}

    void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Mesh"))
		{
			mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
			{
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
					selected_mesh_ = m;
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);
			Parameters& p = parameters_[selected_mesh_];

			if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == p.vertex_position_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						set_vertex_position(*selected_mesh_, attribute);
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (p.vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					set_vertex_position(*selected_mesh_, nullptr);
				
				ImGui::Separator();

				if (ImGui::BeginCombo("Free vertices", p.selected_free_vertices_set_ ? p.selected_free_vertices_set_->name().c_str() : "-- select --"))
				{
					md->template foreach_cells_set<Vertex>([&] (CellsSet<MESH, Vertex>& cs)
					{
						bool is_selected = &cs == p.selected_free_vertices_set_;
						if (ImGui::Selectable(cs.name().c_str(), is_selected))
							p.selected_free_vertices_set_ = &cs;
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
					ImGui::EndCombo();
				}
				if (p.selected_free_vertices_set_)
				{
					ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
					if (ImGui::Button("X##selected_vertices_set"))
						p.selected_free_vertices_set_ = nullptr;
				}

				if (ImGui::BeginCombo("Handle vertices", p.selected_handle_vertices_set_ ? p.selected_handle_vertices_set_->name().c_str() : "-- select --"))
				{
					md->template foreach_cells_set<Vertex>([&] (CellsSet<MESH, Vertex>& cs)
					{
						bool is_selected = &cs == p.selected_handle_vertices_set_;
						if (ImGui::Selectable(cs.name().c_str(), is_selected))
							p.selected_handle_vertices_set_ = &cs;
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
					ImGui::EndCombo();
				}
				if (p.selected_handle_vertices_set_)
				{
					ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
					if (ImGui::Button("X##selected_vertices_set"))
						p.selected_handle_vertices_set_ = nullptr;
				}
			}
		}
		
		ImGui::End();
	}

private:

	MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	MeshProvider<MESH>* mesh_provider_;

	bool dragging_;
	float64 drag_z_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_DEFORMATION_H_
