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

#ifndef CGOGN_MODULE_VOLUME_DEFORMATION_H_
#define CGOGN_MODULE_VOLUME_DEFORMATION_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <GLFW/glfw3.h>

#include <boost/synapse/connect.hpp>
#include <memory>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class VolumeDeformation : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 3, "VolumeDeformation can only be used with meshes of dimension 3");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters() : vertex_position_(nullptr), selected_handle_vertices_set_(nullptr)
		{
		}

		~Parameters()
		{
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		CellsSet<MESH, Vertex>* selected_handle_vertices_set_;
	};

public:
	VolumeDeformation(const App& app)
		: ViewModule(app, "VolumeDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  dragging_(false)
	{
	}
	~VolumeDeformation()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;
	}

	void set_selected_handle_vertices_set(const MESH& m, CellsSet<MESH, Vertex>* set)
	{
		Parameters& p = parameters_[&m];
		p.selected_handle_vertices_set_ = set;
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &VolumeDeformation<MESH>::init_mesh));
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				if (p.vertex_position_ && p.selected_handle_vertices_set_ &&
					p.selected_handle_vertices_set_->size() > 0)
				{
					drag_z_ = 0.0;
					p.selected_handle_vertices_set_->foreach_cell([&](Vertex v) {
						const Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_, v);
						rendering::GLVec4d vec(pos[0], pos[1], pos[2], 1.0);
						vec = view->projection_matrix_d() * view->modelview_matrix_d() * vec;
						vec /= vec[3];
						drag_z_ += (1.0 + vec[2]) / 2.0;
					});
					drag_z_ /= p.selected_handle_vertices_set_->size();
					previous_drag_pos_ = view->unproject(view->previous_mouse_x(), view->previous_mouse_y(), drag_z_);
					dragging_ = true;
				}
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		unused_parameters(view);
		if (key_code == GLFW_KEY_D)
		{
			if (dragging_)
				dragging_ = false;
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (dragging_)
		{
			Parameters& p = parameters_[selected_mesh_];

			rendering::GLVec3d drag_pos = view->unproject(x, y, drag_z_);
			Vec3 t = drag_pos - previous_drag_pos_;
			p.selected_handle_vertices_set_->foreach_cell(
				[&](Vertex v) { value<Vec3>(*selected_mesh_, p.vertex_position_, v) += t; });
			previous_drag_pos_ = drag_pos;

			mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
		}
	}

	void ui_interface() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
			Parameters& p = parameters_[selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_mesh_, attribute);
												});

			if (p.vertex_position_)
			{
				ImGui::Separator();

				if (ImGui::BeginCombo("Handle vertices", p.selected_handle_vertices_set_
															 ? p.selected_handle_vertices_set_->name().c_str()
															 : "-- select --"))
				{
					md.template foreach_cells_set<Vertex>([&](CellsSet<MESH, Vertex>& cs) {
						bool is_selected = &cs == p.selected_handle_vertices_set_;
						if (ImGui::Selectable(cs.name().c_str(), is_selected))
							set_selected_handle_vertices_set(*selected_mesh_, &cs);
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
					ImGui::EndCombo();
				}
				if (p.selected_handle_vertices_set_)
				{
					ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
					if (ImGui::Button("X##selected_handle_vertices_set"))
						set_selected_handle_vertices_set(*selected_mesh_, nullptr);
				}
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	MeshProvider<MESH>* mesh_provider_;

	bool dragging_;
	float64 drag_z_;
	rendering::GLVec3d previous_drag_pos_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_VOLUME_DEFORMATION_H_
