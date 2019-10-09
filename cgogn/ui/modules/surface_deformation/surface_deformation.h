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

#include <Eigen/Sparse>

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
    using Mat3 = geometry::Mat3;
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

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> working_LAPL_;
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> working_BILAPL_;

		bool initialized_;
		bool solver_ready_;

		std::shared_ptr<Attribute<Vec3>> vertex_position_init_;
		std::shared_ptr<Attribute<Vec3>> vertex_diff_coord_;
		std::shared_ptr<Attribute<Vec3>> vertex_bi_diff_coord_;
		std::shared_ptr<Attribute<Mat3>> vertex_rotation_matrix_;
		std::shared_ptr<Attribute<Vec3>> vertex_rotated_diff_coord_;
		std::shared_ptr<Attribute<Vec3>> vertex_rotated_bi_diff_coord_;
		std::shared_ptr<Attribute<uint32>> vertex_index_;

		std::shared_ptr<Attribute<Scalar>> edge_weight_;

		Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>>* solver_;
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
				if (p.vertex_position_ && p.selected_handle_vertices_set_ && p.selected_handle_vertices_set_->size() > 0)
				{
					drag_z_ = 0.0;
					p.selected_handle_vertices_set_->foreach_cell([&] (Vertex v)
					{
						const Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_, v);
						rendering::GLVec4d vec(pos[0], pos[1], pos[2], 1.0);
						vec = view->projection_matrix_d() * view->modelview_matrix_d() * vec;
						vec /= vec[3];
						drag_z_ += (1.0 + vec[2]) / 2.0;
					});
					drag_z_ /= p.selected_handle_vertices_set_->size();
					previous_drag_pos_ = view->unproject(view->previous_mouse_x(), view->previous_mouse_y(), drag_z_);

					for (View* v : linked_views_)
						v->lock_scene_bb();

					dragging_ = true;
				}
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			if (dragging_)
			{
				dragging_ = false;

				for (View* v : linked_views_)
					v->unlock_scene_bb();
			}
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (dragging_)
		{
			Parameters& p = parameters_[selected_mesh_];
			
			rendering::GLVec3d drag_pos = view->unproject(x, y, drag_z_);
			Vec3 t = drag_pos - previous_drag_pos_;
			p.selected_handle_vertices_set_->foreach_cell([&] (Vertex v)
			{
				value<Vec3>(*selected_mesh_, p.vertex_position_, v) += t;
			});
			previous_drag_pos_ = drag_pos;

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
					if (ImGui::Button("X##selected_free_vertices_set"))
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
					if (ImGui::Button("X##selected_handle_vertices_set"))
						p.selected_handle_vertices_set_ = nullptr;
					ImGui::TextUnformatted("Press D to drag the handle");
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
	rendering::GLVec3d previous_drag_pos_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_DEFORMATION_H_
