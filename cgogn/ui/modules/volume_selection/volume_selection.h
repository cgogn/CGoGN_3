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

#ifndef CGOGN_MODULE_VOLUME_SELECTION_H_
#define CGOGN_MODULE_VOLUME_SELECTION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class VolumeSelection : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 3, "VolumeSelection can only be used with meshes of dimension >= 3");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	enum SelectingCell
	{
		VertexSelect
	};

	enum SelectionMethod
	{
		SingleCell
	};

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_scale_factor_(1.0), selected_vertices_set_(nullptr),
			  selecting_cell_(VertexSelect), selection_method_(SingleCell), choosing_vertex_(false)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell(
					[&](Vertex v) { selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v)); });
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;

		rendering::VBO selected_vertices_vbo_;

		CellsSet<MESH, Vertex>* selected_vertices_set_;

		SelectingCell selecting_cell_;
		SelectionMethod selection_method_;

		std::vector<Vertex> picked_vertices_;
		uint32 current_candidate_index_;
		bool current_candidate_state_;
		bool choosing_vertex_;
		int32 clicked_button_;
	};

public:
	VolumeSelection(const App& app)
		: ViewModule(app, "VolumeSelection (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr)
	{
	}

	~VolumeSelection()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6);
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Vertex>>(
				m, [this, m](CellsSet<MESH, Vertex>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, p.vertex_position_.get()) / 6);
			p.update_selected_vertices_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &VolumeSelection<MESH>::init_mesh));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];

			if (p.choosing_vertex_)
				p.choosing_vertex_ = false;
			else if (view->shift_pressed())
			{
				if (p.vertex_position_)
				{
					rendering::GLVec3d near = view->unproject(x, y, 0.0);
					rendering::GLVec3d far = view->unproject(x, y, 1.0);
					Vec3 A{near.x(), near.y(), near.z()};
					Vec3 B{far.x(), far.y(), far.z()};

					switch (p.selection_method_)
					{
					case SingleCell: {
						switch (p.selecting_cell_)
						{
						case VertexSelect:
							if (p.selected_vertices_set_)
							{
								cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B,
														 p.picked_vertices_);
								if (!p.picked_vertices_.empty())
								{
									p.choosing_vertex_ = true;
									p.clicked_button_ = button;
									p.current_candidate_index_ = 0;
									p.current_candidate_state_ = p.selected_vertices_set_->contains(
										p.picked_vertices_[p.current_candidate_index_]);
									switch (button)
									{
									case 0:
										p.selected_vertices_set_->select(
											p.picked_vertices_[p.current_candidate_index_]);
										break;
									case 1:
										p.selected_vertices_set_->unselect(
											p.picked_vertices_[p.current_candidate_index_]);
										break;
									}
									mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_vertices_set_);
								}
							}
							break;
						}
						break;
					}
					}
				}
			}
		}
	}

	void mouse_wheel_event(View* view, int32, int32 dy) override
	{
		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];
			if (p.choosing_vertex_)
			{
				uint32 bak_index = p.current_candidate_index_;
				bool bak_state = p.current_candidate_state_;
				if (dy > 0)
					p.current_candidate_index_ = (p.current_candidate_index_ + 1) % p.picked_vertices_.size();
				else if (dy < 0)
					p.current_candidate_index_ =
						(p.current_candidate_index_ + p.picked_vertices_.size() - 1) % p.picked_vertices_.size();
				p.current_candidate_state_ =
					p.selected_vertices_set_->contains(p.picked_vertices_[p.current_candidate_index_]);
				switch (p.clicked_button_)
				{
				case 0:
					if (!bak_state)
						p.selected_vertices_set_->unselect(p.picked_vertices_[bak_index]);
					p.selected_vertices_set_->select(p.picked_vertices_[p.current_candidate_index_]);
					break;
				case 1:
					if (bak_state)
						p.selected_vertices_set_->select(p.picked_vertices_[bak_index]);
					p.selected_vertices_set_->unselect(p.picked_vertices_[p.current_candidate_index_]);
					break;
				}
				mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_vertices_set_);
				view->stop_event();
			}
		}
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.selecting_cell_ == VertexSelect && p.selected_vertices_set_ && p.selected_vertices_set_->size() > 0 &&
				p.param_point_sprite_->attributes_initialized())
			{
				if (p.choosing_vertex_)
					p.param_point_sprite_->color_ = rendering::GLColor(1, 1, 0, 0.65f);
				else
					p.param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
				p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
				p.param_point_sprite_->release();
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Volume", [&](MESH* m) {
			selected_mesh_ = m;
			mesh_provider_->mesh_data(selected_mesh_)->outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
			Parameters& p = parameters_[selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_mesh_, attribute);
												});

			if (p.vertex_position_)
			{
				ImGui::Separator();
				int* ptr_sel_cell = reinterpret_cast<int*>(&p.selecting_cell_);
				need_update |= ImGui::RadioButton("Vertex", ptr_sel_cell, VertexSelect);
				ImGui::SameLine();

				ImGui::RadioButton("Single", reinterpret_cast<int*>(&p.selection_method_), SingleCell);

				MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);

				if (p.selecting_cell_ == VertexSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_vertices_set_ ? p.selected_vertices_set_->name().c_str()
																		   : "-- select --"))
					{
						md->template foreach_cells_set<Vertex>([&](CellsSet<MESH, Vertex>& cs) {
							bool is_selected = &cs == p.selected_vertices_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_vertices_set_ = &cs;
								p.update_selected_vertices_vbo();
								for (View* v : linked_views_)
									v->request_update();
							}
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						});
						ImGui::EndCombo();
					}
					if (p.selected_vertices_set_)
					{
						ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
						if (ImGui::Button("X##selected_vertices_set"))
							p.selected_vertices_set_ = nullptr;
					}
					if (ImGui::Button("Create set##vertices_set"))
						md->template add_cells_set<Vertex>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	const MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_VOLUME_SELECTION_H_
