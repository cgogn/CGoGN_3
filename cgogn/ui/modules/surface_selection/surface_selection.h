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

#ifndef CGOGN_MODULE_SURFACE_SELECTION_H_
#define CGOGN_MODULE_SURFACE_SELECTION_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/picking.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/rendering/vbo_update.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceSelection : public ViewModule
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Face = typename mesh_traits<MESH>::Face;

	enum SelectingCell
	{
		VertexSelect = 0,
		EdgeSelect,
		FaceSelect
	};

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters() :
			vertex_position_(nullptr),
			vertex_scale_factor_(1.0),
			selected_vertices_set_(nullptr),
			selected_edges_set_(nullptr),
			selected_faces_set_(nullptr),
			selecting_cell_(VertexSelect)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65);
			param_point_sprite_->set_vbos(&selected_vertices_vbo_);

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 0, 0, 0.65);
			param_edge_->width_= 2.0f;
			param_edge_->set_vbos(&selected_edges_vbo_);

			param_flat_ =  rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(1, 0, 0, 0.65);
			param_flat_->back_color_ = rendering::GLColor(1, 0, 0, 0.65);
			param_flat_->double_side_ = true;
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_flat_->set_vbos(&selected_faces_vbo_);
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

public:

		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell([&] (Vertex v)
				{
					selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v));
				});
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		void update_selected_edges_vbo()
		{
			if (selected_edges_set_)
			{
				std::vector<Vec3> selected_edges_position;
				selected_edges_position.reserve(selected_edges_set_->size() * 2);
				selected_edges_set_->foreach_cell([&] (Edge e)
				{
					std::vector<Vertex> vertices = incident_vertices(*mesh_, e);
					selected_edges_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[0]));
					selected_edges_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[1]));
				});
				rendering::update_vbo(selected_edges_position, &selected_edges_vbo_);
			}
		}

		void update_selected_faces_vbo()
		{
			if (selected_faces_set_)
			{
				std::vector<Vec3> selected_faces_position;
				selected_faces_position.reserve(selected_faces_set_->size() * 3); // TODO: manage polygonal faces
				selected_faces_set_->foreach_cell([&] (Face f)
				{
					foreach_incident_vertex(*mesh_, f, [&] (Vertex v) -> bool
					{
						selected_faces_position.push_back(value<Vec3>(*mesh_, vertex_position_, v));
						return true;
					});
				});
				rendering::update_vbo(selected_faces_position, &selected_faces_vbo_);
			}
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;

		rendering::VBO selected_vertices_vbo_;
		rendering::VBO selected_edges_vbo_;
		rendering::VBO selected_faces_vbo_;

		CellsSet<MESH, Vertex>* selected_vertices_set_;
		CellsSet<MESH, Edge>* selected_edges_set_;
		CellsSet<MESH, Face>* selected_faces_set_;

		SelectingCell selecting_cell_;
	};

public:

	SurfaceSelection(const App& app) :
		ViewModule(app, "SurfaceSelection (" + std::string{mesh_traits<MESH>::name} + ")"),
		selected_mesh_(nullptr)
	{}

	~SurfaceSelection()
	{}

private:

	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m] (Attribute<Vec3>* attribute)
				{
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6.0;
						p.update_selected_vertices_vbo();
						p.update_selected_edges_vbo();
						p.update_selected_faces_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}
			)
		);
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Vertex>>(
				m, [this, m] (CellsSet<MESH, Vertex>* set)
				{
					Parameters& p = parameters_[m];
					if (p.vertex_position_)
						p.update_selected_vertices_vbo();
					for (View* v : linked_views_)
						v->request_update();
				}
			)
		);
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Edge>>(
				m, [this, m] (CellsSet<MESH, Edge>* set)
				{
					Parameters& p = parameters_[m];
					if (p.vertex_position_)
						p.update_selected_edges_vbo();
					for (View* v : linked_views_)
						v->request_update();
				}
			)
		);
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Face>>(
				m, [this, m] (CellsSet<MESH, Face>* set)
				{
					Parameters& p = parameters_[m];
					if (p.vertex_position_)
						p.update_selected_faces_vbo();
					for (View* v : linked_views_)
						v->request_update();
				}
			)
		);
	}

public:

	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = geometry::mean_edge_length(m, p.vertex_position_.get()) / 6.0;
			p.update_selected_vertices_vbo();
			p.update_selected_edges_vbo();
			p.update_selected_faces_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
				mesh_provider_, this, &SurfaceSelection<MESH>::init_mesh
			)
		);
	}
	
	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selected_mesh_ && view->shift_pressed())
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);
			Parameters& p = parameters_[selected_mesh_];
			
			if (p.vertex_position_)
			{
				rendering::GLVec3d near = view->unproject(x, y, 0.0);
				rendering::GLVec3d far = view->unproject(x, y, 1.0);
				Vec3 A{ near.x(), near.y(), near.z() };
				Vec3 B{ far.x(), far.y(), far.z() };
				
				if (p.selecting_cell_ == VertexSelect && p.selected_vertices_set_)
				{	
					std::vector<Vertex> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
									
					if (!picked.empty())
					{
						if (button == 0)
							p.selected_vertices_set_->select(picked[0]);
						else if (button == 1)
							p.selected_vertices_set_->unselect(picked[0]);
						p.update_selected_vertices_vbo();
					}
				}
				else if (p.selecting_cell_ == EdgeSelect && p.selected_edges_set_)
				{
					std::vector<Edge> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
									
					if (!picked.empty())
					{
						if (button == 0)
							p.selected_edges_set_->select(picked[0]);
						else if (button == 1)
							p.selected_edges_set_->unselect(picked[0]);
						p.update_selected_edges_vbo();
					}
				}
				else if (p.selecting_cell_ == FaceSelect && p.selected_faces_set_)
				{
					std::vector<Face> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
									
					if (!picked.empty())
					{
						if (button == 0)
							p.selected_faces_set_->select(picked[0]);
						else if (button == 1)
							p.selected_faces_set_->unselect(picked[0]);
						p.update_selected_faces_vbo();
					}
				}
			}
		}
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.selecting_cell_ == VertexSelect &&
				p.selected_vertices_set_ && p.selected_vertices_set_->size() > 0 &&
				p.param_point_sprite_->vao_initialized())
			{
				p.param_point_sprite_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
				p.param_point_sprite_->release();
			}
			else if (p.selecting_cell_ == EdgeSelect &&
				p.selected_edges_set_ && p.selected_edges_set_->size() > 0 &&
				p.param_edge_->vao_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_LINES, 0, p.selected_edges_set_->size() * 2);
				p.param_edge_->release();
			}
			else if (p.selecting_cell_ == FaceSelect &&
				p.selected_faces_set_ && p.selected_faces_set_->size() > 0 &&
				p.param_flat_->vao_initialized())
			{
				p.param_flat_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_TRIANGLES, 0, p.selected_faces_set_->size() * 3); // TODO: manage polygonal faces
				p.param_flat_->release();
			}
		}
	}

    void interface() override
	{
		bool need_update = false;

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
				
				need_update |= ImGui::RadioButton("Vertex", (int*)(&p.selecting_cell_), VertexSelect); ImGui::SameLine();
        		need_update |= ImGui::RadioButton("Edge", (int*)(&p.selecting_cell_), EdgeSelect); ImGui::SameLine();
        		need_update |= ImGui::RadioButton("Face", (int*)(&p.selecting_cell_), FaceSelect);

				MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);
				Parameters& p = parameters_[selected_mesh_];

				if (p.selecting_cell_ == VertexSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_vertices_set_ ? p.selected_vertices_set_->name().c_str() : "-- select --"))
					{
						md->template foreach_cells_set<Vertex>([&] (CellsSet<MESH, Vertex>& cs)
						{
							bool is_selected = &cs == p.selected_vertices_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_vertices_set_ = &cs;
								p.update_selected_vertices_vbo();
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
					if (ImGui::Button("Create##vertices_set"))
						md->template add_cells_set<Vertex>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1, 2.0);
				}
				else if (p.selecting_cell_ == EdgeSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_edges_set_ ? p.selected_edges_set_->name().c_str() : "-- select --"))
					{
						md->template foreach_cells_set<Edge>([&] (CellsSet<MESH, Edge>& cs)
						{
							bool is_selected = &cs == p.selected_edges_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_edges_set_ = &cs;
								p.update_selected_edges_vbo();
							}
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						});
						ImGui::EndCombo();
					}
					if (p.selected_edges_set_)
					{
						ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
						if (ImGui::Button("X##selected_edges_set"))
							p.selected_edges_set_ = nullptr;
					}
					if (ImGui::Button("Create##edges_set"))
						md->template add_cells_set<Edge>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
				}
				else if (p.selecting_cell_ == FaceSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_faces_set_ ? p.selected_faces_set_->name().c_str() : "-- select --"))
					{
						md->template foreach_cells_set<Face>([&] (CellsSet<MESH, Face>& cs)
						{
							bool is_selected = &cs == p.selected_faces_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_faces_set_ = &cs;
								p.update_selected_faces_vbo();
							}
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						});
						ImGui::EndCombo();
					}
					if (p.selected_faces_set_)
					{
						ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
						if (ImGui::Button("X##selected_faces_set"))
							p.selected_faces_set_ = nullptr;
					}
					if (ImGui::Button("Create##faces_set"))
						md->template add_cells_set<Face>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
				}
			}
		}

		ImGui::End();

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

#endif // CGOGN_MODULE_SURFACE_SELECTION_H_
