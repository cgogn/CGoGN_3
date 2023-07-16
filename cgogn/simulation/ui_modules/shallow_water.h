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

#ifndef CGOGN_MODULE_SHALLOW_WATER_H_
#define CGOGN_MODULE_SHALLOW_WATER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/simulation/algos/shallow_water/shallow_water.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace ui
{

using geometry::Vec3;
using geometry::Scalar;

template <typename MESH>
class ShallowWater : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "ShallowWater can only be used with meshes of dimension 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using BoundaryCondition = simulation::shallow_water::BoundaryCondition;

public:
	ShallowWater(const App& app) : Module(app, "ShallowWater (" + std::string{mesh_traits<MESH>::name} + ")")
	{
	}
	~ShallowWater()
	{
	}

	void set_domain(MESH* m)
	{
		domain_ = m;
		domain_initialized_ = false;
		init_domain();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() { update_render_data(); });
	}

	void init_domain()
	{
		if (!domain_)
			domain_initialized_ = false;
		else
		{
			simulation::shallow_water::get_attributes(*domain_, sw_attributes_);
			simulation::shallow_water::init_attributes(*domain_, sw_attributes_, sw_context_);

			vertex_water_position_ = get_attribute<Vec3, Vertex>(*domain_, "water_position");
			if (!vertex_water_position_)
				vertex_water_position_ = add_attribute<Vec3, Vertex>(*domain_, "water_position");

			vertex_water_flux_ = get_attribute<Vec3, Vertex>(*domain_, "water_flux");
			if (!vertex_water_flux_)
				vertex_water_flux_ = add_attribute<Vec3, Vertex>(*domain_, "water_flux");

			vertex_water_position_->copy(sw_attributes_.vertex_position_.get());

			domain_connections_.clear();
			domain_connections_.push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					domain_, [this](Attribute<Vec3>* attribute) {
						if (sw_attributes_.vertex_position_.get() == attribute)
						{
							simulation::shallow_water::domain_geometry_changed(*domain_, sw_attributes_, sw_context_);
							update_render_data(true);
						}
					}));

			running_ = false;
			domain_initialized_ = true;

			update_render_data();
		}
	}

	void start()
	{
		cgogn_message_assert(domain_initialized_, "Domain is not initialized");

		running_ = true;

		launch_thread([this]() {
			while (this->running_)
			{
				simulation::shallow_water::execute_time_step(*domain_, sw_attributes_, sw_context_);
				// if (sw_context_.t_ == sw_context_.t_max_)
				// 	stop();
			}
		});

		app_.start_timer(50, [this]() -> bool { return !running_; });
	}

	void stop()
	{
		cgogn_message_assert(domain_initialized_, "Domain is not initialized");

		running_ = false;
	}

	void update_render_data(bool update_position = false)
	{
		cgogn_message_assert(domain_initialized_, "Domain is not initialized");

		parallel_foreach_cell(*domain_, [&](Vertex v) -> bool {
			Scalar h = 0.0;
			Scalar q = 0.0;
			Scalar r = 0.0;
			uint32 nbf = 0;
			foreach_incident_face(*domain_, v, [&](Face f) -> bool {
				uint32 fidx = index_of(*domain_, f);
				h += (*sw_attributes_.face_h_)[fidx];
				q += (*sw_attributes_.face_q_)[fidx];
				r += (*sw_attributes_.face_r_)[fidx];
				++nbf;
				return true;
			});
			h /= nbf;
			uint32 vidx = index_of(*domain_, v);
			if (update_position)
			{
				const Vec3& p = (*sw_attributes_.vertex_position_)[vidx];
				(*vertex_water_position_)[vidx] = {p[0], p[1], h};
			}
			else
				(*vertex_water_position_)[vidx][2] = h;
			(*vertex_water_flux_)[vidx] = {(q / nbf) / h, (r / nbf) / h, 0.0};
			return true;
		});

		mesh_provider_->emit_attribute_changed(*domain_, vertex_water_position_.get());
		mesh_provider_->emit_attribute_changed(*domain_, vertex_water_flux_.get());
		mesh_provider_->emit_attribute_changed(*domain_, sw_attributes_.face_h_.get());
	}

	void left_panel() override
	{
		if (domain_initialized_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
			MeshData<MESH>& md = mesh_provider_->mesh_data(*domain_);

			if (!running_)
			{
				if (ImGui::Button("Start"))
					start();
				ImGui::SameLine();
				if (ImGui::Button("step"))
				{
					simulation::shallow_water::execute_time_step(*domain_, sw_attributes_, sw_context_);
					update_render_data();
				}
			}
			else
			{
				if (ImGui::Button("Stop"))
					stop();
			}
			ImGui::Text("Simulation time: %f", sw_context_.t_);
			ImGui::Text("Current time step: %f", sw_context_.dt_);

			ImGui::Separator();

			if (ImGui::BeginCombo("Faces", selected_faces_set_ ? selected_faces_set_->name().c_str() : "-- select --"))
			{
				md.template foreach_cells_set<Face>([&](CellsSet<MESH, Face>& cs) {
					bool is_selected = &cs == selected_faces_set_;
					if (ImGui::Selectable(cs.name().c_str(), is_selected))
						selected_faces_set_ = &cs;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_faces_set_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##selected_faces_set"))
					selected_faces_set_ = nullptr;

				ImGui::Separator();

				for (Face f : *selected_faces_set_)
					ImGui::InputDouble(("h" + std::to_string(index_of(*domain_, f))).c_str(),
									   &value<Scalar>(*domain_, sw_attributes_.face_h_, f), 0.01f, 1.0f, "%.3f");
			}

			ImGui::Separator();

			if (ImGui::BeginCombo("Edges", selected_edges_set_ ? selected_edges_set_->name().c_str() : "-- select --"))
			{
				md.template foreach_cells_set<Edge>([&](CellsSet<MESH, Edge>& cs) {
					bool is_selected = &cs == selected_edges_set_;
					if (ImGui::Selectable(cs.name().c_str(), is_selected))
						selected_edges_set_ = &cs;
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (selected_edges_set_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##selected_edges_set"))
					selected_edges_set_ = nullptr;

				ImGui::Separator();

				for (Edge e : *selected_edges_set_)
				{
					BoundaryCondition e_bc_type = value<BoundaryCondition>(*domain_, sw_attributes_.edge_bc_type_, e);
					if (ImGui::BeginCombo(("BC_Type_" + std::to_string(index_of(*domain_, e))).c_str(),
										  simulation::shallow_water::bc_name(e_bc_type).c_str()))
					{
						for (uint32 i = 0; i < 6; ++i)
						{
							BoundaryCondition bc = static_cast<BoundaryCondition>(i);
							bool is_selected = bc == e_bc_type;
							if (ImGui::Selectable(simulation::shallow_water::bc_name(bc).c_str(), is_selected))
								value<BoundaryCondition>(*domain_, sw_attributes_.edge_bc_type_, e) = bc;
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						}
						ImGui::EndCombo();
					}
					ImGui::InputDouble(("BC_Value_" + std::to_string(index_of(*domain_, e))).c_str(),
									   &value<Scalar>(*domain_, sw_attributes_.edge_bc_value_, e), 0.01f, 1.0f, "%.3f");
				}
			}
		}
	}

private:
	MeshProvider<MESH>* mesh_provider_ = nullptr;
	MESH* domain_ = nullptr;

	bool domain_initialized_ = false;
	bool running_ = false;

	std::shared_ptr<Attribute<Vec3>> vertex_water_position_;
	std::shared_ptr<Attribute<Vec3>> vertex_water_flux_;

	simulation::shallow_water::Attributes<MESH> sw_attributes_;
	simulation::shallow_water::Context sw_context_;

	CellsSet<MESH, Face>* selected_faces_set_ = nullptr;
	CellsSet<MESH, Edge>* selected_edges_set_ = nullptr;

	std::vector<std::shared_ptr<boost::synapse::connection>> domain_connections_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHALLOW_WATER_H_
