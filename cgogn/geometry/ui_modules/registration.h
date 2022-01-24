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

#ifndef CGOGN_MODULE_REGISTRATION_H_
#define CGOGN_MODULE_REGISTRATION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/registration.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class Registration : public Module
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	Registration(const App& app) : Module(app, "Registration (" + std::string{mesh_traits<MESH>::name} + ")")
	{
	}
	~Registration()
	{
	}

	void rigid_register_mesh(MESH& source, Attribute<Vec3>* source_position, MESH& target,
							 const Attribute<Vec3>* target_position)
	{
		geometry::rigid_register_mesh(source, source_position, target, target_position);
		mesh_provider_->emit_attribute_changed(source, source_position);
	}

	void non_rigid_register_mesh(MESH& source, std::shared_ptr<Attribute<Vec3>>& source_position, MESH& target,
								 const Attribute<Vec3>* target_position, Scalar fit_to_target, bool relax)
	{
		geometry::non_rigid_register_mesh(source, source_position, target, target_position, fit_to_target, relax);
		mesh_provider_->emit_attribute_changed(source, source_position.get());
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void interface() override
	{
		imgui_mesh_selector(mesh_provider_, selected_source_mesh_, "Source mesh", [&](MESH& m) {
			selected_source_mesh_ = &m;
			selected_source_vertex_position_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
		if (selected_source_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(*selected_source_mesh_, selected_source_vertex_position_,
												"Source position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													selected_source_vertex_position_ = attribute;
												});
		}

		imgui_mesh_selector(mesh_provider_, selected_target_mesh_, "Target mesh", [&](MESH& m) {
			selected_target_mesh_ = &m;
			selected_target_vertex_position_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
		if (selected_target_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(*selected_target_mesh_, selected_target_vertex_position_,
												"Target position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													selected_target_vertex_position_ = attribute;
												});
		}

		if (selected_source_mesh_ && selected_source_vertex_position_ && selected_target_mesh_ &&
			selected_target_vertex_position_)
		{
			if (ImGui::Button("Rigid registration"))
				rigid_register_mesh(*selected_source_mesh_, selected_source_vertex_position_.get(),
									*selected_target_mesh_, selected_target_vertex_position_.get());
			static float fit_to_target = 0.05f;
			static bool relax = false;
			ImGui::SliderFloat("Fit to target", &fit_to_target, 0.0, 10.0);
			if (ImGui::Checkbox("Relax", &relax))
			{
				if (relax)
					fit_to_target = 5.0;
				else
					fit_to_target = 0.05;
			}
			if (ImGui::Button("Non-rigid registration"))
				non_rigid_register_mesh(*selected_source_mesh_, selected_source_vertex_position_,
										*selected_target_mesh_, selected_target_vertex_position_.get(), fit_to_target,
										relax);
		}
	}

private:
	MESH* selected_source_mesh_ = nullptr;
	std::shared_ptr<Attribute<Vec3>> selected_source_vertex_position_ = nullptr;
	MESH* selected_target_mesh_ = nullptr;
	std::shared_ptr<Attribute<Vec3>> selected_target_vertex_position_ = nullptr;
	MeshProvider<MESH>* mesh_provider_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_REGISTRATION_H_
