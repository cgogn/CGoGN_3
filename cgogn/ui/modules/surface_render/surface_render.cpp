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

#include <cgogn/ui/modules/surface_render/surface_render.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/modules/cmap_provider/cmap_provider.h>

#include <cgogn/geometry/algos/length.h>

#include <imgui/imgui.h>

namespace cgogn
{

namespace ui
{

SurfaceRender::SurfaceRender(const App& app) :
	cgogn::ui::Module(app, "SurfaceRender"),
	selected_mesh_(nullptr),
	selected_vertex_position_(nullptr),
	selected_vertex_normal_(nullptr)
{}

SurfaceRender::~SurfaceRender()
{}

void SurfaceRender::init()
{
	cmap_provider_ = static_cast<cgogn::ui::CMapProvider*>(app_.module("CMapProvider"));
	cmap_provider_->foreach_cmap2([this] (Mesh* m, const std::string& name)
	{
		parameters_.emplace(m, Parameters(m));
	});
}

void SurfaceRender::update(const Mesh& m, const Attribute<Vec3>* vertex_position, const Attribute<Vec3>* vertex_normal)
{
	Parameters& p = parameters_[&m];

	for (cgogn::uint32 i = 0; i < 3; ++i)
	{
		p.bb_min_[i] = std::numeric_limits<cgogn::float64>::max();
		p.bb_max_[i] = std::numeric_limits<cgogn::float64>::lowest();
	}
	for (const Vec3& v : *vertex_position)
	{
		for (cgogn::uint32 i = 0; i < 3; ++i)
		{
			if (v[i] < p.bb_min_[i])
				p.bb_min_[i] = v[i];
			if (v[i] > p.bb_max_[i])
				p.bb_max_[i] = v[i];
		}
	}

	p.render_->init_primitives(m, cgogn::rendering::POINTS);
	p.render_->init_primitives(m, cgogn::rendering::LINES);
	p.render_->init_primitives(m, cgogn::rendering::TRIANGLES);

	cgogn::rendering::update_vbo(vertex_position, p.vbo_position_.get());
	if (vertex_normal)
		cgogn::rendering::update_vbo(vertex_normal, p.vbo_normal_.get());

	p.vertex_base_size_ = cgogn::geometry::mean_edge_length(m, vertex_position) / 7.0;

	p.initialized_ = true;
}

void SurfaceRender::draw(cgogn::ui::View* view)
{
	for (auto& [m, p] : parameters_)
	{
		if (!p.initialized_)
			continue;

		Vec3 diagonal = p.bb_max_ - p.bb_min_;
		view->set_scene_radius(diagonal.norm() / 2.0f);
		Vec3 center = (p.bb_max_ + p.bb_min_) / 2.0f;
		view->set_scene_center(center);

		glEnable(GL_DEPTH_TEST);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		const cgogn::rendering::GLMat4& proj_matrix = view->projection_matrix();
		const cgogn::rendering::GLMat4& view_matrix = view->modelview_matrix();

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0f, 2.0f);

		if (p.render_faces_)
		{
			if (p.phong_shading_)
			{
				p.param_phong_->bind(proj_matrix, view_matrix);
				p.render_->draw(cgogn::rendering::TRIANGLES);
				p.param_phong_->release();
			}
			else
			{
				p.param_flat_->bind(proj_matrix, view_matrix);
				p.render_->draw(cgogn::rendering::TRIANGLES);
				p.param_flat_->release();
			}
		}
		
		glDisable(GL_POLYGON_OFFSET_FILL);

		if (p.render_vertices_)
		{
			p.param_point_sprite_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
			p.param_point_sprite_->bind(proj_matrix, view_matrix);
			p.render_->draw(cgogn::rendering::POINTS);
			p.param_point_sprite_->release();
		}

		if (p.render_edges_)
		{
			p.param_edge_->bind(proj_matrix, view_matrix);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			p.render_->draw(cgogn::rendering::LINES);
			glDisable(GL_BLEND);
			p.param_edge_->release();
		}
	}
}

void SurfaceRender::interface()
{
	bool need_update = false;

	ImGui::Begin("Surface render", nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0, 0});

	if (ImGui::ListBoxHeader("Select mesh"))
	{
		cmap_provider_->foreach_cmap2([this] (Mesh* m, const std::string& name)
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
		Parameters& p = parameters_[selected_mesh_];

		std::string selected_vertex_position_name_ = selected_vertex_position_ ? selected_vertex_position_->name() : "-- select --";
		if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
		{
			cgogn::foreach_attribute<Vec3, Vertex>(*selected_mesh_, [this] (Attribute<Vec3>* attribute)
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
			cgogn::foreach_attribute<Vec3, Vertex>(*selected_mesh_, [this] (Attribute<Vec3>* attribute)
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
			if (ImGui::Button("Update data"))
				update(*selected_mesh_, selected_vertex_position_, selected_vertex_normal_);
		}

		ImGui::Separator();
		need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
		need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
		need_update |= ImGui::Checkbox("Faces", &p.render_faces_);
		if (p.render_faces_)
			need_update |= ImGui::Checkbox("Phong shading", &p.phong_shading_);

		if (p.render_faces_)
		{
			if (p.phong_shading_)
			{
				ImGui::Separator();
				ImGui::Text("Phong parameters");
				need_update |= ImGui::ColorEdit3("front color##phong", p.param_phong_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
				ImGui::SameLine();
				need_update |= ImGui::ColorEdit3("back color##phong", p.param_phong_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("spec##phong", &(p.param_phong_->specular_coef_), 10.0f, 1000.0f);
				need_update |= ImGui::Checkbox("double side##phong", &(p.param_phong_->double_side_));
			}
			else
			{
				ImGui::Separator();
				ImGui::Text("Flat parameters");
				need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
				ImGui::SameLine();
				need_update |= ImGui::ColorEdit3("back color##flat", p.param_flat_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_->double_side_));
			}
		}

		if (p.render_edges_)
		{
			ImGui::Separator();
			ImGui::Text("Edges parameters");
			need_update |= ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
			need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
		}

		if (p.render_vertices_)
		{
			ImGui::Separator();
			ImGui::Text("Vertices parameters");
			need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(), ImGuiColorEditFlags_NoInputs);
			need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1, 2.0);
		}
	}

	ImGui::End();

	if (need_update)
		for (cgogn::ui::View* v : linked_views_)
			v->request_update();
}

} // namespace ui

} // namespace cgogn
