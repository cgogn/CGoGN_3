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

#ifndef CGOGN_MODULE_SURFACE_RENDER_H_
#define CGOGN_MODULE_SURFACE_RENDER_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceRender : public Module
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
			vertex_normal_(nullptr),
			render_vertices_(false),
			render_edges_(false),
			render_faces_(true),
			phong_shading_(false),
			vertex_scale_factor_(1.0)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 1);

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 1, 0, 1);
			param_edge_->width_= 2.5f;

			param_flat_ =  rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(0, 0.8f, 0, 1);
			param_flat_->back_color_ = rendering::GLColor(0, 0, 0.8f, 1);
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_phong_ = rendering::ShaderPhong::generate_param();
			param_phong_->front_color_ = rendering::GLColor(0, 0.8f, 0, 1);
			param_phong_->back_color_ = rendering::GLColor(0, 0, 0.8f, 1);
			param_phong_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_phong_->specular_coef_ = 250.0f;
		}

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_normal_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderPhong::Param> param_phong_;

		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool phong_shading_;
		
		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

public:

	SurfaceRender(const App& app) :
		ui::Module(app, "SurfaceRender (" + mesh_traits<MESH>::name + ")"),
		selected_mesh_(nullptr)
	{}

	~SurfaceRender()
	{}

private:

	void init_mesh(MESH* m)
	{
		parameters_.emplace(m, Parameters());
		mesh_connections_[m].push_back(boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(m, [this, m] (Attribute<Vec3>* attribute)
		{
			// Parameters& p = parameters_[m];
			// MeshData<MESH>* md = mesh_provider_->mesh_data(m);
			// if (p.vertex_position_.get() == attribute || p.vertex_normal_.get() == attribute)
			// 	md->update_vbo(attribute);

			for (ui::View* v : linked_views_)
				v->request_update();
		}));
	}

public:

	void init()
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + mesh_traits<MESH>::name + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(mesh_provider_, this, &SurfaceRender<MESH>::init_mesh)
		);
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = geometry::mean_edge_length(m, vertex_position.get()) / 7.0;
			md->update_vbo(vertex_position.get(), true);
		}

		p.param_point_sprite_->set_vbos(md->vbo(p.vertex_position_.get()));
		p.param_edge_->set_vbos(md->vbo(p.vertex_position_.get()));
		p.param_flat_->set_vbos(md->vbo(p.vertex_position_.get()));
		p.param_phong_->set_vbos(md->vbo(p.vertex_position_.get()), md->vbo(p.vertex_normal_.get()));

		for (ui::View* v : linked_views_)
			v->request_update();
	}

	void set_vertex_normal(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_normal)
	{
		Parameters& p = parameters_[&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_normal_ = vertex_normal;
		if (p.vertex_normal_)
			md->update_vbo(vertex_normal.get(), true);
		
		p.param_phong_->set_vbos(md->vbo(p.vertex_position_.get()), md->vbo(p.vertex_normal_.get()));

		for (ui::View* v : linked_views_)
			v->request_update();
	}

protected:

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_faces_)
			{
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0f, 2.0f);

				if (p.phong_shading_ && p.param_phong_->vao_initialized())
				{
					p.param_phong_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES);
					p.param_phong_->release();
				}
				else if (p.param_flat_->vao_initialized())
				{
					p.param_flat_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES);
					p.param_flat_->release();
				}
			
				glDisable(GL_POLYGON_OFFSET_FILL);
			}

			if (p.render_vertices_ && p.param_point_sprite_->vao_initialized())
			{
				p.param_point_sprite_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				md->draw(rendering::POINTS);
				p.param_point_sprite_->release();
			}

			if (p.render_edges_ && p.param_edge_->vao_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				md->draw(rendering::LINES);
				glDisable(GL_BLEND);
				p.param_edge_->release();
			}
		}
	}

    void interface() override
	{
		bool need_update = false;

		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Select mesh"))
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
			}

			if (ImGui::BeginCombo("Normal", p.vertex_normal_ ? p.vertex_normal_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == p.vertex_normal_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						set_vertex_normal(*selected_mesh_, attribute);
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (p.vertex_normal_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##normal"))
					set_vertex_normal(*selected_mesh_, nullptr);
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
			for (ui::View* v : linked_views_)
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

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
