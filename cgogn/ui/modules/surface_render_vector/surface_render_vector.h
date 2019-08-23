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

#ifndef CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_
#define CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceRenderVector : public Module
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters() :
			render_vectors_(true),
			vector_scale_factor_(1.0)
		{
			param_vector_per_vertex_ = rendering::ShaderVectorPerVertex::generate_param();
			param_vector_per_vertex_->color_ = rendering::GLColor(1, 0, 0, 1);
		}

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_vector_;

		std::unique_ptr<rendering::ShaderVectorPerVertex::Param> param_vector_per_vertex_;
		
		bool render_vectors_;
		
		float32 vector_scale_factor_;
		float32 vector_base_size_;
	};

public:

	SurfaceRenderVector(const App& app) :
		ui::Module(app, "SurfaceRenderVector (" + mesh_traits<MESH>::name + ")"),
		selected_mesh_(nullptr)
	{}

	~SurfaceRenderVector()
	{}

private:

	void init_mesh(MESH* m)
	{
		parameters_.emplace(m, Parameters());
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m] (Attribute<Vec3>* attribute)
				{
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
						p.vector_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 2.0;

					for (ui::View* v : linked_views_)
						v->request_update();
				}
			)
		);
	}

public:

	void init()
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + mesh_traits<MESH>::name + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
				mesh_provider_, this, &SurfaceRenderVector<MESH>::init_mesh
			)
		);
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vector_base_size_ = geometry::mean_edge_length(m, vertex_position.get()) / 2.0;
			md->update_vbo(vertex_position.get(), true);
		}

		p.param_vector_per_vertex_->set_vbos(md->vbo(p.vertex_position_.get()), md->vbo(p.vertex_vector_.get()));

		for (ui::View* v : linked_views_)
			v->request_update();
	}

	void set_vertex_vector(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_vector)
	{
		Parameters& p = parameters_[&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_vector_ = vertex_vector;
		if (p.vertex_vector_)
			md->update_vbo(vertex_vector.get(), true);
		
		p.param_vector_per_vertex_->set_vbos(md->vbo(p.vertex_position_.get()), md->vbo(p.vertex_vector_.get()));

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

			if (p.render_vectors_ && p.param_vector_per_vertex_->vao_initialized())
			{
				p.param_vector_per_vertex_->length_ = p.vector_base_size_ * p.vector_scale_factor_;
				p.param_vector_per_vertex_->bind(proj_matrix, view_matrix);
				md->draw(rendering::POINTS);
				p.param_vector_per_vertex_->release();
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

			if (ImGui::BeginCombo("Vector", p.vertex_vector_ ? p.vertex_vector_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == p.vertex_vector_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						set_vertex_vector(*selected_mesh_, attribute);
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
			if (p.vertex_vector_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##vector"))
					set_vertex_vector(*selected_mesh_, nullptr);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Render vectors", &p.render_vectors_);

			ImGui::Separator();
			ImGui::Text("Vector parameters");
			need_update |= ImGui::ColorEdit3("color##vectors", p.param_vector_per_vertex_->color_.data(), ImGuiColorEditFlags_NoInputs);
			need_update |= ImGui::SliderFloat("length##vectors", &(p.vector_scale_factor_), 0.1, 2.0);
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

#endif // CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_
