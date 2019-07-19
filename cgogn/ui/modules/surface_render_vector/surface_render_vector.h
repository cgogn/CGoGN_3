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

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

class App;
class View;
template <typename MESH>
class MeshProvider;

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
			initialized_(false),
			vector_scale_factor_(1.0)
		{
			param_vector_per_vertex_ = rendering::ShaderVectorPerVertex::generate_param();
			param_vector_per_vertex_->color_ = rendering::GLColor(1, 0, 0, 1);
		}

		bool initialized_;

		std::string vertex_position_name_;
		std::string vertex_vector_name_;

		std::unique_ptr<rendering::ShaderVectorPerVertex::Param> param_vector_per_vertex_;
		
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

	void init()
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + mesh_traits<MESH>::name + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
		{
			parameters_.emplace(m, Parameters());
		});
	}

	void update_data(const MESH& m)
	{
		Parameters& p = parameters_[&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		md->update_vbo(p.vertex_position_name_);
		md->update_vbo(p.vertex_vector_name_);

		std::shared_ptr<Attribute<Vec3>> vertex_position = get_attribute<Vec3, Vertex>(m, p.vertex_position_name_);
		p.vector_base_size_ = geometry::mean_edge_length(m, vertex_position.get()) / 3.5;

		p.param_vector_per_vertex_->set_vbos(md->vbo(p.vertex_position_name_), md->vbo(p.vertex_vector_name_));

		p.initialized_ = true;
	}

	void set_vertex_position(const MESH& m, const std::string& vertex_position_name)
	{
		parameters_[&m].vertex_position_name_ = vertex_position_name;
	}

	void set_vertex_vector(const MESH& m, const std::string& vertex_vector_name)
	{
		parameters_[&m].vertex_vector_name_ = vertex_vector_name;
	}

protected:

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			if (!p.initialized_)
				continue;
			
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			p.param_vector_per_vertex_->length_ = p.vector_base_size_ * p.vector_scale_factor_;
			p.param_vector_per_vertex_->bind(proj_matrix, view_matrix);
			md->draw(rendering::POINTS);
			p.param_vector_per_vertex_->release();
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
			Parameters& p = parameters_[selected_mesh_];

			std::vector<Attribute<Vec3>*> vec3_attributes;
			foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (Attribute<Vec3>* attribute)
			{
				vec3_attributes.push_back(attribute);
			});

			std::string selected_vertex_position_name_ = p.vertex_position_name_.empty() ? "-- select --" : p.vertex_position_name_;
			if (ImGui::BeginCombo("Position", selected_vertex_position_name_.c_str()))
			{
				for (Attribute<Vec3>* attribute : vec3_attributes)
				{
					bool is_selected = attribute->name() == p.vertex_position_name_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						p.vertex_position_name_ = attribute->name();
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
			std::string selected_vertex_vector_name_ = p.vertex_vector_name_.empty() ? "-- select --" : p.vertex_vector_name_;
			if (ImGui::BeginCombo("Vector", selected_vertex_vector_name_.c_str()))
			{
				for (Attribute<Vec3>* attribute : vec3_attributes)
				{
					bool is_selected = attribute->name() == p.vertex_vector_name_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						p.vertex_vector_name_ = attribute->name();
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}

			if (!p.vertex_position_name_.empty() && !p.vertex_vector_name_.empty())
			{
				if (ImGui::Button("Update data"))
					update_data(*selected_mesh_);
			}

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
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_
