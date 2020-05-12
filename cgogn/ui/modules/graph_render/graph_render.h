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

#ifndef CGOGN_MODULE_GRAPH_RENDER_H_
#define CGOGN_MODULE_GRAPH_RENDER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/geometry/algos/length.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class GraphRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 1, "GraphRender can only be used with meshes of dimension >= 1");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_position_vbo_(nullptr), vertex_radius_(nullptr),
			  vertex_radius_vbo_(nullptr), render_vertices_(true), render_edges_(true), vertex_scale_factor_(1.0)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = {1, 0.5f, 0, 1};

			param_point_sprite_size_ = rendering::ShaderPointSpriteSize::generate_param();
			param_point_sprite_size_->color_ = {1, 0.5f, 0, 1};

			param_bold_line_ = rendering::ShaderBoldLine::generate_param();
			param_bold_line_->color_ = {1, 1, 1, 1};
			param_bold_line_->width_ = 3.0f;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vertex_position_vbo_;
		std::shared_ptr<Attribute<Scalar>> vertex_radius_;
		rendering::VBO* vertex_radius_vbo_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderPointSpriteSize::Param> param_point_sprite_size_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_bold_line_;

		bool render_vertices_;
		bool render_edges_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

public:
	GraphRender(const App& app)
		: ViewModule(app, "GraphRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~GraphRender()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		for (View* v : linked_views_)
		{
			parameters_[v][m];
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);

			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
						{
							p.vertex_base_size_ =
								float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0);
							if (p.vertex_base_size_ == 0.0)
							{
								MeshData<MESH>* md = mesh_provider_->mesh_data(m);
								p.vertex_base_size_ = float32((md->bb_max_ - md->bb_min_).norm() / 20.0);
							}
						}
						v->request_update();
					}));
		}
	}

public:
	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_position_ == vertex_position)
			return;

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
			p.vertex_position_vbo_ = md->update_vbo(p.vertex_position_.get(), true);
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, p.vertex_position_.get()) / 7.0);
			if (p.vertex_base_size_ == 0.0)
				p.vertex_base_size_ = float32((md->bb_max_ - md->bb_min_).norm() / 20.0);
		}
		else
			p.vertex_position_vbo_ = nullptr;

		p.param_point_sprite_->set_vbos({p.vertex_position_vbo_});
		p.param_point_sprite_size_->set_vbos({p.vertex_position_vbo_, p.vertex_radius_vbo_});
		p.param_bold_line_->set_vbos({p.vertex_position_vbo_});

		v.request_update();
	}

	void set_vertex_radius(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& vertex_radius)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_radius_ == vertex_radius)
			return;

		p.vertex_radius_ = vertex_radius;
		if (p.vertex_radius_)
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
			p.vertex_radius_vbo_ = md->update_vbo(vertex_radius.get(), true);
		}
		else
			p.vertex_radius_vbo_ = nullptr;

		p.param_point_sprite_size_->set_vbos({p.vertex_position_vbo_, p.vertex_radius_vbo_});

		v.request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &GraphRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_edges_ && p.param_bold_line_->attributes_initialized())
			{
				p.param_bold_line_->bind(proj_matrix, view_matrix);
				md->draw(rendering::LINES);
				p.param_bold_line_->release();
			}

			if (p.render_vertices_)
			{
				if (p.param_point_sprite_size_->attributes_initialized())
				{
					p.param_point_sprite_size_->bind(proj_matrix, view_matrix);
					md->draw(rendering::POINTS);
					p.param_point_sprite_size_->release();
				}
				else if (p.param_point_sprite_->attributes_initialized())
				{
					p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
					p.param_point_sprite_->bind(proj_matrix, view_matrix);
					md->draw(rendering::POINTS);
					p.param_point_sprite_->release();
				}
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_, [&](MESH* m) { selected_mesh_ = m; });

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_view_, *selected_mesh_, attribute);
												});

			imgui_combo_attribute<Vertex, Scalar>(*selected_mesh_, p.vertex_radius_, "Radius",
												  [&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
													  set_vertex_radius(*selected_view_, *selected_mesh_, attribute);
												  });

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			if (p.render_vertices_)
			{
				if (p.vertex_radius_)
					need_update |= ImGui::ColorEdit3("Color##vertices", p.param_point_sprite_size_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				else
				{
					need_update |= ImGui::ColorEdit3("Color##vertices", p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("Size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
				}
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			if (p.render_edges_)
			{
				need_update |=
					ImGui::ColorEdit3("Color##edges", p.param_bold_line_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("Width##edges", &(p.param_bold_line_->width_), 1.0f, 10.0f);
			}

			if (need_update)
				for (View* v : linked_views_)
					v->request_update();
		}
	}

private:
	View* selected_view_;
	const MESH* selected_mesh_;
	std::unordered_map<View*, std::unordered_map<const MESH*, Parameters>> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_GRAPH_RENDER_H_
