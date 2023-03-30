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

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class PointCloudRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 0, "Point Cloud Render can only be used with meshes of dimension == 0 ");

	enum AttributePerCell
	{
		GLOBAL = 0,
		PER_VERTEX
	};
	enum ColorType
	{
		VECTOR = 0
	};

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_position_vbo_(nullptr), vertex_radius_(nullptr),
			  vertex_radius_vbo_(nullptr), vertex_color_(nullptr), vertex_color_vbo_(nullptr), render_vertices_(true),
			   color_per_cell_(GLOBAL), color_type_(VECTOR), vertex_scale_factor_(1.0)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = {1.0f, 1.0f, 0.0f, 1.0f};

			param_point_sprite_size_ = rendering::ShaderPointSpriteSize::generate_param();
			param_point_sprite_size_->color_ = {1.0f, 1.0f, 0.0f, 1.0f};

			param_point_sprite_color_ = rendering::ShaderPointSpriteColor::generate_param();

			param_point_sprite_color_size_ = rendering::ShaderPointSpriteColorSize::generate_param();
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vertex_position_vbo_;
		std::shared_ptr<Attribute<Scalar>> vertex_radius_;
		rendering::VBO* vertex_radius_vbo_;
		std::shared_ptr<Attribute<Vec3>> vertex_color_;
		rendering::VBO* vertex_color_vbo_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderPointSpriteSize::Param> param_point_sprite_size_;
		std::unique_ptr<rendering::ShaderPointSpriteColor::Param> param_point_sprite_color_;
		std::unique_ptr<rendering::ShaderPointSpriteColorSize::Param> param_point_sprite_color_size_;

		bool render_vertices_;

		AttributePerCell color_per_cell_;
		ColorType color_type_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

public:
	PointCloudRender(const App& app)
		: ViewModule(app, "PointCloudRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~PointCloudRender()
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
							MeshData<MESH>& md = mesh_provider_->mesh_data(*m);
							p.vertex_base_size_ = float32((md.bb_max_ - md.bb_min_).norm() / 20.0);
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
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.vertex_base_size_ = float32((md.bb_max_ - md.bb_min_).norm() / 20.0);
			p.vertex_position_vbo_ = md.update_vbo(p.vertex_position_.get(), true);
			
		}
		else
			p.vertex_position_vbo_ = nullptr;

		p.param_point_sprite_->set_vbos({p.vertex_position_vbo_});
		p.param_point_sprite_size_->set_vbos({p.vertex_position_vbo_, p.vertex_radius_vbo_});
		p.param_point_sprite_color_->set_vbos({p.vertex_position_vbo_, p.vertex_color_vbo_});
		p.param_point_sprite_color_size_->set_vbos({p.vertex_position_vbo_, p.vertex_color_vbo_, p.vertex_radius_vbo_});

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
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.vertex_radius_vbo_ = md.update_vbo(vertex_radius.get(), true);
		}
		else
			p.vertex_radius_vbo_ = nullptr;

		p.param_point_sprite_size_->set_vbos({p.vertex_position_vbo_, p.vertex_radius_vbo_});
		p.param_point_sprite_color_size_->set_vbos({p.vertex_position_vbo_, p.vertex_color_vbo_, p.vertex_radius_vbo_});

		v.request_update();
	}

	void set_vertex_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_color)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_color_ == vertex_color)
			return;

		p.vertex_color_ = vertex_color;
		if (p.vertex_color_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.vertex_color_vbo_ = md.update_vbo(p.vertex_color_.get(), true);
		}
		else
			p.vertex_color_vbo_ = nullptr;

		p.param_point_sprite_color_->set_vbos({p.vertex_position_vbo_, p.vertex_color_vbo_});
		p.param_point_sprite_color_size_->set_vbos({p.vertex_position_vbo_, p.vertex_color_vbo_, p.vertex_radius_vbo_});

		v.request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &PointCloudRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_vertices_)
			{
				if (p.vertex_radius_)
				{
					switch (p.color_per_cell_)
					{
					case GLOBAL: {
						if (p.param_point_sprite_size_->attributes_initialized())
						{
							p.param_point_sprite_size_->bind(proj_matrix, view_matrix);
							md.draw(rendering::POINTS);
							p.param_point_sprite_size_->release();
						}
					}
					break;
					case PER_VERTEX: {
						if (p.param_point_sprite_color_size_->attributes_initialized())
						{
							p.param_point_sprite_color_size_->bind(proj_matrix, view_matrix);
							md.draw(rendering::POINTS);
							p.param_point_sprite_color_size_->release();
						}
					}
					break;
					}
				}
				else
				{
					switch (p.color_per_cell_)
					{
					case GLOBAL: {
						if (p.param_point_sprite_->attributes_initialized())
						{
							p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
							p.param_point_sprite_->bind(proj_matrix, view_matrix);
							md.draw(rendering::POINTS);
							p.param_point_sprite_->release();
						}
					}
					break;
					case PER_VERTEX: {
						if (p.param_point_sprite_color_->attributes_initialized())
						{
							p.param_point_sprite_color_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
							p.param_point_sprite_color_->bind(proj_matrix, view_matrix);
							md.draw(rendering::POINTS);
							p.param_point_sprite_color_->release();
						}
					}
					break;
					}
				}
			}
		}
	}

	void left_panel() override
	{
		bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Point_Cloud", [&](MESH& m) { selected_mesh_ = &m; });

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
				if (!p.vertex_radius_)
					need_update |= ImGui::SliderFloat("Size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);

				ImGui::TextUnformatted("Colors");
				ImGui::BeginGroup();
				if (ImGui::RadioButton("Global##color", p.color_per_cell_ == GLOBAL))
				{
					p.color_per_cell_ = GLOBAL;
					need_update = true;
				}
				ImGui::SameLine();
				if (ImGui::RadioButton("Per vertex##color", p.color_per_cell_ == PER_VERTEX))
				{
					p.color_per_cell_ = PER_VERTEX;
					need_update = true;
				}
				ImGui::EndGroup();

				if (p.color_per_cell_ == GLOBAL)
				{
					if (ImGui::ColorEdit3("Color##vertices", p.param_point_sprite_->color_.data(),
										  ImGuiColorEditFlags_NoInputs))
					{
						p.param_point_sprite_size_->color_ = p.param_point_sprite_->color_;
						need_update = true;
					}
				}
				else if (p.color_per_cell_ == PER_VERTEX)
				{
					imgui_combo_attribute<Vertex, Vec3>(
						*selected_mesh_, p.vertex_color_, "Attribute##vectorvertexcolor",
						[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
							set_vertex_color(*selected_view_, *selected_mesh_, attribute);
						});
				}
			}

			ImGui::Separator();
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
