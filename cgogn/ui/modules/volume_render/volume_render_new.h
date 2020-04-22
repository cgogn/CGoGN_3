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

#ifndef CGOGN_MODULE_VOLUME_RENDER_H_
#define CGOGN_MODULE_VOLUME_RENDER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_explode_volumes.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/length.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class VolumeRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 3, "VolumeRender can only be used with meshes of dimension >= 3");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Volume = typename mesh_traits<MESH>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), volume_center_(), render_vertices_(true), render_edges_(true),
			  render_volumes_(false), render_volume_lines_(true), vertex_scale_factor_(1.0)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0.5f, 0, 1);

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 1, 1, 1);
			param_edge_->width_ = 1.0f;

			param_volume_ = rendering::ShaderExplodeVolumes::generate_param();
			param_volume_line_ = rendering::ShaderExplodeVolumesLine::generate_param();

			param_volume_line_->explode_ = param_volume_->explode_;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> volume_center_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderExplodeVolumes::Param> param_volume_;
		std::unique_ptr<rendering::ShaderExplodeVolumesLine::Param> param_volume_line_;

		bool render_vertices_;
		bool render_edges_;
		bool render_volumes_;
		bool render_volume_lines_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

public:
	VolumeRender(const App& app)
		: ViewModule(app, "VolumeRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~VolumeRender()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		for (View* v : linked_views_)
		{
			Parameters& p = parameters_[v][m];
			p.volume_center_ = add_attribute<Vec3, Volume>(*m, "volume_center_" + v->name());

			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);

			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, v, m]() {
					Parameters& p = parameters_[v][m];
					MeshData<MESH>* md = mesh_provider_->mesh_data(m);
					if (p.vertex_position_)
					{
						p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0;
						geometry::compute_centroid<Vec3, Volume>(*m, p.vertex_position_.get(), p.volume_center_.get());
						md->update_vbo(p.volume_center_.get(), true);
					}
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						MeshData<MESH>* md = mesh_provider_->mesh_data(m);
						if (p.vertex_position_.get() == attribute)
						{
							p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0;
							geometry::compute_centroid<Vec3, Volume>(*m, p.vertex_position_.get(),
																	 p.volume_center_.get());
							md->update_vbo(p.volume_center_.get(), true);
						}
						v->request_update();
					}));
		}
	}

public:
	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = geometry::mean_edge_length(m, vertex_position.get()) / 7.0;
			md->update_vbo(vertex_position.get(), true);

			geometry::compute_centroid<Vec3, Volume>(m, p.vertex_position_.get(), p.volume_center_.get());
			md->update_vbo(p.volume_center_.get(), true);
		}

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vc = md->vbo(p.volume_center_.get());

		p.param_point_sprite_->set_vbos({vp});
		p.param_edge_->set_vbos({vp});
		p.param_volume_->set_vbos({vp, vc});
		p.param_volume_line_->set_vbos({vp, vc});

		v.request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &VolumeRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_volumes_ && p.param_volume_->vao_initialized())
			{
				p.param_volume_->bind(proj_matrix, view_matrix);
				md->draw(rendering::VOLUMES_FACES, p.vertex_position_);
				p.param_volume_->release();

				if (p.render_volume_lines_ && p.param_volume_line_->vao_initialized())
				{
					p.param_volume_line_->bind(proj_matrix, view_matrix);
					md->draw(rendering::VOLUMES_EDGES);
					p.param_volume_line_->release();
				}
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

		if (ImGui::BeginCombo("View", selected_view_->name().c_str()))
		{
			for (View* v : linked_views_)
			{
				bool is_selected = v == selected_view_;
				if (ImGui::Selectable(v->name().c_str(), is_selected))
					selected_view_ = v;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
		}

		if (ImGui::ListBoxHeader("Mesh"))
		{
			mesh_provider_->foreach_mesh([this](MESH* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
					selected_mesh_ = m;
			});
			ImGui::ListBoxFooter();
		}

		if (selected_view_ && selected_mesh_)
		{
			double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			Parameters& p = parameters_[selected_view_][selected_mesh_];

			if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(
					*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
						bool is_selected = attribute == p.vertex_position_;
						if (ImGui::Selectable(attribute->name().c_str(), is_selected))
							set_vertex_position(*selected_view_, *selected_mesh_, attribute);
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
				ImGui::EndCombo();
			}
			if (p.vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					set_vertex_position(*selected_view_, *selected_mesh_, nullptr);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			need_update |= ImGui::Checkbox("Volumes", &p.render_volumes_);

			if (p.render_vertices_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Vertices parameters");
				need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1, 2.0);
			}

			if (p.render_edges_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Edges parameters");
				need_update |=
					ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
			}

			if (p.render_volumes_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Volumes parameters");
				if (ImGui::SliderFloat("Explode", &(p.param_volume_->explode_), 0.01f, 1.0f))
				{
					need_update = true;
					p.param_volume_line_->explode_ = p.param_volume_->explode_;
				}
				need_update |=
					ImGui::ColorEdit3("Volume color", p.param_volume_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::Checkbox("Volume lines", &p.render_volume_lines_);
				if (p.render_volume_lines_)
					need_update |= ImGui::ColorEdit3("Volume lines color", p.param_volume_line_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
			}
		}

		ImGui::End();

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
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

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
