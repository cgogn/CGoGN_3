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

#ifndef CGOGN_MODULE_VOLUME_RENDER_H_
#define CGOGN_MODULE_VOLUME_RENDER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/tools.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/rendering/shaders/compute_volume_centers.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_explode_volumes.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_color.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_scalar.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/topo_drawer.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class Volume_Render : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "Volume_Render can only be used with meshes of dimension >= 2");

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
			: vertex_position_(nullptr), volume_center_(nullptr), vbo_volume_center_(nullptr), render_topo_(false),
			  render_vertices_(false), render_edges_(false), render_faces_(false), render_volumes_b_(true),
			  render_volumes_(1), render_volumes_line_(true), auto_update_scalar_min_max_(true), gpu_center_(false),
			  edge_blending_(true)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0.5f, 0, 1);

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 1, 1, 1);
			param_edge_->width_ = 1.0f;

			param_flat_ = rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(0, 0.69f, 0.83f, 1);
			param_flat_->back_color_ = rendering::GLColor(0, 1, 0.5f, 1);
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_volumes_ = rendering::ShaderExplodeVolumes::generate_param();
			param_volumes_line_ = rendering::ShaderExplodeVolumesLine::generate_param();
			param_volumes_scalar_ = rendering::ShaderExplodeVolumesScalar::generate_param();
			param_volumes_color_ = rendering::ShaderExplodeVolumesColor::generate_param();

			topo_drawer_ = std::make_unique<rendering::TopoDrawer>();
			topo_renderer_ = topo_drawer_->generate_renderer();
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		inline void update_topo(const MESH& m)
		{
			if (render_topo_ && vertex_position_)
				topo_drawer_->update3D(m, vertex_position_.get());
		}

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> volume_center_;
		std::shared_ptr<Attribute<Vec3>> volume_color_;
		std::shared_ptr<Attribute<Scalar>> volume_scalar_;
		std::unique_ptr<rendering::VBO> vbo_volume_center_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderExplodeVolumes::Param> param_volumes_;
		std::unique_ptr<rendering::ShaderExplodeVolumesScalar::Param> param_volumes_scalar_;
		std::unique_ptr<rendering::ShaderExplodeVolumesColor::Param> param_volumes_color_;
		std::unique_ptr<rendering::ShaderExplodeVolumesLine::Param> param_volumes_line_;
		rendering::VBO* vbo_center_;

		std::unique_ptr<rendering::TopoDrawer> topo_drawer_;
		std::unique_ptr<rendering::TopoDrawer::Renderer> topo_renderer_;
		bool render_topo_;

		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool render_volumes_b_;
		int render_volumes_;
		bool render_volumes_line_;
		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		bool auto_update_scalar_min_max_;
		bool gpu_center_;
		bool edge_blending_;
	};

public:
	Volume_Render(const App& app)
		: ViewModule(app, "Volume_Render (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr), compute_center_engine_(nullptr)
	{
	}

	~Volume_Render()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		for (View* v : linked_views_)
		{
			auto& p = parameters_[v][m];
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);

			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, v, m]() {
					Parameters& p = parameters_[v][m];
					if (p.vertex_position_)
					{
						MeshData<MESH>* md = mesh_provider_->mesh_data(m);
						p.update_topo(*m);
						p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0;
						rendering::MeshRender* render = md->get_render();
						render->set_all_dirty();
						if (!render->is_primitive_uptodate(rendering::VOLUMES_VERTICES))
							render->init_primitives(*m, rendering::VOLUMES_VERTICES, p.vertex_position_.get());
						p.vbo_center_->bind();
						uint32 mi = is_indexed<Volume>(*m) ? m->attribute_containers_[Volume::ORBIT].maximum_index() + 1
														   : nb_cells<Volume>(*m);
						p.vbo_center_->allocate(mi, 3);
						p.vbo_center_->release();
						compute_center_engine_->compute(md->vbo(p.vertex_position_.get()), render, p.vbo_center_);
					}
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
						{
							p.update_topo(*m);
							p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0;
						}
						v->request_update();
					}));

			p.param_volumes_line_->explode_ = p.param_volumes_->explode_;
			p.param_volumes_scalar_->explode_ = p.param_volumes_->explode_;
			p.param_volumes_color_->explode_ = p.param_volumes_->explode_;
			if (p.render_topo_)
			{
				p.topo_drawer_->shrink_v_ = std::min(1.0, p.param_volumes_->explode_ + 0.02);
				p.update_topo(*selected_mesh_);
			}
			v->request_update();
		}
	}

public:
	void set_volume_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& volume_scalar)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.volume_scalar_ = volume_scalar;
		if (p.volume_scalar_)
		{
			md->update_vbo(volume_scalar.get(), true);
		}

		p.param_volumes_scalar_->set_vbos(
			{md->vbo(p.vertex_position_.get()), p.vbo_center_, md->vbo(p.volume_scalar_.get())});
		v.request_update();
	}

	void set_volume_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& volume_color)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.volume_color_ = volume_color;
		if (p.volume_color_)
		{
			md->update_vbo(volume_color.get(), true);
		}

		p.param_volumes_color_->set_vbos(
			{md->vbo(p.vertex_position_.get()), p.vbo_center_, md->vbo(p.volume_color_.get())});
		v.request_update();
	}

	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_position_ = vertex_position;
		p.vbo_center_ = nullptr;
		if (p.vertex_position_)
		{
			p.update_topo(m);

			p.vertex_base_size_ = geometry::mean_edge_length(m, vertex_position.get()) / 7.0;
			md->update_vbo(vertex_position.get(), true);

			p.volume_center_ = cgogn::get_attribute<Vec3, Volume>(m, "center");
			p.gpu_center_ |= (p.volume_center_ == nullptr);

			if (p.volume_center_ && !p.gpu_center_)
			{
				md->update_vbo(p.volume_center_.get(), true);
				p.vbo_center_ = md->vbo(p.volume_center_.get());
			}
		}

		if (p.gpu_center_)
		{
			if (p.vbo_center_ == nullptr)
			{
				p.vbo_volume_center_ = std::make_unique<rendering::VBO>();
				p.vbo_center_ = p.vbo_volume_center_.get();
				p.vbo_center_->bind();
				uint32 mi = is_indexed<Volume>(m) ? m.attribute_containers_[Volume::ORBIT].maximum_index()
												  : nb_cells<Volume>(m);
				p.vbo_center_->allocate(mi, 3);
				p.vbo_center_->release();
			}
			rendering::MeshRender* render = md->get_render();
			if (!render->is_primitive_uptodate(rendering::VOLUMES_VERTICES))
				render->init_primitives(m, rendering::VOLUMES_VERTICES, p.vertex_position_.get());
			compute_center_engine_->compute(md->vbo(p.vertex_position_.get()), render, p.vbo_center_);
		}

		auto* vp = md->vbo(p.vertex_position_.get());
		p.param_point_sprite_->set_vbos({vp});
		p.param_edge_->set_vbos({vp});
		p.param_flat_->set_vbos({vp});
		p.param_volumes_->set_vbos({vp, p.vbo_center_});
		p.param_volumes_line_->set_vbos({vp, p.vbo_center_});

		v.request_update();
	}

protected:
	void update_scalar_min_max_values(Parameters& p)
	{
		Scalar min = std::numeric_limits<float64>::max();
		Scalar max = std::numeric_limits<float64>::lowest();
		for (const Scalar& v : *p.vertex_scalar_)
		{
			if (v < min)
				min = v;
			if (v > max)
				max = v;
		}
		p.param_scalar_per_vertex_->min_value_ = min;
		p.param_scalar_per_vertex_->max_value_ = max;
		p.param_scalar_per_vertex_gouraud_->min_value_ = min;
		p.param_scalar_per_vertex_gouraud_->max_value_ = max;
	}

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &Volume_Render<MESH>::init_mesh));
		compute_center_engine_ = std::make_unique<rendering::ComputeCenterEngine>();
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0f, 1.5f);

			if (p.render_faces_)
			{

				if (p.param_flat_->vao_initialized())
				{
					p.param_flat_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES, p.vertex_position_);
					p.param_flat_->release();
				}
			}

			if (p.render_volumes_b_)
			{
				//				md->update_vbo(p.vertex_position_.get());
				compute_center_engine_->compute(md->vbo(p.vertex_position_.get()), md->get_render(), p.vbo_center_);

				if (p.render_volumes_ == 1)
				{
					if (p.param_volumes_->vao_initialized())
					{
						p.param_volumes_->bind(proj_matrix, view_matrix);
						md->draw(rendering::VOLUMES_FACES, p.vertex_position_);
						p.param_volumes_->release();
					}
				}

				if (p.render_volumes_ == 2 && p.volume_scalar_)
				{
					if (p.param_volumes_scalar_->vao_initialized())
					{
						p.param_volumes_scalar_->bind(proj_matrix, view_matrix);
						md->draw(rendering::VOLUMES_FACES, p.vertex_position_);
						p.param_volumes_scalar_->release();
					}
				}

				if (p.render_volumes_ == 3 && p.volume_color_)
				{
					if (p.param_volumes_color_->vao_initialized())
					{
						p.param_volumes_color_->bind(proj_matrix, view_matrix);
						md->draw(rendering::VOLUMES_FACES, p.vertex_position_);
						p.param_volumes_color_->release();
					}
				}
			}

			glDisable(GL_POLYGON_OFFSET_FILL);

			if (p.render_volumes_b_ && p.render_volumes_line_)
			{
				if (p.param_volumes_->vao_initialized())
				{
					p.param_volumes_line_->bind(proj_matrix, view_matrix);
					md->draw(rendering::VOLUMES_EDGES);
					p.param_volumes_line_->release();
				}
			}

			if (p.render_edges_ && p.param_edge_->vao_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				if (p.edge_blending_)
				{
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				}
				md->draw(rendering::LINES);
				glDisable(GL_BLEND);
				p.param_edge_->release();
			}

			if (p.render_vertices_ && p.param_point_sprite_->vao_initialized())
			{
				p.param_point_sprite_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				md->draw(rendering::POINTS);
				p.param_point_sprite_->release();
			}

			if (p.render_topo_)
			{
				p.topo_renderer_->draw(proj_matrix, view_matrix);
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});
		std::stringstream ss;
		ss << std::setw(6) << std::fixed << std::setprecision(2) << App::fps();
		std::string str_fps = ss.str() + " fps";
		ImGui::Text(str_fps.c_str());

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
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position Attribute",
												[&](const decltype(p.vertex_position_)& att) {
													set_vertex_position(*selected_view_, *selected_mesh_, att);
												});

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			need_update |= ImGui::Checkbox("Faces", &p.render_faces_);

			need_update |= ImGui::Checkbox("Volumes##b", &p.render_volumes_b_);
			if (p.render_volumes_b_)
			{
				ImGui::TextUnformatted("Attribute:");
				ImGui::BeginGroup();
				need_update |= ImGui::RadioButton("None", &p.render_volumes_, 1);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Scalar", &p.render_volumes_, 2);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color", &p.render_volumes_, 3);
				ImGui::EndGroup();
				need_update |= ImGui::Checkbox("VolumesLine", &p.render_volumes_line_);
				ImGui::TextUnformatted("Volume parameters");

				if (ImGui::SliderFloat("explode", &(p.param_volumes_->explode_), 0.01f, 1.0f))
				{
					p.param_volumes_line_->explode_ = p.param_volumes_->explode_;
					p.param_volumes_scalar_->explode_ = p.param_volumes_->explode_;
					p.param_volumes_color_->explode_ = p.param_volumes_->explode_;
					if (p.render_topo_)
					{
						p.topo_drawer_->shrink_v_ = std::min(1.0, p.param_volumes_->explode_ + 0.02);
						p.update_topo(*selected_mesh_);
					}
					need_update = true;
				}

				if (p.render_volumes_ == 1)
				{
					ImGui::Separator();
					need_update |= ImGui::ColorEdit3("color##volume", p.param_volumes_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				}
				if (p.render_volumes_ == 2)
				{
					ImGui::Separator();

					imgui_combo_attribute<Volume, Scalar>(*selected_mesh_, p.volume_scalar_, "Vol. Scalar Attribute",
														  [&](const decltype(p.volume_scalar_)& att) {
															  set_volume_scalar(*selected_view_, *selected_mesh_, att);
														  });
					need_update |= imgui_colormap_interface(p.param_volumes_scalar_->cm_, "vol_scal");
				}

				if (p.render_volumes_ == 3)
				{
					ImGui::Separator();
					imgui_combo_attribute<Volume, Vec3>(*selected_mesh_, p.volume_color_, "Vol. Color Attribute",
														[&](const decltype(p.volume_color_)& att) {
															set_volume_color(*selected_view_, *selected_mesh_, att);
														});
				}

				if (p.render_volumes_line_)
				{
					ImGui::Separator();
					ImGui::TextUnformatted("Volume line parameters");
					need_update |= ImGui::ColorEdit3("color##volume_line", p.param_volumes_line_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				}
			}

			ImGui::Separator();
			if (ImGui::Checkbox("Topo", &p.render_topo_))
			{
				need_update = true;
				p.update_topo(*selected_mesh_);
			}

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
				need_update |= ImGui::Checkbox("Blendinh", &p.edge_blending_);
				need_update |=
					ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
			}

			if (p.render_faces_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Flat parameters");
				need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				if (p.param_flat_->double_side_)
					need_update |= ImGui::ColorEdit3("back color##flat", p.param_flat_->back_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_->double_side_));
			}

			if (p.render_topo_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Volume parameters");
				need_update |=
					ImGui::ColorEdit3("colorDarts", p.topo_drawer_->dart_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |=
					ImGui::ColorEdit3("colorPhi2", p.topo_drawer_->phi2_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |=
					ImGui::ColorEdit3("colorPhi3", p.topo_drawer_->phi3_color_.data(), ImGuiColorEditFlags_NoInputs);
				if (ImGui::SliderFloat("explodeEdges", &(p.topo_drawer_->shrink_e_), 0.01f, 1.0f))
				{
					need_update = true;
					p.update_topo(*selected_mesh_);
				}

				if (ImGui::SliderFloat("explodeFaces", &(p.topo_drawer_->shrink_f_), 0.01f, 1.0f))
				{
					need_update = true;
					p.update_topo(*selected_mesh_);
				}
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
	std::unique_ptr<rendering::ComputeCenterEngine> compute_center_engine_;
};

} // namespace ui

} // namespace cgogn

#endif
