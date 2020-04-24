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

#include <GLFW/glfw3.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/rendering/frame_manipulator.h>
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

#include <cgogn/ui/tools.h>
#include <cgogn/ui/view.h>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class Volume_Render : public ViewModule
{

	enum VBOContent : int32
	{
		Positions = 0,
		Centers,
		VolumeScalar,
		VolumeColor,
		NbVBOs
	};

	enum ShaderIndex : int32
	{
		Points = 0,
		Edges,
		Flat,
		ExplodeVolumes,
		ExplodeVolumesLines,
		ExplodeVolumeScalar,
		ExplodeVolumeColor,
		NbShaders
	};

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
			  render_vertices_(false), render_edges_(true), render_faces_(false), render_volumes_b_(true),
			  render_volumes_style(ExplodeVolumes), render_volumes_line_(true), vertex_scale_factor_(1.0),
			  vertex_base_size_(1.0), auto_update_scalar_min_max_(true), gpu_center_(false), centers_dirty_(true),
			  topo_dirty_(true)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0.5f, 0, 1);
			params_[Points] = param_point_sprite_.get();

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 1, 1, 1);
			param_edge_->width_ = 1.0f;
			params_[Edges] = param_edge_.get();

			param_flat_ = rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(0, 0.69f, 0.83f, 1);
			param_flat_->back_color_ = rendering::GLColor(0, 1, 0.5f, 1);
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			params_[Flat] = param_flat_.get();

			param_volumes_ = rendering::ShaderExplodeVolumes::generate_param();
			params_[ExplodeVolumes] = param_volumes_.get();
			param_volumes_line_ = rendering::ShaderExplodeVolumesLine::generate_param();
			params_[ExplodeVolumesLines] = param_volumes_line_.get();
			param_volumes_scalar_ = rendering::ShaderExplodeVolumesScalar::generate_param();
			params_[ExplodeVolumeScalar] = param_volumes_scalar_.get();
			param_volumes_color_ = rendering::ShaderExplodeVolumesColor::generate_param();
			params_[ExplodeVolumeColor] = param_volumes_color_.get();

			topo_drawer_ = std::make_unique<rendering::TopoDrawer>();
			topo_renderer_ = topo_drawer_->generate_renderer();

			for (auto& c : vbo_cache_)
				c = nullptr;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		inline void update_topo(const MESH& m)
		{
			if (render_topo_ && vertex_position_)
				topo_drawer_->update3D(m, vertex_position_.get());
			topo_dirty_ = false;
		}

		inline void update_positions_params()
		{
			for (int32 si = Points; si < NbShaders; ++si)
				params_[si]->set_vbo(1, vbo_cache_[Positions]);
		}

		inline void update_centers_params()
		{
			for (int32 si = ExplodeVolumes; si <= ExplodeVolumeColor; ++si)
				params_[si]->set_vbo(2, vbo_cache_[Centers]);
		}

		inline void update_gpu_centers(MeshData<MESH>* md)
		{
			if (centers_dirty_)
			{
				const MESH& m = *(md->mesh_);
				centers_dirty_ = false;
				if (vbo_cache_[Centers] == nullptr)
				{
					vbo_volume_center_ = std::make_unique<rendering::VBO>();
					auto* centers = vbo_volume_center_.get();
					vbo_cache_[Centers] = centers;
					centers->bind();
					const CMapBase& map = static_cast<const CMapBase&>(m);
					uint32 mi = is_indexed<Volume>(m) ? map.attribute_containers_[Volume::ORBIT].maximum_index()
													  : nb_cells<Volume>(m);
					centers->allocate(mi, 3);
					centers->release();
				}
			}
		}

		/// attributes
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> volume_center_;
		std::shared_ptr<Attribute<Vec3>> volume_color_;
		std::shared_ptr<Attribute<Scalar>> volume_scalar_;

		/// VBOs
		std::unique_ptr<rendering::VBO> vbo_volume_center_;
		std::array<rendering::VBO*, 4> vbo_cache_;

		/// SHaders
		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderExplodeVolumes::Param> param_volumes_;
		std::unique_ptr<rendering::ShaderExplodeVolumesScalar::Param> param_volumes_scalar_;
		std::unique_ptr<rendering::ShaderExplodeVolumesColor::Param> param_volumes_color_;
		std::unique_ptr<rendering::ShaderExplodeVolumesLine::Param> param_volumes_line_;
		std::array<rendering::ShaderParam*, 7> params_;

		/// topo
		std::unique_ptr<rendering::TopoDrawer> topo_drawer_;
		std::unique_ptr<rendering::TopoDrawer::Renderer> topo_renderer_;

		/// booleans
		bool render_topo_;
		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool render_volumes_b_;
		int render_volumes_style;
		bool render_volumes_line_;
		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		bool auto_update_scalar_min_max_;
		bool gpu_center_;

		bool centers_dirty_;
		bool topo_dirty_;
	};

public:
	Volume_Render(const App& app)
		: ViewModule(app, "Volume_Render (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr), compute_center_engine_(nullptr),
		  frame_manip_(nullptr), show_frame_manip_(false), manipullating_frame_(false)
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
						p.topo_dirty_ = true;
						p.centers_dirty_ = true;
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7);
						rendering::MeshRender* render = md->get_render();
						render->set_all_dirty();
					}
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
						{
							p.topo_dirty_ = true;
							p.centers_dirty_ = true;
							p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7);
						}
						v->request_update();
					}));

			p.param_volumes_line_->explode_ = p.param_volumes_->explode_;
			p.param_volumes_scalar_->explode_ = p.param_volumes_->explode_;
			p.param_volumes_color_->explode_ = p.param_volumes_->explode_;
			if (p.render_topo_)
			{
				p.topo_drawer_->shrink_v_ = std::min(1.0f, p.param_volumes_->explode_ + 0.02f);
				p.topo_dirty_ = true;
			}
			v->request_update();
		}
	}

	void update_gpu_centers(Parameters& p, MeshData<MESH>* md)
	{
		if (p.centers_dirty_)
		{
			const MESH& m = *(md->mesh_);
			p.centers_dirty_ = false;
			if (p.vbo_cache_[Centers] == nullptr)
			{
				p.vbo_volume_center_ = std::make_unique<rendering::VBO>();
				auto* centers = p.vbo_volume_center_.get();
				p.vbo_cache_[Centers] = centers;
				centers->bind();
				uint32 mi = is_indexed<Volume>(m) ? m.attribute_containers_[Volume::ORBIT].maximum_index()
												  : nb_cells<Volume>(m);
				centers->allocate(mi, 3);
				centers->release();
			}
		}

		auto* render = md->get_render();
		compute_center_engine_->check_primitives(render, *(md->mesh_), p.vertex_position_.get());
		compute_center_engine_->compute(p.vbo_cache_[Positions], render, p.vbo_cache_[Centers]);
		p.update_centers_params();
	}

public:
	void set_volume_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& volume_scalar)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.volume_scalar_ = volume_scalar;
		p.vbo_cache_[VolumeScalar] = md->update_vbo(volume_scalar.get(), true);
		p.param_volumes_scalar_->set_vbo(3, p.vbo_cache_[VolumeScalar]);
		v.request_update();
	}

	void set_volume_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& volume_color)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.volume_color_ = volume_color;
		p.vbo_cache_[VolumeColor] = md->update_vbo(volume_color.get(), true);
		p.param_volumes_color_->set_vbo(3, p.vbo_cache_[VolumeColor]);
		v.request_update();
	}

	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			//			p.update_topo(m);

			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, vertex_position.get()) / 7);

			p.vbo_cache_[Positions] = md->update_vbo(vertex_position.get(), true);

			p.volume_center_ = cgogn::get_attribute<Vec3, Volume>(m, "center");
			p.gpu_center_ |= (p.volume_center_ == nullptr);

			if (p.volume_center_ && !p.gpu_center_)
			{
				p.vbo_cache_[Centers] = md->update_vbo(p.volume_center_.get(), true);
			}
			if (p.gpu_center_)
			{
				p.update_gpu_centers(md);
				auto* render = md->get_render();
				compute_center_engine_->check_primitives(render, m, p.vertex_position_.get());
				compute_center_engine_->compute(p.vbo_cache_[Positions], render, p.vbo_cache_[Centers]);
				p.update_centers_params();
			}
		}
		p.update_positions_params();
		p.topo_dirty_ = true;
		p.centers_dirty_ = true;
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
		frame_manip_ = std::make_unique<cgogn::rendering::FrameManipulator>();
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
				if (p.centers_dirty_ && p.gpu_center_)
				{
					p.update_gpu_centers(md);
					auto* render = md->get_render();
					compute_center_engine_->check_primitives(render, *m, p.vertex_position_.get());
					compute_center_engine_->compute(p.vbo_cache_[Positions], render, p.vbo_cache_[Centers]);
					p.update_centers_params();
					view->request_update();
				}

				if (p.render_volumes_style == ExplodeVolumes)
				{
					if (p.param_volumes_->vao_initialized())
					{
						p.param_volumes_->bind(proj_matrix, view_matrix);
						md->draw(rendering::VOLUMES_FACES, p.vertex_position_);
						p.param_volumes_->release();
					}
				}

				if (p.render_volumes_style == ExplodeVolumeScalar && p.volume_scalar_)
				{
					if (p.param_volumes_scalar_->vao_initialized())
					{
						p.param_volumes_scalar_->bind(proj_matrix, view_matrix);
						md->draw(rendering::VOLUMES_FACES, p.vertex_position_);
						p.param_volumes_scalar_->release();
					}
				}

				if (p.render_volumes_style == ExplodeVolumeColor && p.volume_color_)
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
				md->draw(rendering::LINES);
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
				if (p.topo_dirty_)
					p.update_topo(*m);
				p.topo_renderer_->draw(proj_matrix, view_matrix);
			}
			if (show_frame_manip_)
			{
				frame_manip_->draw(true, true, proj_matrix, view_matrix);
				//				*mousePressEvent : *frame_manip_->pick(event->x(), event->y(), P, Q); // P,Q computed
				// ray 				*mouseReleaseEvent : *frame_manip_->release(); 				*mouseMouseEvent:
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		//		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		//		ImGui::SetWindowSize({0, 0});
		std::stringstream ss;
		ss << std::setw(6) << std::fixed << std::setprecision(2) << App::fps();
		std::string str_fps = ss.str() + " fps";
		ImGui::Text(str_fps.c_str());

		imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		need_update |= ImGui::Checkbox("Lock scene bb", &selected_view_->scene_bb_locked_);

		imgui_mesh_selector(mesh_provider_, selected_mesh_, [&](MESH* m) {
			selected_mesh_ = m;
			mesh_provider_->mesh_data(m)->outlined_until_ = App::frame_time_ + 1.0;
		});

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
			if (p.render_vertices_)
			{
				ImGui::SameLine();
				ImGui::TextUnformatted(" > ");
				ImGui::SameLine();
				ImGui::BeginGroup();
				need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
												 ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
				need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
				ImGui::EndGroup();
			}
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			if (p.render_edges_)
			{
				ImGui::SameLine();
				ImGui::TextUnformatted(" > ");
				ImGui::SameLine();
				ImGui::BeginGroup();
				need_update |= ImGui::ColorEdit3("edge color##edges", p.param_edge_->color_.data(),
												 ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
				need_update |= ImGui::SliderFloat("##width_edges", &(p.param_edge_->width_), 1.0f, 10.0f);
				need_update |= ImGui::SliderFloat("##lighted_edges", &(p.param_edge_->lighted_), 0.0f, 1.0f);
				ImGui::EndGroup();
			}
			need_update |= ImGui::Checkbox("Faces", &p.render_faces_);
			if (p.render_faces_)
			{
				ImGui::SameLine();
				ImGui::TextUnformatted(" > ");
				ImGui::SameLine();
				ImGui::BeginGroup();
				need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(),
												 ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
				if (p.param_flat_->double_side_)
					need_update |= ImGui::ColorEdit3("back color##flat", p.param_flat_->back_color_.data(),
													 ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
				need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_->double_side_));
				ImGui::EndGroup();
			}

			need_update |= ImGui::Checkbox("Volumes##b", &p.render_volumes_b_);
			if (p.render_volumes_b_)
			{
				if (ImGui::Checkbox("FrameManip", &show_frame_manip_))
				{
					auto [pmin, pmax] = mesh_provider_->meshes_bb();
					Scalar sz = (pmax - pmin).norm() / 5;
					Vec3 pos = 0.8 * pmin + 0.2 * pmax;
					frame_manip_->set_size(sz);
					frame_manip_->set_position(pos);
				}

				int32* ptrVS = reinterpret_cast<int32*>(&p.render_volumes_style);

				ImGui::TextUnformatted("Attribute:");
				ImGui::BeginGroup();
				need_update |= ImGui::RadioButton("None", ptrVS, ExplodeVolumes);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Scalar", ptrVS, ExplodeVolumeScalar);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color", ptrVS, ExplodeVolumeColor);
				ImGui::EndGroup();

				switch (p.render_volumes_style)
				{
				case ExplodeVolumes:
					need_update |= ImGui::ColorEdit3("color##volume", p.param_volumes_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					break;
				case ExplodeVolumeScalar:
					imgui_combo_attribute<Volume, Scalar>(*selected_mesh_, p.volume_scalar_, "Vol. Scalar Attribute",
														  [&](const decltype(p.volume_scalar_)& att) {
															  set_volume_scalar(*selected_view_, *selected_mesh_, att);
														  });
					need_update |= imgui_colormap_interface(p.param_volumes_scalar_->cm_, "vol_scal");
					break;
				case ExplodeVolumeColor:
					imgui_combo_attribute<Volume, Vec3>(*selected_mesh_, p.volume_color_, "Vol. Color Attribute",
														[&](const decltype(p.volume_color_)& att) {
															set_volume_color(*selected_view_, *selected_mesh_, att);
														});
					break;
				}

				if (ImGui::SliderFloat("explode", &(p.param_volumes_->explode_), 0.01f, 1.0f))
				{
					p.param_volumes_line_->explode_ = p.param_volumes_->explode_;
					p.param_volumes_scalar_->explode_ = p.param_volumes_->explode_;
					p.param_volumes_color_->explode_ = p.param_volumes_->explode_;
					if (p.render_topo_)
					{
						p.topo_drawer_->shrink_v_ = std::min(1.0f, p.param_volumes_->explode_ + 0.02f);
						p.topo_dirty_ = true;
					}
					need_update = true;
				}

				ImGui::Separator();
				need_update |= ImGui::Checkbox("VolumesLine", &p.render_volumes_line_);
				if (p.render_volumes_line_)
				{
					ImGui::TextUnformatted("Volume line parameters");
					need_update |= ImGui::ColorEdit3("color##volume_line", p.param_volumes_line_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				}

				ImGui::Separator();
				if (ImGui::Checkbox("Topo", &p.render_topo_))
				{
					need_update = true;
					if (p.render_topo_)
					{
						p.topo_drawer_->shrink_v_ = std::min(1.0f, p.param_volumes_->explode_ + 0.02f);
						p.topo_dirty_ = true;
					}
					p.update_topo(*selected_mesh_);
				}
			}

			if (p.render_topo_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Topo parameters");
				need_update |=
					ImGui::ColorEdit3("colorDarts", p.topo_drawer_->dart_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |=
					ImGui::ColorEdit3("colorPhi2", p.topo_drawer_->phi2_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |=
					ImGui::ColorEdit3("colorPhi3", p.topo_drawer_->phi3_color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("width##topo", &(p.topo_renderer_->width_), 1.0f, 5.0f);
				if (ImGui::SliderFloat("explodeEdges", &(p.topo_drawer_->shrink_e_), 0.01f, 1.0f))
				{
					need_update = true;
					p.topo_dirty_ = true;
				}

				if (ImGui::SliderFloat("explodeFaces", &(p.topo_drawer_->shrink_f_), 0.01f, 1.0f))
				{
					need_update = true;
					p.topo_dirty_ = true;
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

	inline void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_F)
		{
			manipullating_frame_ = true;
		}
	}

	inline void key_release_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_F)
		{
			manipullating_frame_ = false;
			view->stop_event();
		}
	}
	inline void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (manipullating_frame_)
		{
			auto [P, Q] = view->pixel_ray(x, y);
			frame_manip_->pick(x, y, P, Q); // P,Q computed ray
			view->stop_event();
			view->request_update();
		}
	}
	inline void mouse_release_event(View* view, int32 button, int32 x, int32 y) override
	{
		unused_parameters(view, button, x, y);
		if (manipullating_frame_)
		{
			frame_manip_->release();
			view->stop_event();
			view->request_update();
		}
	}

	inline void mouse_move_event(View* view, int32 buttons, int32 x, int32 y) override
	{
		unused_parameters(view, buttons);
		if (manipullating_frame_ && buttons & 3)
		{
			frame_manip_->drag(buttons & 2, x, y);
			view->stop_event();
			view->request_update();
		}
	}

private:
	View* selected_view_;
	const MESH* selected_mesh_;
	const MeshData<MESH>* selected_mesh_data_;
	std::unordered_map<View*, std::unordered_map<const MESH*, Parameters>> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
	std::unique_ptr<rendering::ComputeCenterEngine> compute_center_engine_;
	std::unique_ptr<rendering::FrameManipulator> frame_manip_;
	bool show_frame_manip_;
	bool manipullating_frame_;
};

} // namespace ui

} // namespace cgogn

#endif
