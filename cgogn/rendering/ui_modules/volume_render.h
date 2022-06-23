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

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/frame_manipulator.h>
// #include <cgogn/rendering/shaders/compute_volume_centers.h>
#include <cgogn/rendering/shaders/outliner.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_explode_volumes.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_color.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_scalar.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/length.h>

#include <GLFW/glfw3.h>
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

	enum AttributePerCell
	{
		GLOBAL = 0,
		PER_VOLUME
	};
	enum ColorType
	{
		SCALAR = 0,
		VECTOR
	};

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
			: vertex_position_(nullptr), vertex_position_vbo_(nullptr), vertex_clipping_position_(nullptr),
			  volume_clipping_position_(nullptr), volume_clipping_position_vbo_(nullptr), volume_scalar_(nullptr),
			  volume_scalar_vbo_(nullptr), volume_color_(nullptr), volume_color_vbo_(nullptr), volume_center_(nullptr),
			  volume_center_vbo_(nullptr), render_vertices_(false), render_edges_(false), render_volumes_(true),
			  render_volume_lines_(true), color_per_cell_(GLOBAL), color_type_(SCALAR), vertex_scale_factor_(1.0),
			  auto_update_volume_scalar_min_max_(true), clipping_plane_(false), clip_only_volumes_(true),
			  show_frame_manipulator_(false), manipulating_frame_(false)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0.5f, 0, 1);

			param_bold_line_ = rendering::ShaderBoldLine::generate_param();
			param_bold_line_->color_ = {1.0f, 1.0f, 1.0f, 1.0f};
			param_bold_line_->width_ = 2.0f;

			param_volume_ = rendering::ShaderExplodeVolumes::generate_param();
			param_volume_->color_ = {0.4f, 0.8f, 1.0f, 1.0f};

			param_volume_line_ = rendering::ShaderExplodeVolumesLine::generate_param();
			param_volume_line_->color_ = {0.0f, 0.0f, 0.0f, 1.0f};
			param_volume_line_->explode_ = param_volume_->explode_;

			param_volume_color_ = rendering::ShaderExplodeVolumesColor::generate_param();
			param_volume_color_->explode_ = param_volume_->explode_;

			param_volume_scalar_ = rendering::ShaderExplodeVolumesScalar::generate_param();
			param_volume_scalar_->explode_ = param_volume_->explode_;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vertex_position_vbo_;

		std::shared_ptr<Attribute<Vec3>> vertex_clipping_position_;
		rendering::VBO* vertex_clipping_position_vbo_;

		std::shared_ptr<Attribute<Vec3>> volume_clipping_position_;
		rendering::VBO* volume_clipping_position_vbo_;

		std::shared_ptr<Attribute<Scalar>> volume_scalar_;
		rendering::VBO* volume_scalar_vbo_;
		std::shared_ptr<Attribute<Vec3>> volume_color_;
		rendering::VBO* volume_color_vbo_;

		std::shared_ptr<Attribute<Vec3>> volume_center_;
		rendering::VBO* volume_center_vbo_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_bold_line_;
		std::unique_ptr<rendering::ShaderExplodeVolumes::Param> param_volume_;
		std::unique_ptr<rendering::ShaderExplodeVolumesLine::Param> param_volume_line_;
		std::unique_ptr<rendering::ShaderExplodeVolumesColor::Param> param_volume_color_;
		std::unique_ptr<rendering::ShaderExplodeVolumesScalar::Param> param_volume_scalar_;

		bool render_vertices_;
		bool render_edges_;
		bool render_volumes_;
		bool render_volume_lines_;

		AttributePerCell color_per_cell_;
		ColorType color_type_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;

		bool auto_update_volume_scalar_min_max_;

		bool clipping_plane_;
		bool clip_only_volumes_;
		rendering::FrameManipulator frame_manipulator_;
		bool show_frame_manipulator_;
		bool manipulating_frame_;
	};

public:
	VolumeRender(const App& app)
		: ViewModule(app, "VolumeRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
		outline_engine_ = rendering::Outliner::instance();
		// compute_volume_center_engine_ = std::make_unique<rendering::ComputeVolumeCenterEngine>();
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

			p.volume_center_ = add_attribute<Vec3, Volume>(*m, "__volume_center");
			p.volume_clipping_position_ = add_attribute<Vec3, Volume>(*m, "__volume_clipping_position");

			std::shared_ptr<Attribute<Vec3>> vertex_position = get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);

			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, v, m]() {
					Parameters& p = parameters_[v][m];
					if (p.vertex_position_)
					{
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0);
						update_volume_center(*v, *m);
						update_volume_clipping_position(*v, *m);
					}
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
						{
							p.vertex_base_size_ =
								float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0);
							update_volume_center(*v, *m);
							update_volume_clipping_position(*v, *m);
						}
						v->request_update();
					}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Scalar>>(
					m, [this, v, m](Attribute<Scalar>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.volume_scalar_.get() == attribute)
						{
							if (p.auto_update_volume_scalar_min_max_)
								update_volume_scalar_min_max_values(p);
						}
						v->request_update();
					}));
		}
	}

public:
	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>& md = mesh_provider_->mesh_data(m);

		if (p.vertex_position_ == vertex_position)
			return;

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_position_vbo_ = md.update_vbo(vertex_position.get(), true);
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, vertex_position.get()) / 7.0);
			update_volume_center(v, m);

			if (!p.vertex_clipping_position_)
			{
				p.vertex_clipping_position_ = p.vertex_position_;
				p.vertex_clipping_position_vbo_ = md.update_vbo(p.vertex_clipping_position_.get(), true);
				update_volume_clipping_position(v, m);
			}
		}
		else
		{
			p.vertex_clipping_position_ = nullptr;

			p.vertex_position_vbo_ = nullptr;
			p.vertex_clipping_position_vbo_ = nullptr;
			p.volume_clipping_position_vbo_ = nullptr;
		}

		p.param_point_sprite_->set_vbos({p.vertex_position_vbo_, p.vertex_clipping_position_vbo_});
		p.param_bold_line_->set_vbos({p.vertex_position_vbo_, p.vertex_clipping_position_vbo_});

		p.param_volume_->set_vbos({p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_clipping_position_vbo_});
		p.param_volume_line_->set_vbos({p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_clipping_position_vbo_});
		p.param_volume_color_->set_vbos(
			{p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_color_vbo_, p.volume_clipping_position_vbo_});
		p.param_volume_scalar_->set_vbos(
			{p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_scalar_vbo_, p.volume_clipping_position_vbo_});

		Scalar size = (md.bb_max_ - md.bb_min_).norm() / 25;
		Vec3 position = 0.2 * md.bb_min_ + 0.8 * md.bb_max_;
		p.frame_manipulator_.set_size(size);
		p.frame_manipulator_.set_position(position);

		v.request_update();
	}

	void set_vertex_clipping_position(View& v, const MESH& m,
									  const std::shared_ptr<Attribute<Vec3>>& vertex_clipping_position)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>& md = mesh_provider_->mesh_data(m);

		if (p.vertex_clipping_position_ == vertex_clipping_position)
			return;

		p.vertex_clipping_position_ = vertex_clipping_position;
		if (p.vertex_clipping_position_)
		{
			p.vertex_clipping_position_vbo_ = md.update_vbo(p.vertex_clipping_position_.get(), true);
			update_volume_clipping_position(v, m);
		}
		else if (p.vertex_position_)
		{
			p.vertex_clipping_position_ = p.vertex_position_;
			p.vertex_clipping_position_vbo_ = md.update_vbo(p.vertex_clipping_position_.get(), true);
			update_volume_clipping_position(v, m);
		}
		else
		{
			p.vertex_clipping_position_vbo_ = nullptr;
			p.volume_clipping_position_vbo_ = nullptr;
		}

		p.param_point_sprite_->set_vbos({p.vertex_position_vbo_, p.vertex_clipping_position_vbo_});
		p.param_bold_line_->set_vbos({p.vertex_position_vbo_, p.vertex_clipping_position_vbo_});

		p.param_volume_->set_vbos({p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_clipping_position_vbo_});
		p.param_volume_line_->set_vbos({p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_clipping_position_vbo_});
		p.param_volume_color_->set_vbos(
			{p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_color_vbo_, p.volume_clipping_position_vbo_});
		p.param_volume_scalar_->set_vbos(
			{p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_scalar_vbo_, p.volume_clipping_position_vbo_});

		v.request_update();
	}

	void set_volume_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& volume_color)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.volume_color_ == volume_color)
			return;

		p.volume_color_ = volume_color;
		if (p.volume_color_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.volume_color_vbo_ = md.update_vbo(p.volume_color_.get(), true);
		}
		else
			p.volume_color_vbo_ = nullptr;

		p.param_volume_color_->set_vbos({p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_color_vbo_});

		v.request_update();
	}

	void set_volume_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& volume_scalar)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.volume_scalar_ == volume_scalar)
			return;

		p.volume_scalar_ = volume_scalar;
		if (p.volume_scalar_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.volume_scalar_vbo_ = md.update_vbo(p.volume_scalar_.get(), true);
			if (p.auto_update_volume_scalar_min_max_)
				update_volume_scalar_min_max_values(p);
		}
		else
		{
			p.volume_scalar_vbo_ = nullptr;
			p.param_volume_scalar_->color_map_.min_value_ = 0.0f;
			p.param_volume_scalar_->color_map_.max_value_ = 1.0f;
		}

		p.param_volume_scalar_->set_vbos({p.vertex_position_vbo_, p.volume_center_vbo_, p.volume_scalar_vbo_});
	}

protected:
	void update_volume_scalar_min_max_values(Parameters& p)
	{
		Scalar min = std::numeric_limits<float64>::max();
		Scalar max = std::numeric_limits<float64>::lowest();
		for (const Scalar& v : *p.volume_scalar_)
		{
			if (v < min)
				min = v;
			if (v > max)
				max = v;
		}
		p.param_volume_scalar_->color_map_.min_value_ = min;
		p.param_volume_scalar_->color_map_.max_value_ = max;
	}

	void update_volume_center(View& v, const MESH& m)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		if (p.vertex_position_)
		{
			geometry::compute_centroid<Vec3, Volume>(m, p.vertex_position_.get(), p.volume_center_.get());
			p.volume_center_vbo_ = md.update_vbo(p.volume_center_.get(), true);
		}
	}

	void update_volume_clipping_position(View& v, const MESH& m)
	{

		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>& md = mesh_provider_->mesh_data(m);

		if (p.vertex_clipping_position_)
		{
			geometry::compute_centroid<Vec3, Volume>(m, p.vertex_clipping_position_.get(),
													 p.volume_clipping_position_.get());
			p.volume_clipping_position_vbo_ = md.update_vbo(p.volume_clipping_position_.get(), true);
		}
	}

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &VolumeRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_volumes_)
			{
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0f, 1.5f);

				switch (p.color_per_cell_)
				{
				case GLOBAL: {
					if (p.param_volume_->attributes_initialized())
					{
						p.param_volume_->bind(proj_matrix, view_matrix);
						md.draw(rendering::VOLUMES_FACES, p.vertex_position_);
						p.param_volume_->release();
					}
				}
				break;
				case PER_VOLUME: {
					switch (p.color_type_)
					{
					case SCALAR: {
						if (p.param_volume_scalar_->attributes_initialized())
						{
							p.param_volume_scalar_->bind(proj_matrix, view_matrix);
							md.draw(rendering::VOLUMES_FACES, p.vertex_position_);
							p.param_volume_scalar_->release();
						}
					}
					break;
					case VECTOR: {
						if (p.param_volume_color_->attributes_initialized())
						{
							p.param_volume_color_->bind(proj_matrix, view_matrix);
							md.draw(rendering::VOLUMES_FACES, p.vertex_position_);
							p.param_volume_color_->release();
						}
					}
					break;
					}
				}
				break;
				}

				glDisable(GL_POLYGON_OFFSET_FILL);

				if (p.render_volume_lines_ && p.param_volume_line_->attributes_initialized())
				{
					p.param_volume_line_->bind(proj_matrix, view_matrix);
					md.draw(rendering::VOLUMES_EDGES);
					p.param_volume_line_->release();
				}
			}

			if (p.render_edges_ && p.param_bold_line_->attributes_initialized())
			{
				p.param_bold_line_->bind(proj_matrix, view_matrix);
				md.draw(rendering::LINES);
				p.param_bold_line_->release();
			}

			if (p.render_vertices_ && p.param_point_sprite_->attributes_initialized())
			{
				p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				md.draw(rendering::POINTS);
				p.param_point_sprite_->release();
			}

			if (p.show_frame_manipulator_)
				p.frame_manipulator_.draw(true, true, proj_matrix, view_matrix);

			float64 remain = md.outlined_until_ - App::frame_time_;
			if (remain > 0 && p.vertex_position_vbo_)
			{
				rendering::GLColor color{0.9f, 0.9f, 0.1f, 1};
				color *= float(remain * 2);
				if (!md.is_primitive_uptodate(rendering::TRIANGLES))
					md.init_primitives(rendering::TRIANGLES);
				outline_engine_->draw(p.vertex_position_vbo_, md.mesh_render(), proj_matrix, view_matrix, color);
			}
		}
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_C)
		{
			if (view == selected_view_ && selected_mesh_)
			{
				Parameters& p = parameters_[selected_view_][selected_mesh_];
				if (p.show_frame_manipulator_)
					p.manipulating_frame_ = true;
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_C)
		{
			if (view == selected_view_ && selected_mesh_)
			{
				Parameters& p = parameters_[selected_view_][selected_mesh_];
				p.manipulating_frame_ = false;
			}
		}
	}

	void mouse_press_event(View* view, int32, int32 x, int32 y) override
	{
		if (view == selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];
			if (p.manipulating_frame_)
			{
				auto [P, Q] = view->pixel_ray(x, y);
				p.frame_manipulator_.pick(x, y, P, Q);
				view->request_update();
			}
		}
	}

	void mouse_release_event(View* view, int32, int32, int32) override
	{
		if (view == selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];
			p.frame_manipulator_.release();
			view->request_update();
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (view == selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];
			bool leftpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_LEFT);
			bool rightpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_RIGHT);
			if (p.manipulating_frame_ && (rightpress || leftpress))
			{
				p.frame_manipulator_.drag(leftpress, x, y);
				if (p.clipping_plane_)
				{
					Vec3 position;
					p.frame_manipulator_.get_position(position);
					Vec3 axis_z;
					p.frame_manipulator_.get_axis(rendering::FrameManipulator::Zt, axis_z);
					float32 d = -(position.dot(axis_z));
					rendering::GLVec4 plane = rendering::construct_GLVec4(axis_z.x(), axis_z.y(), axis_z.z(), d);
					p.param_volume_->plane_clip_ = plane;
					p.param_volume_line_->plane_clip_ = plane;
					p.param_volume_color_->plane_clip_ = plane;
					p.param_volume_scalar_->plane_clip_ = plane;
				}
				view->stop_event();
				view->request_update();
			}
		}
	}

	void left_panel() override
	{
		bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Volume", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_view_, *selected_mesh_, attribute);
												});

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			if (p.render_vertices_)
			{
				need_update |= ImGui::ColorEdit3("Color##vertices", p.param_point_sprite_->color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("Size##vertices", &p.vertex_scale_factor_, 0.1f, 2.0f);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			if (p.render_edges_)
			{
				need_update |=
					ImGui::ColorEdit3("Color##edges", p.param_bold_line_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("Width##edges", &p.param_bold_line_->width_, 1.0f, 10.0f);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Volumes", &p.render_volumes_);
			if (p.render_volumes_)
			{
				if (ImGui::SliderFloat("Explode", &p.param_volume_->explode_, 0.01f, 1.0f))
				{
					need_update = true;
					p.param_volume_line_->explode_ = p.param_volume_->explode_;
					p.param_volume_color_->explode_ = p.param_volume_->explode_;
					p.param_volume_scalar_->explode_ = p.param_volume_->explode_;
				}

				need_update |= ImGui::Checkbox("Volume lines", &p.render_volume_lines_);
				if (p.render_volume_lines_)
					need_update |= ImGui::ColorEdit3("Volume lines color", p.param_volume_line_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);

				need_update |= ImGui::Checkbox("Apply clipping plane", &p.clipping_plane_);
				if (p.clipping_plane_)
				{
					imgui_combo_attribute<Vertex, Vec3>(
						*selected_mesh_, p.vertex_clipping_position_, "Clipping Position",
						[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
							set_vertex_clipping_position(*selected_view_, *selected_mesh_, attribute);
						});

					need_update |= ImGui::Checkbox("Clip only volumes", &p.clip_only_volumes_);

					Vec3 position;
					p.frame_manipulator_.get_position(position);
					Vec3 axis_z;
					p.frame_manipulator_.get_axis(rendering::FrameManipulator::Zt, axis_z);
					float32 d = -(position.dot(axis_z));
					rendering::GLVec4 plane = rendering::construct_GLVec4(axis_z.x(), axis_z.y(), axis_z.z(), d);

					if (p.clip_only_volumes_)
					{
						p.param_point_sprite_->plane_clip_ = {0, 0, 0, 0};
						p.param_bold_line_->plane_clip_ = {0, 0, 0, 0};
					}
					else
					{
						p.param_point_sprite_->plane_clip_ = plane;
						p.param_bold_line_->plane_clip_ = plane;
					}
					p.param_volume_->plane_clip_ = plane;
					p.param_volume_line_->plane_clip_ = plane;
					p.param_volume_color_->plane_clip_ = plane;
					p.param_volume_scalar_->plane_clip_ = plane;
				}
				else
				{
					p.param_point_sprite_->plane_clip_ = {0, 0, 0, 0};
					p.param_bold_line_->plane_clip_ = {0, 0, 0, 0};
					p.param_volume_->plane_clip_ = {0, 0, 0, 0};
					p.param_volume_line_->plane_clip_ = {0, 0, 0, 0};
					p.param_volume_color_->plane_clip_ = {0, 0, 0, 0};
					p.param_volume_scalar_->plane_clip_ = {0, 0, 0, 0};
				}

				need_update |= ImGui::Checkbox("Show clipping plane", &p.show_frame_manipulator_);
				if (p.show_frame_manipulator_)
					ImGui::TextUnformatted("Press C to manipulate the plane");

				ImGui::TextUnformatted("Colors");
				ImGui::BeginGroup();
				if (ImGui::RadioButton("Global##color", p.color_per_cell_ == GLOBAL))
				{
					p.color_per_cell_ = GLOBAL;
					need_update = true;
				}
				ImGui::SameLine();
				if (ImGui::RadioButton("Per volume##color", p.color_per_cell_ == PER_VOLUME))
				{
					p.color_per_cell_ = PER_VOLUME;
					need_update = true;
				}

				if (p.color_per_cell_ == GLOBAL)
				{
					need_update |=
						ImGui::ColorEdit3("Volume color", p.param_volume_->color_.data(), ImGuiColorEditFlags_NoInputs);
				}
				else if (p.color_per_cell_ == PER_VOLUME)
				{
					ImGui::BeginGroup();
					if (ImGui::RadioButton("Scalar", p.color_type_ == SCALAR))
					{
						p.color_type_ = SCALAR;
						need_update = true;
					}
					ImGui::SameLine();
					if (ImGui::RadioButton("Vector", p.color_type_ == VECTOR))
					{
						p.color_type_ = VECTOR;
						need_update = true;
					}
					ImGui::EndGroup();

					if (p.color_type_ == SCALAR)
					{
						imgui_combo_attribute<Volume, Scalar>(
							*selected_mesh_, p.volume_scalar_, "Attribute##scalarvolumecolor",
							[&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
								set_volume_scalar(*selected_view_, *selected_mesh_, attribute);
							});
						need_update |=
							ImGui::InputFloat("Scalar min##volumecolor", &p.param_volume_scalar_->color_map_.min_value_,
											  0.01f, 1.0f, "%.3f");
						need_update |=
							ImGui::InputFloat("Scalar max##volumecolor", &p.param_volume_scalar_->color_map_.max_value_,
											  0.01f, 1.0f, "%.3f");
						if (ImGui::Checkbox("Auto update min/max##volumecolor", &p.auto_update_volume_scalar_min_max_))
						{
							if (p.auto_update_volume_scalar_min_max_)
							{
								update_volume_scalar_min_max_values(p);
								need_update = true;
							}
						}
					}
					else if (p.color_type_ == VECTOR)
					{
						imgui_combo_attribute<Volume, Vec3>(
							*selected_mesh_, p.volume_color_, "Attribute##vectorvolumecolor",
							[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
								set_volume_color(*selected_view_, *selected_mesh_, attribute);
							});
					}
				}
				ImGui::EndGroup();
			}

			float64 remain = mesh_provider_->mesh_data(*selected_mesh_).outlined_until_ - App::frame_time_;
			if (remain > 0)
				need_update = true;

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

	rendering::Outliner* outline_engine_;
	// std::unique_ptr<rendering::ComputeVolumeCenterEngine> compute_volume_center_engine_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
