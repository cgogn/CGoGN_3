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

#ifndef CGOGN_MODULE_SURFACE_SELECTION_H_
#define CGOGN_MODULE_SURFACE_SELECTION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>

#include <boost/synapse/connect.hpp>

#include <algorithm>
#include <unordered_map>

#undef near
#undef far

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceSelection : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceSelection can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Volume = typename mesh_traits<MESH>::Volume;

	enum SelectingCell
	{
		VertexSelect = 0,
		EdgeSelect,
		FaceSelect
	};

	enum SelectionMethod
	{
		SingleCell = 0,
		WithinSphere,
		FlatArea,
		ConnectedComponent
	};

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  selected_vertices_set_(nullptr), selected_edges_set_(nullptr), selected_faces_set_(nullptr),
			  selecting_cell_(VertexSelect), selection_method_(SingleCell)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_edge_->width_ = 2.0f;
			param_edge_->set_vbos({&selected_edges_vbo_});

			param_flat_ = rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_flat_->back_color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_flat_->double_side_ = true;
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_flat_->set_vbos({&selected_faces_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell(
					[&](Vertex v) { selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v)); });
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		void update_selected_edges_vbo()
		{
			if (selected_edges_set_)
			{
				std::vector<Vec3> selected_edges_position;
				selected_edges_position.reserve(selected_edges_set_->size() * 2);
				selected_edges_set_->foreach_cell([&](Edge e) {
					std::vector<Vertex> vertices = incident_vertices(*mesh_, e);
					selected_edges_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[0]));
					selected_edges_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[1]));
				});
				rendering::update_vbo(selected_edges_position, &selected_edges_vbo_);
			}
		}

		void update_selected_faces_vbo()
		{
			if (selected_faces_set_)
			{
				std::vector<Vec3> selected_faces_position;
				selected_faces_position.reserve(selected_faces_set_->size() * 3); // TODO: manage polygonal faces
				selected_faces_set_->foreach_cell([&](Face f) {
					foreach_incident_vertex(*mesh_, f, [&](Vertex v) -> bool {
						selected_faces_position.push_back(value<Vec3>(*mesh_, vertex_position_, v));
						return true;
					});
				});
				rendering::update_vbo(selected_faces_position, &selected_faces_vbo_);
			}
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;
		float32 angle_threshold_;

		rendering::VBO selected_vertices_vbo_;
		rendering::VBO selected_edges_vbo_;
		rendering::VBO selected_faces_vbo_;

		CellsSet<MESH, Vertex>* selected_vertices_set_;
		CellsSet<MESH, Edge>* selected_edges_set_;
		CellsSet<MESH, Face>* selected_faces_set_;

		SelectingCell selecting_cell_;
		SelectionMethod selection_method_;
	};

public:
	SurfaceSelection(const App& app)
		: ViewModule(app, "SurfaceSelection (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr)
	{
		param_point_sprite_selecting_sphere_ = rendering::ShaderPointSprite::generate_param();
		param_point_sprite_selecting_sphere_->color_ = rendering::GLColor(1, 1, 0, 0.25f);
		param_point_sprite_selecting_sphere_->set_vbos({&selecting_sphere_vbo_});

		param_edge_selecting_edge_ = rendering::ShaderBoldLine::generate_param();
		param_edge_selecting_edge_->color_ = rendering::GLColor(1, 1, 0, 0.25f);
		param_edge_selecting_edge_->width_ = 4.0f;
		param_edge_selecting_edge_->set_vbos({&selecting_edge_vbo_});

		param_flat_selecting_face_ = rendering::ShaderFlat::generate_param();
		param_flat_selecting_face_->front_color_ = rendering::GLColor(1, 1, 0, 0.25f);
		param_flat_selecting_face_->back_color_ = rendering::GLColor(1, 1, 0, 0.25f);
		param_flat_selecting_face_->double_side_ = true;
		param_flat_selecting_face_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
		param_flat_selecting_face_->set_vbos({&selecting_face_vbo_});
	}

	~SurfaceSelection()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6);
						p.update_selected_vertices_vbo();
						p.update_selected_edges_vbo();
						p.update_selected_faces_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Vertex>>(
				m, [this, m](CellsSet<MESH, Vertex>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Edge>>(
				m, [this, m](CellsSet<MESH, Edge>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_edges_set_ == set && p.vertex_position_)
					{
						p.update_selected_edges_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Face>>(
				m, [this, m](CellsSet<MESH, Face>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_faces_set_ == set && p.vertex_position_)
					{
						p.update_selected_faces_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, p.vertex_position_.get()) / 6); // 6 ???
			p.update_selected_vertices_vbo();
			p.update_selected_edges_vbo();
			p.update_selected_faces_vbo();
		}
		else
			selecting_ = false;

		for (View* v : linked_views_)
			v->request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SurfaceSelection<MESH>::init_mesh));
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (selected_mesh_ && view->shift_pressed())
		{
			Parameters& p = parameters_[selected_mesh_];

			if (p.vertex_position_ && (p.selecting_cell_ == VertexSelect && p.selected_vertices_set_ ||
									   p.selecting_cell_ == EdgeSelect && p.selected_edges_set_ ||
									   p.selecting_cell_ == FaceSelect && p.selected_faces_set_))
			{
				selecting_ = true;
				has_selecting_vertex_ = false;
				has_selecting_edge_ = false;
				has_selecting_face_ = false;
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		if (!view->shift_pressed())
		{
			selecting_ = false;
			has_selecting_vertex_ = false;
			has_selecting_edge_ = false;
			has_selecting_face_ = false;

			for (View* v : linked_views_)
				v->request_update();
		}
	}

	void mouse_wheel_event(View* view, int32 dx, int32 dy) override
	{
		if (selecting_)
		{
			Parameters& p = parameters_[selected_mesh_];
			if (p.selection_method_ == WithinSphere)
			{
				p.sphere_scale_factor_ += dy > 0 ? 1 : -1;
				p.sphere_scale_factor_ = std::clamp(p.sphere_scale_factor_, 10.0f, 100.0f);
				view->stop_event();
				view->request_update();
			}
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (selecting_)
		{
			Parameters& p = parameters_[selected_mesh_];

			rendering::GLVec3d near = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near.x(), near.y(), near.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			switch (p.selection_method_)
			{
			case SingleCell:
			case FlatArea:
			case ConnectedComponent: {
				switch (p.selecting_cell_)
				{
				case VertexSelect: {
					std::vector<Vertex> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
					if (!picked.empty())
					{
						selecting_vertex_ = picked[0];
						has_selecting_vertex_ = true;
						std::vector<Vec3> selecting_point = {
							value<Vec3>(*selected_mesh_, p.vertex_position_, selecting_vertex_)};
						rendering::update_vbo(selecting_point, &selecting_sphere_vbo_);
					}
					else
						has_selecting_vertex_ = false;
				}
				break;
				case EdgeSelect: {
					std::vector<Edge> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
					if (!picked.empty())
					{
						selecting_edge_ = picked[0];
						has_selecting_edge_ = true;
						std::vector<Vertex> vertices = incident_vertices(*selected_mesh_, selecting_edge_);
						std::vector<Vec3> selecting_segment = {
							value<Vec3>(*selected_mesh_, p.vertex_position_, vertices[0]),
							value<Vec3>(*selected_mesh_, p.vertex_position_, vertices[1])};
						rendering::update_vbo(selecting_segment, &selecting_edge_vbo_);
					}
					else
						has_selecting_edge_ = false;
				}
				break;
				case FaceSelect: {
					std::vector<Face> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
					if (!picked.empty())
					{
						selecting_face_ = picked[0];
						has_selecting_face_ = true;
						std::vector<Vec3> selecting_polygon;
						foreach_incident_vertex(*selected_mesh_, selecting_face_, [&](Vertex v) -> bool {
							selecting_polygon.push_back(value<Vec3>(*selected_mesh_, p.vertex_position_, v));
							return true;
						});
						rendering::update_vbo(selecting_polygon, &selecting_face_vbo_);
					}
					else
						has_selecting_face_ = false;
				}
				break;
				}
			}
			break;
			case WithinSphere: {
				std::vector<Vertex> picked;
				cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
				if (!picked.empty())
				{
					selecting_vertex_ = picked[0];
					has_selecting_vertex_ = true;
					std::vector<Vec3> selecting_point = {
						value<Vec3>(*selected_mesh_, p.vertex_position_, selecting_vertex_)};
					rendering::update_vbo(selecting_point, &selecting_sphere_vbo_);
				}
				else
					has_selecting_vertex_ = false;
			}
			break;
			}

			for (View* v : linked_views_)
				v->request_update();
		}
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selecting_)
		{
			Parameters& p = parameters_[selected_mesh_];

			switch (p.selection_method_)
			{
			case SingleCell: {
				switch (p.selecting_cell_)
				{
				case VertexSelect:
					if (has_selecting_vertex_)
					{
						switch (button)
						{
						case 0:
							p.selected_vertices_set_->select(selecting_vertex_);
							break;
						case 1:
							p.selected_vertices_set_->unselect(selecting_vertex_);
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
					}
					break;
				case EdgeSelect:
					if (has_selecting_edge_)
					{
						switch (button)
						{
						case 0:
							p.selected_edges_set_->select(selecting_edge_);
							break;
						case 1:
							p.selected_edges_set_->unselect(selecting_edge_);
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_edges_set_);
					}
					break;
				case FaceSelect:
					if (has_selecting_face_)
					{
						switch (button)
						{
						case 0:
							p.selected_faces_set_->select(selecting_face_);
							break;
						case 1:
							p.selected_faces_set_->unselect(selecting_face_);
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_faces_set_);
					}
					break;
				}
			}
			break;
			case WithinSphere: {
				if (has_selecting_vertex_)
				{
					CellCache<MESH> cache =
						geometry::within_sphere(*selected_mesh_, selecting_vertex_,
												p.vertex_base_size_ * p.sphere_scale_factor_, p.vertex_position_.get());
					switch (p.selecting_cell_)
					{
					case VertexSelect:
						switch (button)
						{
						case 0:
							foreach_cell(cache, [&p](Vertex v) -> bool {
								p.selected_vertices_set_->select(v);
								return true;
							});
							break;
						case 1:
							foreach_cell(cache, [&p](Vertex v) -> bool {
								p.selected_vertices_set_->unselect(v);
								return true;
							});
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
						break;
					case EdgeSelect:
						switch (button)
						{
						case 0:
							foreach_cell(cache, [&p](Edge e) -> bool {
								p.selected_edges_set_->select(e);
								return true;
							});
							break;
						case 1:
							foreach_cell(cache, [&p](Edge e) -> bool {
								p.selected_edges_set_->unselect(e);
								return true;
							});
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_edges_set_);
						break;
					case FaceSelect:
						switch (button)
						{
						case 0:
							foreach_cell(cache, [&p](Face f) -> bool {
								p.selected_faces_set_->select(f);
								return true;
							});
							break;
						case 1:
							foreach_cell(cache, [&p](Face f) -> bool {
								p.selected_faces_set_->unselect(f);
								return true;
							});
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_faces_set_);
						break;
					}
				}
			}
			break;
			case FlatArea: {
				switch (p.selecting_cell_)
				{
				case VertexSelect:
					if (has_selecting_vertex_)
					{
						std::vector<Vertex> vertices = geometry::within_normal_angle_threshold<Vertex>(
							*selected_mesh_, selecting_vertex_, p.angle_threshold_, p.vertex_position_.get());
						switch (button)
						{
						case 0:
							for (Vertex v : vertices)
								p.selected_vertices_set_->select(v);
							break;
						case 1:
							for (Vertex v : vertices)
								p.selected_vertices_set_->unselect(v);
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
					}
					break;
				case EdgeSelect:
					if (has_selecting_edge_)
					{
						std::vector<Edge> edges = geometry::within_normal_angle_threshold<Edge>(
							*selected_mesh_, selecting_edge_, p.angle_threshold_, p.vertex_position_.get());
						switch (button)
						{
						case 0:
							for (Edge e : edges)
								p.selected_edges_set_->select(e);
							break;
						case 1:
							for (Edge e : edges)
								p.selected_edges_set_->unselect(e);
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_edges_set_);
					}
					break;
				case FaceSelect:
					if (has_selecting_face_)
					{
						std::vector<Face> faces = geometry::within_normal_angle_threshold<Face>(
							*selected_mesh_, selecting_face_, p.angle_threshold_, p.vertex_position_.get());
						switch (button)
						{
						case 0:
							for (Face f : faces)
								p.selected_faces_set_->select(f);
							break;
						case 1:
							for (Face f : faces)
								p.selected_faces_set_->unselect(f);
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_faces_set_);
					}
					break;
				}
			}
			break;
			case ConnectedComponent: {
				switch (p.selecting_cell_)
				{
				case VertexSelect:
					if (has_selecting_vertex_)
					{
						Volume vol = incident_volumes(*selected_mesh_, selecting_vertex_)[0];
						switch (button)
						{
						case 0:
							foreach_incident_vertex(*selected_mesh_, vol, [&](Vertex v) -> bool {
								p.selected_vertices_set_->select(v);
								return true;
							});
							break;
						case 1:
							foreach_incident_vertex(*selected_mesh_, vol, [&](Vertex v) -> bool {
								p.selected_vertices_set_->unselect(v);
								return true;
							});
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
					}
					break;
				case EdgeSelect:
					if (has_selecting_edge_)
					{
						Volume vol = incident_volumes(*selected_mesh_, selecting_edge_)[0];
						switch (button)
						{
						case 0:
							foreach_incident_edge(*selected_mesh_, vol, [&](Edge e) -> bool {
								p.selected_edges_set_->select(e);
								return true;
							});
							break;
						case 1:
							foreach_incident_edge(*selected_mesh_, vol, [&](Edge e) -> bool {
								p.selected_edges_set_->unselect(e);
								return true;
							});
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_edges_set_);
					}
					break;
				case FaceSelect:
					if (has_selecting_face_)
					{
						Volume vol = incident_volumes(*selected_mesh_, selecting_face_)[0];
						switch (button)
						{
						case 0:
							foreach_incident_face(*selected_mesh_, vol, [&](Face f) -> bool {
								p.selected_faces_set_->select(f);
								return true;
							});
							break;
						case 1:
							foreach_incident_face(*selected_mesh_, vol, [&](Face f) -> bool {
								p.selected_faces_set_->unselect(f);
								return true;
							});
							break;
						}
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_faces_set_);
					}
					break;
				}
			}
			break;
			}
		}
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			// draw selected cells

			if (p.selecting_cell_ == VertexSelect && p.selected_vertices_set_ && p.selected_vertices_set_->size() > 0 &&
				p.param_point_sprite_->attributes_initialized())
			{
				p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
				p.param_point_sprite_->release();
			}
			else if (p.selecting_cell_ == EdgeSelect && p.selected_edges_set_ && p.selected_edges_set_->size() > 0 &&
					 p.param_edge_->attributes_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_LINES, 0, p.selected_edges_set_->size() * 2);
				p.param_edge_->release();
			}
			else if (p.selecting_cell_ == FaceSelect && p.selected_faces_set_ && p.selected_faces_set_->size() > 0 &&
					 p.param_flat_->attributes_initialized())
			{
				p.param_flat_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_TRIANGLES, 0, p.selected_faces_set_->size() * 3); // TODO: manage polygonal faces
				p.param_flat_->release();
			}

			// draw selection helpers

			if (selecting_)
			{
				if (p.selection_method_ == WithinSphere && has_selecting_vertex_)
				{
					param_point_sprite_selecting_sphere_->point_size_ = p.vertex_base_size_ * p.sphere_scale_factor_;
					param_point_sprite_selecting_sphere_->bind(proj_matrix, view_matrix);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glDrawArrays(GL_POINTS, 0, 1);
					glDisable(GL_BLEND);
					param_point_sprite_selecting_sphere_->release();
				}
				else // SingleCell & ConnectedComponent
				{
					if (p.selecting_cell_ == VertexSelect && has_selecting_vertex_)
					{
						param_point_sprite_selecting_sphere_->point_size_ =
							p.vertex_base_size_ * p.vertex_scale_factor_;
						param_point_sprite_selecting_sphere_->bind(proj_matrix, view_matrix);
						glEnable(GL_BLEND);
						glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
						glDrawArrays(GL_POINTS, 0, 1);
						glDisable(GL_BLEND);
						param_point_sprite_selecting_sphere_->release();
					}
					if (p.selecting_cell_ == EdgeSelect && has_selecting_edge_)
					{
						param_edge_selecting_edge_->bind(proj_matrix, view_matrix);
						glDrawArrays(GL_LINES, 0, 2);
						param_edge_selecting_edge_->release();
					}
					if (p.selecting_cell_ == FaceSelect && has_selecting_face_)
					{
						param_flat_selecting_face_->bind(proj_matrix, view_matrix);
						glDrawArrays(GL_TRIANGLES, 0, 3); // TODO: manage polygonal faces
						param_flat_selecting_face_->release();
					}
				}
			}
		}
	}

	void left_panel() override
	{
		bool need_update = false;

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			// float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
			Parameters& p = parameters_[selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_mesh_, attribute);
												});

			if (p.vertex_position_)
			{
				ImGui::Separator();
				int* ptr_sel_cell = reinterpret_cast<int*>(&p.selecting_cell_);
				need_update |= ImGui::RadioButton("Vertex", ptr_sel_cell, VertexSelect);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Edge", ptr_sel_cell, EdgeSelect);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Face", ptr_sel_cell, FaceSelect);

				ImGui::RadioButton("Single", reinterpret_cast<int*>(&p.selection_method_), SingleCell);
				ImGui::SameLine();
				ImGui::RadioButton("Sphere", reinterpret_cast<int*>(&p.selection_method_), WithinSphere);
				ImGui::SameLine();
				ImGui::RadioButton("CC", reinterpret_cast<int*>(&p.selection_method_), ConnectedComponent);
				ImGui::SameLine();
				ImGui::RadioButton("Flat", reinterpret_cast<int*>(&p.selection_method_), FlatArea);

				if (p.selection_method_ == WithinSphere)
					ImGui::SliderFloat("Sphere radius", &(p.sphere_scale_factor_), 10.0f, 100.0f);
				if (p.selection_method_ == FlatArea)
					ImGui::SliderFloat("Angle threshold", &(p.angle_threshold_), 0.0f, M_PI);

				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

				if (p.selecting_cell_ == VertexSelect)
				{
					if (ImGui::Button("Create set##vertices_set"))
						md.template add_cells_set<Vertex>();
					imgui_combo_cells_set(md, p.selected_vertices_set_, "Sets", [&](CellsSet<MESH, Vertex>* cs) {
						p.selected_vertices_set_ = cs;
						p.update_selected_vertices_vbo();
						need_update = true;
					});
					if (p.selected_vertices_set_)
					{
						ImGui::Text("(nb elements: %d)", p.selected_vertices_set_->size());
						if (ImGui::Button("Clear##vertices_set"))
						{
							p.selected_vertices_set_->clear();
							mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
						}
					}
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
				}
				else if (p.selecting_cell_ == EdgeSelect)
				{
					if (ImGui::Button("Create set##edges_set"))
						md.template add_cells_set<Edge>();
					imgui_combo_cells_set(md, p.selected_edges_set_, "Sets", [&](CellsSet<MESH, Edge>* cs) {
						p.selected_edges_set_ = cs;
						p.update_selected_edges_vbo();
						need_update = true;
					});
					if (p.selected_edges_set_)
					{
						ImGui::Text("(nb elements: %d)", p.selected_edges_set_->size());
						if (ImGui::Button("Clear##edges_set"))
						{
							p.selected_edges_set_->clear();
							mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_edges_set_);
						}
					}
					ImGui::TextUnformatted("Drawing parameters");
					need_update |=
						ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
				}
				else if (p.selecting_cell_ == FaceSelect)
				{
					if (ImGui::Button("Create set##faces_set"))
						md.template add_cells_set<Face>();
					imgui_combo_cells_set(md, p.selected_faces_set_, "Sets", [&](CellsSet<MESH, Face>* cs) {
						p.selected_faces_set_ = cs;
						p.update_selected_faces_vbo();
						need_update = true;
					});
					if (p.selected_faces_set_)
					{
						ImGui::Text("(nb elements: %d)", p.selected_faces_set_->size());
						if (ImGui::Button("Clear##faces_set"))
						{
							p.selected_faces_set_->clear();
							mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_faces_set_);
						}
					}
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	const MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;

	bool selecting_ = false;

	std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_selecting_sphere_;
	rendering::VBO selecting_sphere_vbo_;
	bool has_selecting_vertex_ = false;
	Vertex selecting_vertex_;
	std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_selecting_edge_;
	rendering::VBO selecting_edge_vbo_;
	bool has_selecting_edge_ = false;
	Edge selecting_edge_;
	std::unique_ptr<rendering::ShaderFlat::Param> param_flat_selecting_face_;
	rendering::VBO selecting_face_vbo_;
	bool has_selecting_face_ = false;
	Face selecting_face_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_SELECTION_H_
