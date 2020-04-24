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

#ifndef CGOGN_MODULE_SURFACE_RENDER_H_
#define CGOGN_MODULE_SURFACE_RENDER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_color_per_vertex.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_flat_color_per_face.h>
#include <cgogn/rendering/shaders/shader_flat_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_phong_color_per_face.h>
#include <cgogn/rendering/shaders/shader_phong_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_scalar_per_vertex.h>

#include <cgogn/geometry/algos/length.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

enum AttributePerCell
{
	GLOBAL = 0,
	PER_VERTEX,
	PER_FACE
};
enum ColorType
{
	SCALAR = 0,
	VECTOR
};

template <typename MESH>
class SurfaceRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceRender can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Volume = typename mesh_traits<MESH>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_normal_(nullptr), vertex_scalar_(nullptr), vertex_color_(nullptr),
			  face_scalar_(nullptr), face_color_(nullptr), render_vertices_(false), render_edges_(false),
			  render_faces_(true), normal_per_cell_(PER_FACE), color_per_cell_(GLOBAL), color_type_(SCALAR),
			  vertex_scale_factor_(1.0), auto_update_vertex_scalar_min_max_(true)
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

			param_flat_scalar_per_vertex_ = rendering::ShaderFlatScalarPerVertex::generate_param();
			param_flat_scalar_per_vertex_->min_value_ = 0.0f;
			param_flat_scalar_per_vertex_->max_value_ = 1.0f;
			param_flat_scalar_per_vertex_->color_map_ = rendering::BWR;

			param_flat_color_per_vertex_ = rendering::ShaderFlatColorPerVertex::generate_param();

			param_flat_color_per_face_ = rendering::ShaderFlatColorPerFace::generate_param();
			param_flat_color_per_face_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_flat_scalar_per_face_ = rendering::ShaderFlatScalarPerFace::generate_param();
			param_flat_scalar_per_face_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_flat_scalar_per_face_->min_value_ = 0.0f;
			param_flat_scalar_per_face_->max_value_ = 1.0f;
			param_flat_scalar_per_face_->color_map_ = rendering::BWR;

			param_phong_ = rendering::ShaderPhong::generate_param();
			param_phong_->front_color_ = rendering::GLColor(0, 0.69f, 0.83f, 1);
			param_phong_->back_color_ = rendering::GLColor(0, 1, 0.5f, 1);
			param_phong_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_phong_->specular_coef_ = 250.0f;

			param_phong_scalar_per_vertex_ = rendering::ShaderPhongScalarPerVertex::generate_param();
			param_phong_scalar_per_vertex_->min_value_ = 0.0f;
			param_phong_scalar_per_vertex_->max_value_ = 1.0f;
			param_phong_scalar_per_vertex_->color_map_ = rendering::BWR;

			param_phong_color_per_vertex_ = rendering::ShaderPhongColorPerVertex::generate_param();
			param_phong_color_per_vertex_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_phong_color_per_face_ = rendering::ShaderPhongColorPerFace::generate_param();
			param_phong_color_per_face_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_phong_scalar_per_face_ = rendering::ShaderPhongScalarPerFace::generate_param();
			param_phong_scalar_per_face_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_phong_scalar_per_face_->min_value_ = 0.0f;
			param_phong_scalar_per_face_->max_value_ = 1.0f;
			param_phong_scalar_per_face_->color_map_ = rendering::BWR;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_normal_;
		std::shared_ptr<Attribute<Scalar>> vertex_scalar_;
		std::shared_ptr<Attribute<Vec3>> vertex_color_;
		std::shared_ptr<Attribute<Scalar>> face_scalar_;
		std::shared_ptr<Attribute<Vec3>> face_color_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderFlatScalarPerVertex::Param> param_flat_scalar_per_vertex_;
		std::unique_ptr<rendering::ShaderFlatColorPerVertex::Param> param_flat_color_per_vertex_;
		std::unique_ptr<rendering::ShaderFlatScalarPerFace::Param> param_flat_scalar_per_face_;
		std::unique_ptr<rendering::ShaderFlatColorPerFace::Param> param_flat_color_per_face_;
		std::unique_ptr<rendering::ShaderPhong::Param> param_phong_;
		std::unique_ptr<rendering::ShaderPhongScalarPerVertex::Param> param_phong_scalar_per_vertex_;
		std::unique_ptr<rendering::ShaderPhongColorPerVertex::Param> param_phong_color_per_vertex_;
		std::unique_ptr<rendering::ShaderPhongScalarPerFace::Param> param_phong_scalar_per_face_;
		std::unique_ptr<rendering::ShaderPhongColorPerFace::Param> param_phong_color_per_face_;

		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;

		AttributePerCell normal_per_cell_;
		AttributePerCell color_per_cell_;
		ColorType color_type_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;

		bool auto_update_vertex_scalar_min_max_;
		bool auto_update_face_scalar_min_max_;
	};

public:
	SurfaceRender(const App& app)
		: ViewModule(app, "SurfaceRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~SurfaceRender()
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
				boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, v, m]() {
					Parameters& p = parameters_[v][m];
					if (p.vertex_position_)
						p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0;
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
							p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7.0;
						v->request_update();
					}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Scalar>>(
					m, [this, v, m](Attribute<Scalar>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_scalar_.get() == attribute && p.auto_update_vertex_scalar_min_max_)
							update_vertex_scalar_min_max_values(p);
						if (p.face_scalar_.get() == attribute && p.auto_update_face_scalar_min_max_)
							update_face_scalar_min_max_values(p);
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
		}

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vn = md->vbo(p.vertex_normal_.get());
		auto* vs = md->vbo(p.vertex_scalar_.get());
		auto* vc = md->vbo(p.vertex_color_.get());
		auto* fs = md->vbo(p.face_scalar_.get());
		auto* fc = md->vbo(p.face_color_.get());

		p.param_point_sprite_->set_vbos({vp});
		p.param_edge_->set_vbos({vp});
		p.param_flat_->set_vbos({vp});
		p.param_flat_scalar_per_vertex_->set_vbos({vp, vs});
		p.param_flat_color_per_vertex_->set_vbos({vp, vc});
		p.param_flat_scalar_per_face_->set_vbos({vp, fs});
		p.param_flat_color_per_face_->set_vbos({vp, fc});
		p.param_phong_->set_vbos({vp, vn});
		p.param_phong_scalar_per_vertex_->set_vbos({vp, vn, vs});
		p.param_phong_color_per_vertex_->set_vbos({vp, vn, vc});
		p.param_phong_scalar_per_face_->set_vbos({vp, vn, fs});
		p.param_phong_color_per_face_->set_vbos({vp, vn, fc});

		v.request_update();
	}

	void set_vertex_normal(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_normal)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_normal_ = vertex_normal;
		if (p.vertex_normal_)
			md->update_vbo(vertex_normal.get(), true);

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vn = md->vbo(p.vertex_normal_.get());
		auto* vs = md->vbo(p.vertex_scalar_.get());
		auto* vc = md->vbo(p.vertex_color_.get());
		auto* fs = md->vbo(p.face_scalar_.get());
		auto* fc = md->vbo(p.face_color_.get());

		p.param_phong_->set_vbos({vp, vn});
		p.param_phong_scalar_per_vertex_->set_vbos({vp, vn, vs});
		p.param_phong_color_per_vertex_->set_vbos({vp, vn, vc});
		p.param_phong_scalar_per_face_->set_vbos({vp, vn, fs});
		p.param_phong_color_per_face_->set_vbos({vp, vn, fc});

		v.request_update();
	}

	void set_vertex_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& vertex_scalar)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_scalar_ = vertex_scalar;
		if (p.vertex_scalar_)
		{
			md->update_vbo(vertex_scalar.get(), true);
			if (p.auto_update_vertex_scalar_min_max_)
				update_vertex_scalar_min_max_values(p);
		}
		else
		{
			p.param_flat_scalar_per_vertex_->min_value_ = 0.0f;
			p.param_flat_scalar_per_vertex_->max_value_ = 1.0f;
			p.param_phong_scalar_per_vertex_->min_value_ = 0.0f;
			p.param_phong_scalar_per_vertex_->max_value_ = 1.0f;
		}

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vn = md->vbo(p.vertex_normal_.get());
		auto* vs = md->vbo(p.vertex_scalar_.get());

		p.param_flat_scalar_per_vertex_->set_vbos({vp, vs});
		p.param_phong_scalar_per_vertex_->set_vbos({vp, vn, vs});

		v.request_update();
	}

	void set_vertex_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_color)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_color_ = vertex_color;
		if (p.vertex_color_)
			md->update_vbo(vertex_color.get(), true);

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vn = md->vbo(p.vertex_normal_.get());
		auto* vc = md->vbo(p.vertex_color_.get());

		p.param_flat_color_per_vertex_->set_vbos({vp, vc});
		p.param_phong_color_per_vertex_->set_vbos({vp, vn, vc});

		v.request_update();
	}

	void set_face_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& face_scalar)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.face_scalar_ = face_scalar;
		if (p.face_scalar_)
		{
			md->update_vbo(face_scalar.get(), true);
			if (p.auto_update_face_scalar_min_max_)
				update_face_scalar_min_max_values(p);
		}
		else
		{
			p.param_flat_scalar_per_face_->min_value_ = 0.0f;
			p.param_flat_scalar_per_face_->max_value_ = 1.0f;
		}

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vn = md->vbo(p.vertex_normal_.get());
		auto* fs = md->vbo(p.face_scalar_.get());

		p.param_flat_scalar_per_face_->set_vbos({vp, fs});
		p.param_phong_scalar_per_face_->set_vbos({vp, vn, fs});

		v.request_update();
	}

	void set_face_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& face_color)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.face_color_ = face_color;
		if (p.face_color_)
			md->update_vbo(face_color.get(), true);

		auto* vp = md->vbo(p.vertex_position_.get());
		auto* vn = md->vbo(p.vertex_normal_.get());
		auto* fc = md->vbo(p.face_color_.get());

		p.param_flat_color_per_face_->set_vbos({vp, fc});
		p.param_phong_color_per_face_->set_vbos({vp, vn, fc});

		v.request_update();
	}

protected:
	void update_vertex_scalar_min_max_values(Parameters& p)
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
		p.param_flat_scalar_per_vertex_->min_value_ = min;
		p.param_flat_scalar_per_vertex_->max_value_ = max;
		p.param_phong_scalar_per_vertex_->min_value_ = min;
		p.param_phong_scalar_per_vertex_->max_value_ = max;
	}

	void update_face_scalar_min_max_values(Parameters& p)
	{
		Scalar min = std::numeric_limits<float64>::max();
		Scalar max = std::numeric_limits<float64>::lowest();
		for (const Scalar& v : *p.face_scalar_)
		{
			if (v < min)
				min = v;
			if (v > max)
				max = v;
		}
		p.param_flat_scalar_per_face_->min_value_ = min;
		p.param_flat_scalar_per_face_->max_value_ = max;
		p.param_phong_scalar_per_face_->min_value_ = min;
		p.param_phong_scalar_per_face_->max_value_ = max;
	}

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SurfaceRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_faces_)
			{
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0f, 2.0f);

				switch (p.normal_per_cell_)
				{
				case PER_VERTEX: {
					switch (p.color_per_cell_)
					{
					case GLOBAL: {
						if (p.param_phong_->vao_initialized())
						{
							p.param_phong_->bind(proj_matrix, view_matrix);
							md->draw(rendering::TRIANGLES, p.vertex_position_);
							p.param_phong_->release();
						}
					}
					break;
					case PER_VERTEX: {
						switch (p.color_type_)
						{
						case SCALAR: {
							if (p.param_phong_scalar_per_vertex_->vao_initialized())
							{
								p.param_phong_scalar_per_vertex_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES, p.vertex_position_);
								p.param_phong_scalar_per_vertex_->release();
							}
						}
						break;
						case VECTOR: {
							if (p.param_phong_color_per_vertex_->vao_initialized())
							{
								p.param_phong_color_per_vertex_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES, p.vertex_position_);
								p.param_phong_color_per_vertex_->release();
							}
						}
						break;
						}
					}
					break;
					case PER_FACE: {
						switch (p.color_type_)
						{
						case SCALAR: {
							if (p.param_phong_scalar_per_face_->vao_initialized())
							{
								p.param_phong_scalar_per_face_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
								p.param_phong_scalar_per_face_->release();
							}
						}
						break;
						case VECTOR: {
							if (p.param_phong_color_per_face_->vao_initialized())
							{
								p.param_phong_color_per_face_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
								p.param_phong_color_per_face_->release();
							}
						}
						break;
						}
					}
					break;
					}
				}
				break;
				case PER_FACE: {
					switch (p.color_per_cell_)
					{
					case GLOBAL: {
						if (p.param_flat_->vao_initialized())
						{
							p.param_flat_->bind(proj_matrix, view_matrix);
							md->draw(rendering::TRIANGLES, p.vertex_position_);
							p.param_flat_->release();
						}
					}
					break;
					case PER_VERTEX: {
						switch (p.color_type_)
						{
						case SCALAR: {
							if (p.param_flat_scalar_per_vertex_->vao_initialized())
							{
								p.param_flat_scalar_per_vertex_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES, p.vertex_position_);
								p.param_flat_scalar_per_vertex_->release();
							}
						}
						break;
						case VECTOR: {
							if (p.param_flat_color_per_vertex_->vao_initialized())
							{
								p.param_flat_color_per_vertex_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES, p.vertex_position_);
								p.param_flat_color_per_vertex_->release();
							}
						}
						break;
						}
					}
					break;
					case PER_FACE: {
						switch (p.color_type_)
						{
						case SCALAR: {
							if (p.param_flat_scalar_per_face_->vao_initialized())
							{
								p.param_flat_scalar_per_face_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
								p.param_flat_scalar_per_face_->release();
							}
						}
						break;
						case VECTOR: {
							if (p.param_flat_color_per_face_->vao_initialized())
							{
								p.param_flat_color_per_face_->bind(proj_matrix, view_matrix);
								md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
								p.param_flat_color_per_face_->release();
							}
						}
						break;
						}
					}
					break;
					}
				}
				break;
				}

				glDisable(GL_POLYGON_OFFSET_FILL);
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

			if (p.render_vertices_ && p.param_point_sprite_->vao_initialized())
			{
				p.param_point_sprite_->size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				md->draw(rendering::POINTS);
				p.param_point_sprite_->release();
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

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, "Position", p.vertex_position_,
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_view_, *selected_mesh_, attribute);
												});

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			if (p.render_vertices_)
			{
				need_update |= ImGui::ColorEdit3("Color##vertices", p.param_point_sprite_->color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("Size##vertices", &(p.vertex_scale_factor_), 0.1, 2.0);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);
			if (p.render_edges_)
			{
				need_update |=
					ImGui::ColorEdit3("Color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("Width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
			}

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Faces", &p.render_faces_);
			if (p.render_faces_)
			{
				ImGui::TextUnformatted("Normals");
				ImGui::BeginGroup();
				if (ImGui::RadioButton("Per vertex##normal", p.normal_per_cell_ == PER_VERTEX))
				{
					p.normal_per_cell_ = PER_VERTEX;
					need_update = true;
				}
				ImGui::SameLine();
				if (ImGui::RadioButton("Per face##normal", p.normal_per_cell_ == PER_FACE))
				{
					p.normal_per_cell_ = PER_FACE;
					need_update = true;
				}
				ImGui::EndGroup();

				if (p.normal_per_cell_ == PER_VERTEX)
				{
					imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, "Attribute##normal", p.vertex_normal_,
														[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
															set_vertex_normal(*selected_view_, *selected_mesh_,
																			  attribute);
														});
				}
				else if (p.normal_per_cell_ == PER_FACE)
				{
				}

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
				ImGui::SameLine();
				if (ImGui::RadioButton("Per face##color", p.color_per_cell_ == PER_FACE))
				{
					p.color_per_cell_ = PER_FACE;
					need_update = true;
				}
				ImGui::EndGroup();

				if (p.color_per_cell_ == GLOBAL)
				{
					if (ImGui::ColorEdit3("Front color", p.param_flat_->front_color_.data(),
										  ImGuiColorEditFlags_NoInputs))
					{
						p.param_phong_->front_color_ = p.param_flat_->front_color_;
						need_update = true;
					}
					if (ImGui::ColorEdit3("Back color", p.param_flat_->back_color_.data(),
										  ImGuiColorEditFlags_NoInputs))
					{
						p.param_phong_->back_color_ = p.param_flat_->back_color_;
						need_update = true;
					}
				}
				else if (p.color_per_cell_ == PER_VERTEX)
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
						imgui_combo_attribute<Vertex, Scalar>(
							*selected_mesh_, "Attribute##scalarvertexcolor", p.vertex_scalar_,
							[&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
								set_vertex_scalar(*selected_view_, *selected_mesh_, attribute);
							});
						if (ImGui::InputFloat("Scalar min##vertexcolor", &p.param_flat_scalar_per_vertex_->min_value_,
											  0.01f, 1.0f, "%.3f"))
						{
							p.param_phong_scalar_per_vertex_->min_value_ = p.param_flat_scalar_per_vertex_->min_value_;
							need_update = true;
						}
						if (ImGui::InputFloat("Scalar max##vertexcolor", &p.param_flat_scalar_per_vertex_->max_value_,
											  0.01f, 1.0f, "%.3f"))
						{
							p.param_phong_scalar_per_vertex_->max_value_ = p.param_flat_scalar_per_vertex_->max_value_;
							need_update = true;
						}
						if (ImGui::Checkbox("Auto update min/max##vertexcolor", &p.auto_update_vertex_scalar_min_max_))
						{
							if (p.auto_update_vertex_scalar_min_max_)
							{
								update_vertex_scalar_min_max_values(p);
								need_update = true;
							}
						}
					}
					else if (p.color_type_ == VECTOR)
					{
						imgui_combo_attribute<Vertex, Vec3>(
							*selected_mesh_, "Attribute##vectorvertexcolor", p.vertex_color_,
							[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
								set_vertex_color(*selected_view_, *selected_mesh_, attribute);
							});
					}
				}
				else if (p.color_per_cell_ == PER_FACE)
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
						imgui_combo_attribute<Face, Scalar>(
							*selected_mesh_, "Attribute##scalarfacecolor", p.face_scalar_,
							[&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
								set_face_scalar(*selected_view_, *selected_mesh_, attribute);
							});
						if (ImGui::InputFloat("Scalar min##facecolor", &p.param_flat_scalar_per_face_->min_value_,
											  0.01f, 1.0f, "%.3f"))
						{
							p.param_phong_scalar_per_face_->min_value_ = p.param_flat_scalar_per_face_->min_value_;
							need_update = true;
						}
						if (ImGui::InputFloat("Scalar max##facecolor", &p.param_flat_scalar_per_face_->max_value_,
											  0.01f, 1.0f, "%.3f"))
						{
							p.param_phong_scalar_per_face_->max_value_ = p.param_flat_scalar_per_face_->max_value_;
							need_update = true;
						}
						if (ImGui::Checkbox("Auto update min/max##facecolor", &p.auto_update_face_scalar_min_max_))
						{
							if (p.auto_update_face_scalar_min_max_)
							{
								update_face_scalar_min_max_values(p);
								need_update = true;
							}
						}
					}
					else if (p.color_type_ == VECTOR)
					{
						imgui_combo_attribute<Face, Vec3>(*selected_mesh_, "Attribute##vectorfacecolor", p.face_color_,
														  [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
															  set_face_color(*selected_view_, *selected_mesh_,
																			 attribute);
														  });
					}
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
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
