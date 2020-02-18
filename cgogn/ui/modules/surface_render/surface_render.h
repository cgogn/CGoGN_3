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

#ifndef CGOGN_MODULE_SURFACE_RENDER_H_
#define CGOGN_MODULE_SURFACE_RENDER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_flat_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_flat_color_per_face.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_phong_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_phong_color_per_face.h>

#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_scalar_per_vertex.h>

#include <cgogn/rendering/shaders/compute_volume_centers.h>


#include <cgogn/geometry/algos/length.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace rendering
{

static GLColor WHITE{1,1,1,1};
static GLColor BLACK{0,0,0,1};
static GLColor BLUE{0,0.69f,0.83f,1};
static GLColor GREEN{0,1,0.5f,1};
static GLColor GREY1{0.1f,0.1f,0.1f,1};
static GLColor YELLOW{0.1f,0.1f,0.1f,1};
static GLColor ORANGE{1,0.5f,0,1};

template <class SHADER, typename ...Args>
std::unique_ptr< typename SHADER::Param> generate_shader_param(Args&& ...args)
{
	std::unique_ptr< typename SHADER::Param> p = SHADER::generate_param();
	p->fill(args...);
	return p;
}

}


namespace ui
{

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

	enum RenderFaceStyle : int32
	{Nonoe=0,Flat,Flat_scalar_per_vertex,Flat_color_per_vertex,Flat_scalar_per_face,Flat_color_per_face,
	Phong,Phong_scalar_per_vertex,Phong_color_per_vertex,Phong_scalar_per_face,Phong_color_per_face};

	struct Parameters
	{
		inline Parameters() :
			vertex_position_(nullptr), vbo_vertex_position_(nullptr),
			vertex_normal_(nullptr), vbo_vertex_normal_(nullptr),
			vertex_scalar_(nullptr), vbo_vertex_scalar_(nullptr),
			vertex_color_(nullptr), vbo_vertex_color_(nullptr),
			face_scalar_(nullptr), vbo_face_scalar_(nullptr),
			face_color_(nullptr), vbo_face_color_(nullptr),
			render_vertices_(false),
			render_edges_(false),
			render_faces_style_(1),
			vertex_scale_factor_(1.0), auto_update_scalar_min_max_(true)
		{
			rendering::GLVec3 LightPosition{10,100,1000};
			using namespace rendering;
			param_point_sprite_ = generate_shader_param<ShaderPointSprite>(ORANGE,GREY1,LightPosition);
			param_edge_ = generate_shader_param<ShaderBoldLine>(WHITE,2.0f,true);
			param_flat_ = generate_shader_param<ShaderFlat>(BLUE,GREEN,GREY1,LightPosition,true);
			param_flat_spf = generate_shader_param<ShaderFlatScalarPerFace>(GREY1,LightPosition,true);
			param_flat_cpf = generate_shader_param<ShaderFlatColorPerFace>(GREY1,LightPosition,true);
			param_phong_ = generate_shader_param<ShaderPhong>(BLUE,GREEN,GREY1,WHITE,250.0f,LightPosition,true);
			param_phong_spf = generate_shader_param<ShaderPhongScalarPerFace>(GREY1,WHITE,250.0f,LightPosition,true);
			param_phong_cpf = generate_shader_param<ShaderPhongColorPerFace>(GREY1,WHITE,250.0f,LightPosition,true);
		}

		inline void update_position()
		{
			param_point_sprite_->set_vbos({vbo_vertex_position_});
			param_edge_->set_vbos({vbo_vertex_position_});
			param_flat_->set_vbos({vbo_vertex_position_});
			param_flat_spf->set_vbos({vbo_vertex_position_,vbo_face_scalar_});
			param_flat_cpf->set_vbos({vbo_vertex_position_,vbo_face_color_});
			param_phong_->set_vbos({vbo_vertex_position_, vbo_vertex_normal_});
			param_phong_spf->set_vbos({vbo_vertex_position_, vbo_vertex_normal_,vbo_face_scalar_});
			param_phong_cpf->set_vbos({vbo_vertex_position_, vbo_vertex_normal_,vbo_face_color_});
		}

		inline void update_normal()
		{
			param_phong_->set_vbos({vbo_vertex_position_, vbo_vertex_normal_});
			param_phong_spf->set_vbos({vbo_vertex_position_, vbo_vertex_normal_,vbo_face_scalar_});
			param_phong_cpf->set_vbos({vbo_vertex_position_, vbo_vertex_normal_,vbo_face_color_});
		}

		inline void update_face_scalar()
		{
			param_flat_spf->set_vbos({vbo_vertex_position_,vbo_face_scalar_});
			param_phong_spf->set_vbos({vbo_vertex_position_, vbo_vertex_normal_,vbo_face_scalar_});
		}

		inline void update_face_color()
		{
			param_flat_cpf->set_vbos({vbo_vertex_position_,vbo_face_color_});
			param_phong_cpf->set_vbos({vbo_vertex_position_, vbo_vertex_normal_,vbo_face_color_});
		}


		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vbo_vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_normal_;
		rendering::VBO* vbo_vertex_normal_;

		std::shared_ptr<Attribute<Vec3>> vertex_color_;
		rendering::VBO* vbo_vertex_color_;
		std::shared_ptr<Attribute<Scalar>> vertex_scalar_;
		rendering::VBO* vbo_vertex_scalar_;

		std::shared_ptr<Attribute<Vec3>> face_color_;
		rendering::VBO* vbo_face_color_;
		std::shared_ptr<Attribute<Scalar>> face_scalar_;
		rendering::VBO* vbo_face_scalar_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderFlatScalarPerFace::Param> param_flat_spf;
		std::unique_ptr<rendering::ShaderFlatColorPerFace::Param> param_flat_cpf;
		std::unique_ptr<rendering::ShaderPhong::Param> param_phong_;
		std::unique_ptr<rendering::ShaderPhongColorPerFace::Param> param_phong_cpf;
		std::unique_ptr<rendering::ShaderPhongScalarPerFace::Param> param_phong_spf;

//		std::unique_ptr<rendering::ShaderScalarPerVertex::Param> param_scalar_per_vertex_;
//		std::unique_ptr<rendering::ShaderScalarPerVertexGouraud::Param> param_scalar_per_vertex_gouraud_;

		bool render_vertices_;
		bool render_edges_;
		int32 render_faces_style_;


		float32 vertex_scale_factor_;
		float32 vertex_base_size_;

		bool auto_update_scalar_min_max_;
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
						if (p.vertex_scalar_.get() == attribute && p.auto_update_scalar_min_max_)
							update_scalar_min_max_values(p);
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
			p.vbo_vertex_position_ = md->update_vbo(vertex_position.get(), true);
		}

		p.update_position();

		v.request_update();
	}

	void set_vertex_normal(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_normal)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_normal_ = vertex_normal;
		p.vbo_vertex_normal_ = md->update_vbo(vertex_normal.get(), true);
		p.update_normal();

		v.request_update();
	}

//	void set_vertex_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& vertex_scalar)
//	{
//		Parameters& p = parameters_[&v][&m];
//		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
//	}

	void set_face_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& face_scal)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.face_scalar_ = face_scal;
		p.vbo_face_scalar_ = md->update_vbo(face_scal.get(), true);
		p.update_face_scalar();
		v.request_update();
	}

	void set_face_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& face_col)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.face_color_ = face_col;
		p.vbo_face_color_ = md->update_vbo(face_col.get(), true);
		p.update_face_color();
		v.request_update();
	}


protected:
	void update_scalar_min_max_values(Parameters& p)
	{
//		Scalar min = std::numeric_limits<float64>::max();
//		Scalar max = std::numeric_limits<float64>::lowest();
//		for (const Scalar& v : *p.vertex_scalar_)
//		{
//			if (v < min)
//				min = v;
//			if (v > max)
//				max = v;
//		}
//		p.param_scalar_per_vertex_->min_value_ = min;
//		p.param_scalar_per_vertex_->max_value_ = max;
//		p.param_scalar_per_vertex_gouraud_->min_value_ = min;
//		p.param_scalar_per_vertex_gouraud_->max_value_ = max;
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
//		for (auto& c : parameters_[view])
		{
//			MESH* m = c.first;
//			Parameters& p = c.second;

			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0f, 2.0f);

			switch (p.render_faces_style_)
			{
			case Flat:
				if (p.param_flat_->vao_initialized())
				{
					p.param_flat_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES, p.vertex_position_);
					p.param_flat_->release();
				}
			break;
			case Flat_scalar_per_face:
				if (p.param_flat_spf->vao_initialized())
				{
					p.param_flat_spf->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
					p.param_flat_spf->release();
				}
			break;

			case Flat_color_per_face:
				if (p.param_flat_cpf->vao_initialized())
				{
					p.param_flat_cpf->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
					p.param_flat_cpf->release();
				}
			break;
//			case Flat_scalar_per_vertex:
//				if (p.param_flat_spv->vao_initialized())
//				{
//					p.param_flat_spv->bind(proj_matrix, view_matrix);
//					md->draw(rendering::TRIANGLES, p.vertex_position_);
//					p.param_flat_spv->release();
//				}
//			break;
//			case Flat_color_per_vertex:
//				if (p.param_flat_spv->vao_initialized())
//				{
//					p.param_flat_spv->bind(proj_matrix, view_matrix);
//					md->draw(rendering::TRIANGLES, p.vertex_position_);
//					p.param_flat_spv->release();
//				}
//			break;

			case Phong:
				if (p.param_phong_->vao_initialized())
				{
					p.param_phong_->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES, p.vertex_position_);
					p.param_phong_->release();
				}
			break;
			case Phong_scalar_per_face:
				if (p.param_phong_spf->vao_initialized())
				{
					p.param_phong_spf->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
					p.param_phong_spf->release();
				}
			break;

			case Phong_color_per_face:
				if (p.param_phong_cpf->vao_initialized())
				{
					p.param_phong_cpf->bind(proj_matrix, view_matrix);
					md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
					p.param_phong_cpf->release();
				}
			break;
//			case Phong_scalar_per_vertex:
//				if (p.param_phong_spv->vao_initialized())
//				{
//					p.param_phong_spv->bind(proj_matrix, view_matrix);
//					md->draw(rendering::TRIANGLES, p.vertex_position_);
//					p.param_phong_spv->release();
//				}
//			break;
//			case Phong_color_per_vertex:
//				if (p.param_phong_spv->vao_initialized())
//				{
//					p.param_phong_spv->bind(proj_matrix, view_matrix);
//					md->draw(rendering::TRIANGLES, p.vertex_position_);
//					p.param_phong_spv->release();
//				}
//			break;
			}
			glDisable(GL_POLYGON_OFFSET_FILL);

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
			Parameters& p = parameters_[selected_view_][selected_mesh_];


			imgui_combo_attribute<Vertex,Vec3>(*selected_mesh_,p.vertex_position_, "Position", [&] ()
			{
				set_vertex_position(*selected_view_, *selected_mesh_, p.vertex_position_);
			});

			imgui_combo_attribute<Vertex,Vec3>(*selected_mesh_,p.vertex_normal_, "Normal", [&] ()
			{
				set_vertex_normal(*selected_view_, *selected_mesh_, p.vertex_normal_);
			});

			imgui_combo_attribute<Vertex,Scalar>(*selected_mesh_,p.vertex_scalar_, "Vertex Scalar", [&] ()
			{
//				set_vertex_scalar(*selected_view_, *selected_mesh_, p.vertex_scalar_);
			});

			imgui_combo_attribute<Vertex,Vec3>(*selected_mesh_,p.vertex_color_, "Vertex Color", [&] ()
			{
//				set_vertex_color(*selected_view_, *selected_mesh_, p.vertex_color_);
			});

			imgui_combo_attribute<Face,Scalar>(*selected_mesh_,p.face_scalar_, "Face Scalar", [&] ()
			{
				set_face_scalar(*selected_view_, *selected_mesh_, p.face_scalar_);
			});

			imgui_combo_attribute<Face,Vec3>(*selected_mesh_,p.face_color_, "Face Color", [&] ()
			{
				set_face_color(*selected_view_, *selected_mesh_, p.face_color_);
			});


			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);


			ImGui::BeginGroup();
			ImGui::TextUnformatted("Faces");
			need_update |= ImGui::RadioButton("Flat", &p.render_faces_style_,1);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Scalar/Vertex", &p.render_faces_style_, 2);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Color/Vertex", &p.render_faces_style_, 3);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Scalar/Face", &p.render_faces_style_, 4);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Color/Face", &p.render_faces_style_, 5);
			need_update |= ImGui::RadioButton("Phong", &p.render_faces_style_,6);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Scalar/Vertex", &p.render_faces_style_, 7);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Color/Vertex", &p.render_faces_style_, 8);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Scalar/Face", &p.render_faces_style_, 9);ImGui::SameLine();
			need_update |= ImGui::RadioButton("Color/Face", &p.render_faces_style_,10);
			ImGui::EndGroup();

			switch(p.render_faces_style_)
			{
			case Flat:
				ImGui::Separator();
				ImGui::TextUnformatted("Flat parameters");
				need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				if (p.param_flat_->double_side_)
					need_update |= ImGui::ColorEdit3("back color##flat", p.param_flat_->back_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_->double_side_));
			break;
			case Flat_scalar_per_face:
				need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_spf->double_side_));
			break;
			case Flat_color_per_face:
				need_update |= ImGui::Checkbox("double side##flat", &(p.param_flat_cpf->double_side_));
			break;
			case Phong:
				ImGui::Separator();
				ImGui::TextUnformatted("Phong parameters");
				need_update |= ImGui::ColorEdit3("front color##phong", p.param_phong_->front_color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				if (p.param_phong_->double_side_)
					need_update |= ImGui::ColorEdit3("back color##phong", p.param_phong_->back_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				need_update |=
					ImGui::SliderFloat("spec##phong", &(p.param_phong_->specular_coef_), 10.0f, 1000.0f);
				need_update |= ImGui::Checkbox("double side##phong", &(p.param_phong_->double_side_));
			break;
			case Phong_scalar_per_face:
				need_update |= ImGui::Checkbox("double side##phong", &(p.param_phong_spf->double_side_));
			break;
			case Phong_color_per_face:
				need_update |= ImGui::Checkbox("double side##phong", &(p.param_phong_cpf->double_side_));
			break;
			default:
			break;
			}

			if (p.render_edges_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Edges parameters");
				need_update |=
					ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);
			}

			if (p.render_vertices_)
			{
				ImGui::Separator();
				ImGui::TextUnformatted("Vertices parameters");
				need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
												 ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1, 2.0);
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
