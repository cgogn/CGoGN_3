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
#include <cgogn/ui/tools.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_flat_color_per_face.h>
#include <cgogn/rendering/shaders/shader_flat_color_per_vertex.h>
#include <cgogn/rendering/shaders/shader_flat_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_flat_scalar_per_vertex.h>
#include <cgogn/rendering/shaders/shader_no_illum.h>
#include <cgogn/rendering/shaders/shader_no_illum_color_per_face.h>
#include <cgogn/rendering/shaders/shader_no_illum_color_per_vertex.h>
#include <cgogn/rendering/shaders/shader_no_illum_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_no_illum_scalar_per_vertex.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_phong_color_per_face.h>
#include <cgogn/rendering/shaders/shader_phong_color_per_vertex.h>
#include <cgogn/rendering/shaders/shader_phong_scalar_per_face.h>
#include <cgogn/rendering/shaders/shader_phong_scalar_per_vertex.h>

#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_scalar_per_vertex.h>
#include <cgogn/rendering/topo_drawer.h>

#include <cgogn/rendering/shaders/compute_normals.h>

#include <cgogn/geometry/algos/length.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

enum SurfaceRenderStyle : int32
{
	RS_None = -1,
	RS_Points,
	RS_Lines,
	RS_NoIllum,
	RS_NoIllum_scalar_per_vertex,
	RS_NoIllum_color_per_vertex,
	RS_NoIllum_scalar_per_face,
	RS_NoIllum_color_per_face,
	RS_Flat,
	RS_Flat_scalar_per_vertex,
	RS_Flat_color_per_vertex,
	RS_Flat_scalar_per_face,
	RS_Flat_color_per_face,
	RS_Phong,
	RS_Phong_scalar_per_vertex,
	RS_Phong_color_per_vertex,
	RS_Phong_scalar_per_face,
	RS_Phong_color_per_face
};

inline SurfaceRenderStyle& operator++(SurfaceRenderStyle& d)
{
	++*reinterpret_cast<int32*>(&d);
	return d;
}

template <typename MESH>
class SurfaceRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceRender can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using AttributeGen = typename mesh_traits<MESH>::AttributeGen;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Volume = typename mesh_traits<MESH>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	constexpr bool is_per_face(SurfaceRenderStyle rs)
	{
		return ((rs > 1) && ((rs % 5) < 2));
	}

	enum VBOContent : int32
	{
		VBO_Position = 0,
		VBO_Normal,
		VBO_Scalar_per_vertex,
		VBO_Vec3_per_vertex,
		VBO_Scalar_per_face,
		VBO_Vec3_per_face,
	};

	// using ILSurfaceRenderStyles = std::initializer_list<SurfaceRenderStyle>;
	using ILSurfaceRenderStyles = std::vector<SurfaceRenderStyle>;

	static constexpr ILSurfaceRenderStyles style_with_vbo(VBOContent v)
	{
		switch (v)
		{
		case VBO_Position:
			return ILSurfaceRenderStyles{RS_Points,
										 RS_Lines,
										 RS_NoIllum,
										 RS_NoIllum_scalar_per_vertex,
										 RS_NoIllum_color_per_vertex,
										 RS_NoIllum_scalar_per_face,
										 RS_NoIllum_color_per_face,
										 RS_Flat,
										 RS_Flat_scalar_per_vertex,
										 RS_Flat_color_per_vertex,
										 RS_Flat_scalar_per_face,
										 RS_Flat_color_per_face,
										 RS_Phong,
										 RS_Phong_scalar_per_vertex,
										 RS_Phong_color_per_vertex,
										 RS_Phong_scalar_per_face,
										 RS_Phong_color_per_face};
		case VBO_Normal:
			return ILSurfaceRenderStyles{RS_Phong, RS_Phong_scalar_per_vertex, RS_Phong_color_per_vertex,
										 RS_Phong_scalar_per_face, RS_Phong_color_per_face};
		case VBO_Scalar_per_vertex:
			return ILSurfaceRenderStyles{RS_NoIllum_scalar_per_vertex, RS_Flat_scalar_per_vertex,
										 RS_Phong_scalar_per_vertex};
		case VBO_Vec3_per_vertex:
			return ILSurfaceRenderStyles{RS_NoIllum_color_per_vertex, RS_Flat_color_per_vertex,
										 RS_Phong_color_per_vertex};
		case VBO_Scalar_per_face:
			return ILSurfaceRenderStyles{RS_NoIllum_scalar_per_face, RS_Flat_scalar_per_face, RS_Phong_scalar_per_face};
		case VBO_Vec3_per_face:
			return ILSurfaceRenderStyles{RS_NoIllum_color_per_face, RS_Flat_color_per_face, RS_Phong_color_per_face};
		default:
			return ILSurfaceRenderStyles{};
		}
	}

	struct Parameters
	{
		using Shader_types =
			std::tuple<rendering::ShaderPointSprite, rendering::ShaderBoldLine, rendering::ShaderNoIllum,
					   rendering::ShaderNoIllumScalarPerVertex, rendering::ShaderNoIllumColorPerVertex,
					   rendering::ShaderNoIllumScalarPerFace, rendering::ShaderNoIllumColorPerFace,
					   rendering::ShaderFlat, rendering::ShaderFlatScalarPerVertex, rendering::ShaderFlatColorPerVertex,
					   rendering::ShaderFlatScalarPerFace, rendering::ShaderFlatColorPerFace, rendering::ShaderPhong,
					   rendering::ShaderPhongScalarPerVertex, rendering::ShaderPhongColorPerVertex,
					   rendering::ShaderPhongScalarPerFace, rendering::ShaderPhongColorPerFace>;

		template <int N>
		using type_of_param = typename std::tuple_element_t<N, Shader_types>::Param;

		template <int N>
		using type_of_shader = std::tuple_element_t<N, Shader_types>;

		template <int FIRST, int LAST>
		auto gen_params(const rendering::PossibleParameters& pp) -> std::enable_if_t<(FIRST >= LAST)>
		{
			params_[LAST] = type_of_shader<LAST>::generate_param();
			params_[LAST]->pick_parameters(pp);
		}

		template <int FIRST, int LAST>
		auto gen_params(const rendering::PossibleParameters& pp) -> std::enable_if_t<(FIRST < LAST)>
		{
			params_[FIRST] = type_of_shader<FIRST>::generate_param();
			params_[FIRST]->pick_parameters(pp);
			gen_params<FIRST + 1, LAST>(pp);
		}

		inline Parameters()
			: vbo_normal_(nullptr), render_topo_(false), render_vertices_(false), render_edges_(false),
			  render_faces_(true), render_faces_style_(RS_Flat), vertex_scale_factor_(1.0), lw_(5.0f),
			  auto_update_scalar_min_max_(true)
		{
			static rendering::GLColor WHITE{1, 1, 1, 1};
			//			static rendering::GLColor BLACK{0, 0, 0, 1};
			static rendering::GLColor BLUE{0, 0.69f, 0.83f, 1};
			static rendering::GLColor GREEN{0, 1, 0.5f, 1};
			static rendering::GLColor GREY1{0.1f, 0.1f, 0.1f, 1};
			//			static rendering::GLColor YELLOW{0.1f, 0.1f, 0.1f, 1};
			static rendering::GLColor ORANGE{1, 0.5f, 0, 1};

			rendering::PossibleParameters pp = {ORANGE,							  // color
												GREY1,							  // ambiant
												BLUE,							  // front
												GREEN,							  // back
												WHITE,							  // spec
												250.0f,							  // spec coef
												rendering::GLVec3(10, 100, 1000), // Light position
												true,							  // double side
												5.0f,							  // width
												1.0f,							  // size
												0.9f,
												true,
												rendering::GLVec4(0, 0, 0, 0),
												rendering::GLVec4(0, 0, 0, 0)};

			gen_params<RS_Points, RS_Points>(pp);
			pp.color_ = WHITE;
			gen_params<RS_Lines, RS_Lines>(pp);
			pp.color_ = GREEN;
			gen_params<RS_NoIllum, RS_Phong_color_per_face>(pp);

			for (auto& v : vbo_cache_)
				v = nullptr;

			topo_drawer_ = std::make_unique<rendering::TopoDrawer>();
			topo_renderer_ = topo_drawer_->generate_renderer();
		}

		inline void update_param_using(VBOContent v)
		{
			switch (v)
			{
			case VBO_Position:
				for (auto st : style_with_vbo(v))
					params_[st]->set_vbo(1u, vbo_cache_[v]);
				break;
			case VBO_Normal:
				for (auto st : style_with_vbo(v))
					params_[st]->set_vbo(2u, vbo_cache_[v]);
				break;
			case VBO_Scalar_per_vertex:
			case VBO_Vec3_per_vertex:
			case VBO_Scalar_per_face:
			case VBO_Vec3_per_face:
				for (auto st : style_with_vbo(v))
				{
					GLuint attrib = (st < RS_Phong) ? 2u : 3u;
					params_[st]->set_vbo(attrib, vbo_cache_[v]);
				}
				break;
			default:
				break;
			}
		}

		template <int32 RS>
		type_of_param<RS>& param_typed()
		{
			return *static_cast<type_of_param<RS>*>(params_[RS].get());
		}

		inline void update_topo(const MESH& m)
		{
			if (render_topo_ && vertex_position_)
				topo_drawer_->update2D(m, vertex_position_.get());
			topo_dirty_ = false;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec3>> vertex_normal_;
		std::shared_ptr<Attribute<Scalar>> vertex_scalar_;
		std::shared_ptr<Attribute<Scalar>> face_scalar_;
		std::shared_ptr<Attribute<Vec3>> vertex_color_;
		std::shared_ptr<Attribute<Vec3>> face_color_;
		std::unique_ptr<rendering::VBO> vbo_normal_;
		std::array<rendering::VBO*, 6> vbo_cache_;
		std::array<std::unique_ptr<rendering::ShaderParam>, 18> params_;
		/// topo
		std::unique_ptr<rendering::TopoDrawer> topo_drawer_;
		std::unique_ptr<rendering::TopoDrawer::Renderer> topo_renderer_;

		bool render_topo_;
		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool render_faces_smooth_;
		SurfaceRenderStyle render_faces_style_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 lw_;

		bool auto_update_scalar_min_max_;

		bool topo_dirty_;
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
					{
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 7);
						p.topo_dirty_ = true;
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
							p.topo_dirty_ = true;
						}
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
			//			p.update_topo(m);
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, vertex_position.get()) / 7);
			p.vbo_cache_[VBO_Position] = md->update_vbo(vertex_position.get(), true);
		}
		p.update_param_using(VBO_Position);
		p.topo_dirty_ = true;

		v.request_update();
	}

	void set_vertex_normal(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_normal)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);

		p.vertex_normal_ = vertex_normal;
		p.vbo_cache_[VBO_Normal] = md->update_vbo(vertex_normal.get(), true);
		if ((p.vbo_cache_[VBO_Normal] == nullptr) && (p.vbo_cache_[VBO_Position] != nullptr))
		{
			if (p.vbo_normal_ == nullptr)
				p.vbo_normal_ = std::make_unique<rendering::VBO>();
			p.vbo_cache_[VBO_Normal] = p.vbo_normal_.get();
			p.vbo_cache_[VBO_Normal]->bind();
			p.vbo_cache_[VBO_Normal]->allocate(p.vbo_cache_[VBO_Position]->size(), 3);
			p.vbo_cache_[VBO_Normal]->release();
			rendering::MeshRender* mr = md->get_render();
			if (!mr->is_primitive_uptodate(rendering::TRIANGLES))
				mr->init_primitives(m, rendering::TRIANGLES, p.vertex_position_.get());
			compute_normal_engine->compute(p.vbo_cache_[VBO_Position], md->get_render(), p.vbo_cache_[VBO_Normal]);
		}
		std::cout << *(p.vbo_cache_[VBO_Normal]) << std::endl;
		p.update_param_using(VBO_Normal);
		v.request_update();
	}

	void set_vertex_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& vertex_scal)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.vertex_scalar_ = vertex_scal;
		p.vbo_cache_[VBO_Scalar_per_vertex] = md->update_vbo(vertex_scal.get(), true);
		p.update_param_using(VBO_Scalar_per_vertex);
		v.request_update();
	}

	void set_vertex_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_col)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.vertex_color_ = vertex_col;
		p.vbo_cache_[VBO_Vec3_per_vertex] = md->update_vbo(vertex_col.get(), true);
		p.update_param_using(VBO_Vec3_per_vertex);
		v.request_update();
	}

	void set_face_scalar(View& v, const MESH& m, const std::shared_ptr<Attribute<Scalar>>& face_scal)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.face_scalar_ = face_scal;
		p.vbo_cache_[VBO_Scalar_per_face] = md->update_vbo(face_scal.get(), true);
		p.update_param_using(VBO_Scalar_per_face);
		v.request_update();
	}

	void set_face_color(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& face_col)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.face_color_ = face_col;
		p.vbo_cache_[VBO_Vec3_per_face] = md->update_vbo(face_col.get(), true);
		p.update_param_using(VBO_Vec3_per_face);
		v.request_update();
	}

protected:
	void update_scalar_min_max_values(Parameters& p)
	{
		unused_parameters(p);
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

		compute_normal_engine = rendering::ComputeNormalEngine::generate();
	}

	void draw(View* view) override
	{
		auto& ps = parameters_[view];
		for (auto& [m, p] : ps)
		//		for (auto& c : parameters_[view])
		{
			//			MESH* m = c.first;
			//			Parameters& p = c.second;

			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.render_faces_)
			{
				SurfaceRenderStyle rs = p.render_faces_style_;
				//				rendering::ShaderParam* param = p.params_[rs].get();
				if (p.params_[rs]->vao_initialized())
				{
					glEnable(GL_POLYGON_OFFSET_FILL);
					glPolygonOffset(p.lw_ / 5.0f, p.lw_ / 3.0f);
					p.params_[rs]->bind(proj_matrix, view_matrix);

					if (is_per_face(rs))
						md->draw(rendering::TRIANGLES_TB, p.vertex_position_);
					else
						md->draw(rendering::TRIANGLES, p.vertex_position_);

					glDisable(GL_POLYGON_OFFSET_FILL);
					p.params_[rs]->release();
				}
			}
			//			rendering::ShaderParam* param = p.params_[RS_Points].get();
			if (p.render_vertices_ && p.params_[RS_Points]->vao_initialized())
			{
				p.template param_typed<RS_Points>().size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.params_[RS_Points]->bind(proj_matrix, view_matrix);
				md->draw(rendering::POINTS);
				p.params_[RS_Points]->release();
			}

			//			param = p.params_[RS_Lines];
			if (p.render_edges_ && p.params_[RS_Lines]->vao_initialized())
			{
				p.params_[RS_Lines]->bind(proj_matrix, view_matrix);
				//				glEnable(GL_BLEND);
				//				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				md->draw(rendering::LINES);
				//				glDisable(GL_BLEND);
				p.params_[RS_Lines]->release();
			}

			if (p.render_topo_)
			{
				if (p.topo_dirty_)
					p.update_topo(*m);
				p.topo_renderer_->draw(proj_matrix, view_matrix);
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		//		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		//		ImGui::SetWindowSize({0, 0});

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

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const decltype(p.vertex_position_)& att) {
													set_vertex_position(*selected_view_, *selected_mesh_, att);
												});

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_normal_, "Normal",
												[&](const decltype(p.vertex_normal_)& att) {
													set_vertex_normal(*selected_view_, *selected_mesh_, att);
												});

			imgui_combo_attribute<Vertex, Scalar>(*selected_mesh_, p.vertex_scalar_, "Vertex Scalar",
												  [&](const decltype(p.vertex_scalar_)& att) {
													  set_vertex_scalar(*selected_view_, *selected_mesh_, att);
												  });

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, p.vertex_color_, "Vertex Color",
				[&](const decltype(p.vertex_color_)& att) { set_vertex_color(*selected_view_, *selected_mesh_, att); });

			imgui_combo_attribute<Face, Scalar>(
				*selected_mesh_, p.face_scalar_, "Face Scalar",
				[&](const decltype(p.face_scalar_)& att) { set_face_scalar(*selected_view_, *selected_mesh_, att); });

			imgui_combo_attribute<Face, Vec3>(
				*selected_mesh_, p.face_color_, "Face Color",
				[&](const decltype(p.face_color_)& att) { set_face_color(*selected_view_, *selected_mesh_, att); });

			ImGui::Separator();
			need_update |= ImGui::Checkbox("Vertices", &p.render_vertices_);
			need_update |= ImGui::Checkbox("Edges", &p.render_edges_);

			ImGui::BeginGroup();
			ImGui::TextUnformatted("Faces");

			need_update |= ImGui::Checkbox("Faces", &p.render_faces_);
			if (p.render_faces_)
			{
				int32* ptr_style = reinterpret_cast<int32*>(&p.render_faces_style_);

				ImGui::TextUnformatted("NoIllum:");
				need_update |= ImGui::RadioButton("Simple##NoIllum", ptr_style, RS_NoIllum);
				need_update |= ImGui::RadioButton("Scalar/Vertex##NoIllum", ptr_style, RS_NoIllum_scalar_per_vertex);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color/Vertex##NoIllum", ptr_style, RS_NoIllum_color_per_vertex);
				need_update |= ImGui::RadioButton("Scalar/Face##NoIllum", ptr_style, RS_NoIllum_scalar_per_face);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color/Face##NoIllum", ptr_style, RS_NoIllum_color_per_face);
				ImGui::TextUnformatted("Flat:");
				need_update |= ImGui::RadioButton("Simple##Flat", ptr_style, RS_Flat);
				need_update |= ImGui::RadioButton("Scalar/Vertex##Flat", ptr_style, RS_Flat_scalar_per_vertex);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color/Vertex##Flat", ptr_style, RS_Flat_color_per_vertex);
				need_update |= ImGui::RadioButton("Scalar/Face##Flat", ptr_style, RS_Flat_scalar_per_face);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color/Face##Flat", ptr_style, RS_Flat_color_per_face);
				ImGui::TextUnformatted("Phong:");
				need_update |= ImGui::RadioButton("Simple##Phong", ptr_style, RS_Phong);
				need_update |= ImGui::RadioButton("Scalar/Vertex##Phong", ptr_style, RS_Phong_scalar_per_vertex);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color/Vertex##Phong", ptr_style, RS_Phong_color_per_vertex);
				need_update |= ImGui::RadioButton("Scalar/Face##Phong", ptr_style, RS_Phong_scalar_per_face);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Color/Face##Phong", ptr_style, RS_Phong_color_per_face);
				ImGui::EndGroup();

				switch (p.render_faces_style_)
				{
				case RS_NoIllum: {
					auto& param = p.template param_typed<RS_NoIllum>();
					ImGui::Separator();
					ImGui::TextUnformatted("NoIllum parameters");
					need_update |=
						ImGui::ColorEdit3("color##no_illum", param.color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::Checkbox("double side##no_illum", &(param.double_side_));
				}
				break;
				case RS_NoIllum_scalar_per_vertex: {
					auto& param = p.template param_typed<RS_NoIllum_scalar_per_vertex>();
					need_update |= ImGui::Checkbox("double side##no_illum", &(param.double_side_));
					need_update |= imgui_colormap_interface(param.cm_, "no_illum_spv");
					break;
				}
				case RS_NoIllum_color_per_vertex:

					need_update |= ImGui::Checkbox(
						"double side##no_illum", &(p.template param_typed<RS_NoIllum_color_per_vertex>().double_side_));
					break;

				case RS_NoIllum_scalar_per_face: {
					auto& param = p.template param_typed<RS_NoIllum_scalar_per_face>();
					need_update |= ImGui::Checkbox("double side##no_illum", &(param.double_side_));
					need_update |= imgui_colormap_interface(param.cm_, "no_illum_spf");
					break;
				}
				case RS_NoIllum_color_per_face:
					need_update |= ImGui::Checkbox("double side##no_illum",
												   &(p.template param_typed<RS_NoIllum_color_per_face>().double_side_));
					break;

				case RS_Flat: {
					auto& param = p.template param_typed<RS_Flat>();
					ImGui::Separator();
					ImGui::TextUnformatted("Flat parameters");
					need_update |=
						ImGui::ColorEdit3("front color##flat", param.front_color_.data(), ImGuiColorEditFlags_NoInputs);
					if (param.double_side_)
						need_update |= ImGui::ColorEdit3("back color##flat", param.back_color_.data(),
														 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::Checkbox("double side##flat", &(param.double_side_));
				}
				break;
				case RS_Flat_scalar_per_vertex: {
					auto& param = p.template param_typed<RS_Flat_scalar_per_vertex>();
					need_update |= ImGui::Checkbox("double side##flat", &(param.double_side_));
					need_update |= imgui_colormap_interface(param.cm_, "flat_spv");
					break;
				}
				case RS_Flat_color_per_vertex:
					need_update |= ImGui::Checkbox("double side##flat",
												   &(p.template param_typed<RS_Flat_color_per_vertex>().double_side_));
					break;

				case RS_Flat_scalar_per_face: {
					auto& param = p.template param_typed<RS_Flat_scalar_per_face>();
					need_update |= ImGui::Checkbox("double side##flat", &(param.double_side_));
					need_update |= imgui_colormap_interface(param.cm_, "flat_spf");
					break;
				}
				case RS_Flat_color_per_face:
					need_update |= ImGui::Checkbox("double side##flat",
												   &(p.template param_typed<RS_Flat_color_per_face>().double_side_));
					break;

				case RS_Phong: {
					auto& param = p.template param_typed<RS_Phong>();
					ImGui::Separator();
					ImGui::TextUnformatted("Phong parameters");
					need_update |= ImGui::ColorEdit3("front color##phong", param.front_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					if (param.double_side_)
						need_update |= ImGui::ColorEdit3("back color##phong", param.back_color_.data(),
														 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("spec##phong", &(param.specular_coef_), 10.0f, 1000.0f);
					need_update |= ImGui::Checkbox("double side##phong", &(param.double_side_));
					break;
				}

				case RS_Phong_color_per_face: {
					auto& param = p.template param_typed<RS_Phong_color_per_face>();
					ImGui::Separator();
					ImGui::TextUnformatted("Phong parameters");
					need_update |= ImGui::SliderFloat("spec##phong", &(param.specular_coef_), 10.0f, 1000.0f);
					need_update |= ImGui::Checkbox("double side##phong", &(param.double_side_));
					break;
				}

				case RS_Phong_scalar_per_face: {
					auto& param = p.template param_typed<RS_Phong_scalar_per_face>();
					ImGui::Separator();
					ImGui::TextUnformatted("Phong parameters");
					need_update |= ImGui::SliderFloat("spec##phong", &(param.specular_coef_), 10.0f, 1000.0f);
					need_update |= ImGui::Checkbox("double side##phong", &(param.double_side_));
					need_update |= imgui_colormap_interface(param.cm_, "phong_spf");
					break;
				}
				case RS_Phong_color_per_vertex: {
					auto& param = p.template param_typed<RS_Phong_color_per_vertex>();
					ImGui::Separator();
					ImGui::TextUnformatted("Phong parameters");
					need_update |= ImGui::SliderFloat("spec##phong", &(param.specular_coef_), 10.0f, 1000.0f);
					need_update |= ImGui::Checkbox("double side##phong", &(param.double_side_));
					break;
				}

				case RS_Phong_scalar_per_vertex: {
					auto& param = p.template param_typed<RS_Phong_scalar_per_vertex>();
					ImGui::Separator();
					ImGui::TextUnformatted("Phong parameters");
					need_update |= ImGui::SliderFloat("spec##phong", &(param.specular_coef_), 10.0f, 1000.0f);
					need_update |= ImGui::Checkbox("double side##phong", &(param.double_side_));
					need_update |= imgui_colormap_interface(param.cm_, "phong_spf");
					break;
				}

				default:
					break;
				}
			}

			if (p.render_edges_)
			{
				auto& param = p.template param_typed<RS_Lines>();
				ImGui::Separator();
				ImGui::TextUnformatted("Edges parameters");
				need_update |= ImGui::ColorEdit3("color##edges", param.color_.data(), ImGuiColorEditFlags_NoInputs);
				if (ImGui::SliderFloat("width##edges", &(param.width_), 3.0f, 15.0f))
				{
					p.lw_ = param.width_;
					need_update = true;
				}
			}

			if (p.render_vertices_)
			{
				auto& param = p.template param_typed<RS_Points>();
				ImGui::Separator();
				ImGui::TextUnformatted("Vertices parameters");
				need_update |= ImGui::ColorEdit3("color##vertices", param.color_.data(), ImGuiColorEditFlags_NoInputs);
				need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
			}

			ImGui::Separator();
			if (ImGui::Checkbox("Topo", &p.render_topo_))
			{
				need_update = true;
				if (p.render_topo_)
				{
					ImGui::Separator();
					ImGui::TextUnformatted("Topo parameters");
					need_update |= ImGui::ColorEdit3("colorDarts", p.topo_drawer_->dart_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::ColorEdit3("colorPhi2", p.topo_drawer_->phi2_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::ColorEdit3("colorPhi3", p.topo_drawer_->phi3_color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					if (ImGui::SliderFloat("explodeEdges", &(p.topo_drawer_->shrink_e_), 0.01f, 1.0f))
						p.topo_dirty_ = true;

					if (ImGui::SliderFloat("explodeFaces", &(p.topo_drawer_->shrink_f_), 0.01f, 1.0f))
						p.topo_dirty_ = true;
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
	rendering::ComputeNormalEngine* compute_normal_engine;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
