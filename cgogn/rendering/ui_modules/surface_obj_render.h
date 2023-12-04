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

#ifndef CGOGN_MODULE_SURFACE_OBJ_RENDER_H_
#define CGOGN_MODULE_SURFACE_OBJ_RENDER_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>
#include <cgogn/geometry/types/vector_traits.h>


#include <cgogn/rendering/shaders/shader_obj_flat_texture.h>
#include <cgogn/rendering/shaders/shader_obj_meshuv.h>
#include <cgogn/rendering/shaders/shader_mesh_2d_edges.h>
#include <cgogn/rendering/texture.h>

#include <boost/synapse/connect.hpp>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec2;

template <typename MESH>
class SurfaceObjRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceRender can only be used with meshes of dimension >= 2");


	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	// using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_position_vbo_(nullptr), vertex_tc_(nullptr), vertex_tc_vbo_(nullptr),
			  draw_flatten_(false), draw_bound_pos_(false), draw_bound_tc_(false)
		{
			param_textured_ = rendering::ShaderObjFlatTexture::generate_param();
			param_flatten_ = rendering::ShaderObjMeshUV::generate_param();
			param_boundary_edges_ = rendering::ShaderMesh2DEdges::generate_param();
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::unique_ptr<rendering::ShaderObjMeshUV::Param> param_flatten_;
		std::unique_ptr<rendering::ShaderObjFlatTexture::Param> param_textured_;
		std::unique_ptr<rendering::ShaderMesh2DEdges::Param> param_boundary_edges_;		

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vertex_position_vbo_;

		std::shared_ptr<Attribute<Vec2>> vertex_tc_;
		rendering::VBO* vertex_tc_vbo_;

		std::array<rendering::EBO, 2> tri_ebos_;
		bool draw_flatten_;
		std::array<rendering::EBO, 3> boundary_ebos_;
		bool draw_bound_pos_;
		bool draw_bound_tc_;
	};

public:
	SurfaceObjRender(const App& app)
		: ViewModule(app, "SurfaceRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_pos_(nullptr)
	{
	}

	~SurfaceObjRender()
	{
	}

private:

	void update_ebo(Parameters& p, const MESH* m)
	{
		//const MESH* mtc = tc_mesh(m);
		//const MESH* mno = no_mesh(m);
		auto pm = attached_meshes_[m];

		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Face = typename mesh_traits<MESH>::Face;
		using Edge = typename mesh_traits<MESH>::Edge;
		std::vector<uint32> table_pos_indices;
		table_pos_indices.reserve(8192);
		std::vector<uint32> table_tc_indices;
		table_tc_indices.reserve(8192);
		std::vector<Vertex> vertices;
		vertices.reserve(32u);
		foreach_cell(*m, [&](Face f) -> bool {
			vertices.clear();
			append_incident_vertices(*m, f, vertices);
			for (uint32 i = 1; i < uint32(vertices.size()) - 1; ++i)
			{
				table_pos_indices.push_back(index_of(*m, vertices[0]));
				table_pos_indices.push_back(index_of(*m, vertices[i]));
				table_pos_indices.push_back(index_of(*m, vertices[i + 1]));
				table_tc_indices.push_back(index_of(*pm.first, vertices[0]));
				table_tc_indices.push_back(index_of(*pm.first, vertices[i]));
				table_tc_indices.push_back(index_of(*pm.first, vertices[i + 1]));
			}
			return true;
		});

		if (!p.tri_ebos_[0].is_created())
			p.tri_ebos_[0].create();
		p.tri_ebos_[0].bind();
		p.tri_ebos_[0].allocate(table_pos_indices.size());
		uint32* ptr =p.tri_ebos_[0].lock_pointer();
		std::memcpy(ptr, table_pos_indices.data(), sizeof(uint32) * table_pos_indices.size());
		p.tri_ebos_[0].set_name("EBO_Pos");
		p.tri_ebos_[0].release_pointer();
		p.tri_ebos_[0].release();

		if (!p.tri_ebos_[1].is_created())
			p.tri_ebos_[1].create();
		p.tri_ebos_[1].bind();
		p.tri_ebos_[1].allocate(table_pos_indices.size());
		ptr = p.tri_ebos_[1].lock_pointer();
		std::memcpy(ptr, table_tc_indices.data(), sizeof(uint32) * table_pos_indices.size());
		p.tri_ebos_[1].set_name("EBO_TC");
		p.tri_ebos_[1].release_pointer();
		p.tri_ebos_[1].release();

		if (!p.tri_ebos_[1].is_created())
			p.tri_ebos_[1].create();
		p.tri_ebos_[1].bind();
		p.tri_ebos_[1].allocate(table_pos_indices.size());
		ptr = p.tri_ebos_[1].lock_pointer();
		std::memcpy(ptr, table_tc_indices.data(), sizeof(uint32) * table_pos_indices.size());
		p.tri_ebos_[1].set_name("EBO_TC");
		p.tri_ebos_[1].release_pointer();
		p.tri_ebos_[1].release();


		std::vector<uint32> table_pos_boundary_indices;
		std::vector<uint32> table_tc_boundary_indices;
		std::vector<uint32> table_no_boundary_indices;
		table_pos_boundary_indices.reserve(8192);
		table_tc_boundary_indices.reserve(8192);
		table_no_boundary_indices.reserve(8192);


		//foreach_cell(*m, [&](Edge e) -> bool {
		//	vertices.clear();
		//	append_incident_vertices(*m, e, vertices);
		//		table_pos_boundary_indices.push_back(index_of(*m, vertices[0]));
		//		table_pos_boundary_indices.push_back(index_of(*m, vertices[1]));
		//	return true;
		//});
	
		const MESH& map_pos = *m;
		const MESH& map_tc = *pm.first;
		const MESH& map_no = *pm.second;

		foreach_cell(map_pos, [&](Edge e) -> bool {
			if (is_boundary(map_pos, e.dart_) || is_boundary(map_pos, phi2(map_pos, e.dart_)))
			{
				table_pos_boundary_indices.push_back(index_of(map_tc, Vertex(e.dart_)));
				table_pos_boundary_indices.push_back(index_of(map_tc, Vertex(phi1(map_pos, e.dart_))));
			}
			return true;
		});

		foreach_cell(map_tc, [&](Edge e) -> bool {
			if (is_boundary(map_tc, e.dart_) || is_boundary(map_tc, phi2(map_tc, e.dart_)))
			{
				table_tc_boundary_indices.push_back(index_of(map_tc, Vertex(e.dart_)));
				table_tc_boundary_indices.push_back(index_of(map_tc, Vertex(phi1(map_tc, e.dart_))));
			}
			return true;
		});

		foreach_cell(map_no, [&](Edge e) -> bool {
			if (is_boundary(map_no, e.dart_) || is_boundary(map_no, phi2(map_no, e.dart_)))
			{
				table_no_boundary_indices.push_back(index_of(map_tc, Vertex(e.dart_)));
				table_no_boundary_indices.push_back(index_of(map_tc, Vertex(phi1(map_no, e.dart_))));
			}
			return true;
		});

		if (!p.boundary_ebos_[0].is_created())
			p.boundary_ebos_[0].create();
		p.boundary_ebos_[0].bind();
		p.boundary_ebos_[0].allocate(table_pos_boundary_indices.size());
		ptr = p.boundary_ebos_[0].lock_pointer();
		std::memcpy(ptr, table_pos_boundary_indices.data(), sizeof(uint32) * table_pos_boundary_indices.size());
		p.boundary_ebos_[0].set_name("EBO_Bound_pos");
		p.boundary_ebos_[0].release_pointer();
		p.boundary_ebos_[0].release();

		if (!p.boundary_ebos_[1].is_created())
			p.boundary_ebos_[1].create();
		p.boundary_ebos_[1].bind();
		p.boundary_ebos_[1].allocate(table_tc_boundary_indices.size());
		ptr = p.boundary_ebos_[1].lock_pointer();
		std::memcpy(ptr, table_tc_boundary_indices.data(), sizeof(uint32) * table_tc_boundary_indices.size());
		p.boundary_ebos_[1].set_name("EBO_TC_pos");
		p.boundary_ebos_[1].release_pointer();
		p.boundary_ebos_[1].release();

		if (!p.boundary_ebos_[2].is_created())
			p.boundary_ebos_[2].create();
		p.boundary_ebos_[2].bind();
		p.boundary_ebos_[2].allocate(table_no_boundary_indices.size());
		ptr = p.boundary_ebos_[2].lock_pointer();
		std::memcpy(ptr, table_no_boundary_indices.data(), sizeof(uint32) * table_no_boundary_indices.size());
		p.boundary_ebos_[2].set_name("EBO_No_pos");
		p.boundary_ebos_[2].release_pointer();
		p.boundary_ebos_[2].release();
	}

	void init_mesh(MESH* m)
	{
		const std::string& bname = mesh_provider_->mesh_name(*m);
		std::string tc_name = bname + "_tc";
		std::string no_name = bname + "_no";
		std::pair<MESH*, MESH*> att = std::make_pair<MESH*, MESH*>(nullptr, nullptr);
		mesh_provider_->foreach_mesh([&](MESH& mm, const std::string& name) {
			if (name == tc_name)
				att.first = &mm;
			if (name == no_name)
				att.second = &mm;
		});
		attached_meshes_[m] = att;
		

		for (View* v : linked_views_)
		{
			Parameters& p = parameters_[v][m];
			update_ebo(p, m);
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);
		}
	}



public:
	void load_texture(const std::string& img_name)
	{
		rendering::GLImage img(img_name);
		tex_->load(img);
	}

	void load_texture(const rendering::GLImage& img)
	{
		tex_->load(img);
	}

	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_position_ == vertex_position)
			return;

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.vertex_position_vbo_ = md.update_vbo(p.vertex_position_.get(), true);

			auto pm = attached_meshes_[&m];
			MeshData<MESH>& md_tc = mesh_provider_->mesh_data(*pm.first);
			p.vertex_tc_ = cgogn::get_attribute<Vec2, Vertex>(*md_tc.mesh_, "position");
			p.vertex_tc_vbo_ = md_tc.update_vbo(p.vertex_tc_.get(), true);
		}
		else
		{
			p.vertex_position_vbo_ = nullptr;
			p.vertex_tc_vbo_ = nullptr;
		}

		p.param_textured_->set_vbos({p.vertex_position_vbo_, p.vertex_tc_vbo_});
		p.param_flatten_->set_vbos({p.vertex_tc_vbo_});
		p.param_boundary_edges_->set_vbos({p.vertex_tc_vbo_});
		v.request_update();	
	}

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SurfaceObjRender<MESH>::init_mesh));
		tex_ = std::make_shared<rendering::Texture2D>(
			std::vector<std::pair<GLenum, GLint>>{
			{GL_TEXTURE_MIN_FILTER, GL_LINEAR}, {GL_TEXTURE_MAG_FILTER, GL_NEAREST}, {GL_TEXTURE_WRAP_S, GL_REPEAT},
													   {GL_TEXTURE_WRAP_T, GL_REPEAT}});
	}

	void draw(View* view) override
	{
		for (auto& [m, pp] : parameters_[view])
		{
			Parameters& p = pp;
			MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.draw_flatten_)
			{
				rendering::GLVec2 ratio = (view->viewport_width() > view->viewport_height())
						? rendering::GLVec2(float32(view->viewport_height()) / view->viewport_width(), 1.0f)
						: rendering::GLVec2(1.0f, view->viewport_width() / float32(view->viewport_height()));
				ratio *= 0.95f;
				if (p.param_flatten_->attributes_initialized())
				{
					p.param_flatten_->ratio_ = ratio;
					glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
					p.tri_ebos_[1].bind_texture_buffer(10);
					p.param_flatten_->bind();
					glDrawArraysInstanced(GL_TRIANGLES, 0, 3, p.tri_ebos_[0].size() / 3);
					p.param_flatten_->release();
					p.tri_ebos_[1].release_texture_buffer(10);
					glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
				}

				if (p.draw_bound_pos_ && (p.param_boundary_edges_->attributes_initialized()))
				{
					glDisable(GL_DEPTH_TEST);
					p.boundary_ebos_[0].bind_texture_buffer(10);
					p.param_boundary_edges_->ratio_ = ratio;
					p.param_boundary_edges_->color_ = rendering::GLColor(1, 0, 0, 1);
					p.param_boundary_edges_->bind();
					glDrawArrays(GL_LINES, 0, p.boundary_ebos_[0].size());
					p.param_boundary_edges_->release();
					p.boundary_ebos_[0].release_texture_buffer(10);
					glDisable(GL_DEPTH_TEST);
				}

				if (p.draw_bound_tc_ && (p.param_boundary_edges_->attributes_initialized()))
				{
					glDisable(GL_DEPTH_TEST);
					p.boundary_ebos_[1].bind_texture_buffer(10);
					p.param_boundary_edges_->ratio_ = ratio;
					p.param_boundary_edges_->color_ = rendering::GLColor(0, 1, 1, 1);
					p.param_boundary_edges_->bind();
					glDrawArrays(GL_LINES, 0, p.boundary_ebos_[1].size());
					p.param_boundary_edges_->release();
					p.boundary_ebos_[1].release_texture_buffer(10);
					glDisable(GL_DEPTH_TEST);
				}
			}
			else
			{
				if (p.param_textured_->attributes_initialized())
				{
					p.param_textured_->texture_ = tex_;
					p.tri_ebos_[0].bind_texture_buffer(10);
					p.tri_ebos_[1].bind_texture_buffer(11);
					p.param_textured_->bind(proj_matrix, view_matrix);
					glDrawArraysInstanced(GL_TRIANGLES, 0, 3, p.tri_ebos_[0].size() / 3);
					p.param_textured_->release();
					p.tri_ebos_[1].release_texture_buffer(11);
					p.tri_ebos_[0].release_texture_buffer(10);
				}
			}
		}
	}



	void left_panel() override
	{
		bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_pos_, "Surf Obj Position", [&](MESH& m) {
			selected_mesh_pos_ = &m;
			const std::string& bname = mesh_provider_->mesh_name(m);
			std::string tc_name = bname + "_tc";
//			std::string no_name = bname + "_no";
			mesh_provider_->foreach_mesh([&](MESH& mm, const std::string& name) {
				if (name == tc_name)
					selected_mesh_tc_ = &mm;
//				if (name == no_name)
//					selected_mesh_no_ = &mm;
			});

		});

		if (selected_view_ && selected_mesh_pos_ && selected_mesh_tc_ /* && selected_mesh_no_*/)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_pos_];
			if (ImGui::Checkbox("draw_flatten", &p.draw_flatten_))
				need_update = true;
			if (p.draw_flatten_)
			{
				if (ImGui::Checkbox("draw_boundary position", &p.draw_bound_pos_))
					need_update = true;
				if (ImGui::Checkbox("draw_boundary tex coord", &p.draw_bound_tc_))
					need_update = true;
			}
			if (ImGui::Checkbox("draw_param", &p.param_textured_->draw_param_))
				need_update = true;

			if (need_update)
				for (View* v : linked_views_)
					v->request_update();
		}
	}

private:
	View* selected_view_;
	const MESH* selected_mesh_pos_;
	const MESH* selected_mesh_tc_;
//	const MESH* selected_mesh_no_;
	std::unordered_map<View*, std::unordered_map<const MESH*, Parameters>> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
	std::shared_ptr<rendering::Texture2D> tex_;
	std::unordered_map<const MESH*, std::pair<MESH*,MESH*>> attached_meshes_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
