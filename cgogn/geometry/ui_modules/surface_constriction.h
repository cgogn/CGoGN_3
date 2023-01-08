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

#ifndef CGOGN_MODULE_SURFACE_CONSTRICTION_H_
#define CGOGN_MODULE_SURFACE_CONSTRICTION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/vbo_update.h>

#include <GLFW/glfw3.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/geodesic.h>
#include <cgogn/geometry/algos/intrinsic_triangulation.h>
#include <boost/synapse/connect.hpp>
#include <unordered_set>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceConstriction : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2,
				  "SurfaceConstriction can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Volume = typename mesh_traits<MESH>::Volume;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceConstriction(const App& app)
		: ViewModule(app, "SurfaceConstriction (" + std::string{mesh_traits<MESH>::name} + ")")
	{
		param_intr_ = rendering::ShaderBoldLine::generate_param();
		param_intr_->color_ = rendering::GLColor(1, 1, 1, 0.5f);
		param_intr_->width_ = 3.0f;
		param_intr_->set_vbos({&intr_vbo_});

		param_edge_ = rendering::ShaderBoldLine::generate_param();
		param_edge_->color_ = rendering::GLColor(1, 0, 0, 1);
		param_edge_->width_ = 5.0f;
		param_edge_->set_vbos({&edges_vbo_});
	}

	~SurfaceConstriction()
	{
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		selected_vertex_position_ = vertex_position;
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void draw(View* view) override
	{
		if (intr_traced_.size() > 0)
		{
			param_intr_->bind(view->projection_matrix(), view->modelview_matrix());
			glDrawArrays(GL_LINES, 0, intr_vbo_.size());
			param_intr_->release();
		}

		if (edges_vbo_.size() > 0)
		{
			param_edge_->bind(view->projection_matrix(), view->modelview_matrix());
			glDrawArrays(GL_LINES, 0, edges_vbo_.size());
			param_edge_->release();
		}
	}

	/**
	 * Update geodesic path VBO
	*/
	void update_vbo()
	{
		rendering::update_vbo(intr->edge_list_trace(geodesic_path_), &edges_vbo_);
		need_update = true;
	}

	/**
	 * Show complete intrinsic triangulation
	*/
	void update_intr_traced_set()
	{
		intr_traced_.clear();
		if (show_intr && intr)
		{
			foreach_cell(intr->getMesh(), [&](Edge e) -> bool {
				std::list<Edge> oneEdgeList;
				oneEdgeList.push_back(e);
				std::vector<Vec3> toappend = intr->edge_list_trace(oneEdgeList);
				intr_traced_.insert(intr_traced_.end(), toappend.begin(), toappend.end());
				return true;
			});
		}
		rendering::update_vbo(intr_traced_, &intr_vbo_);
		need_update = true;
	}

	/**
	 * Construct a path with vertices
	*/
	void update_path_from_vertex_set()
	{
		geodesic_path_.clear();
		if (selected_vertices_set_ == nullptr)
			return;
		
		std::vector<Vertex> points_list;	// could not be ordered, implementation defined
		selected_vertices_set_->foreach_cell([&] (Vertex v) -> bool {
			points_list.push_back(v);
			return true;
		});

		int size = points_list.size();
		for (int i=1; i < size; ++i)
		{
			std::list<Edge> segment = geometry::find_path(*selected_mesh_, points_list[i-1], points_list[i]);
			geodesic_path_.splice(geodesic_path_.end(), segment);
		}
		if (cyclic && size > 0)
		{
			std::list<Edge> segment = geometry::find_path(*selected_mesh_, points_list[size-1], points_list[0]);
			geodesic_path_.splice(geodesic_path_.end(), segment);
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			if (selected_vertex_position_)
			{
				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

				imgui_combo_cells_set(md, selected_vertices_set_, "Source vertices",
										[&](CellsSet<MESH, Vertex>* cs) { selected_vertices_set_ = cs; });

				if (ImGui::Button("Compute path")) {
					update_path_from_vertex_set();
					intr = std::make_shared<geometry::IntrinsicTriangulation>(*selected_mesh_, selected_vertex_position_);
					update_intr_traced_set();
					update_vbo();
				}

				if (ImGui::Button("Compute geodesic") || ImGui::InputInt("Flip out iterations", &flip_out_iteration)) {
					update_path_from_vertex_set();
					intr = std::make_shared<geometry::IntrinsicTriangulation>(*selected_mesh_, selected_vertex_position_);
					geometry::geodesic_path(*intr, geodesic_path_, flip_out_iteration, cyclic);
					update_intr_traced_set();
					update_vbo();
				}

				if (ImGui::Checkbox("Closed loop", &cyclic)) {
					update_path_from_vertex_set();
					intr = std::make_shared<geometry::IntrinsicTriangulation>(*selected_mesh_, selected_vertex_position_);
					update_vbo();
				}

				if (ImGui::Checkbox("Show intrinsic mesh", &show_intr)) {
					update_intr_traced_set();
				}
			}
		}

		if (need_update)
		{
			for (View* v : linked_views_)
				v->request_update();
			need_update = false;
		}
	}

private:
	bool need_update = false;
	bool show_intr = false;
	bool cyclic = false;
	int flip_out_iteration = 100;
	MESH* selected_mesh_ = nullptr;
	std::shared_ptr<geometry::IntrinsicTriangulation> intr;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	CellsSet<MESH, Vertex>* selected_vertices_set_ = nullptr;

	std::list<Edge> geodesic_path_;
	std::vector<Vec3> intr_traced_;
	rendering::VBO edges_vbo_;
	rendering::VBO intr_vbo_;
	std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<rendering::ShaderBoldLine::Param> param_intr_;

	MeshProvider<MESH>* mesh_provider_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_CONSTRICTION_H_
