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
		param_edge_ = rendering::ShaderBoldLine::generate_param();
		param_edge_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
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
		if (edges_vbo_.size() > 0)
		{
			param_edge_->bind(view->projection_matrix(), view->modelview_matrix());
			glDrawArrays(GL_LINES, 0, edges_vbo_.size());
			param_edge_->release();
		}
	}

	/**
	 * Show geodesic path
	*/
	void update_vbo()
	{
		if (show_topology)
			rendering::update_vbo(intr->edge_list_topology(geodesic_path_), &edges_vbo_);
		else
			rendering::update_vbo(intr->edge_list_trace(geodesic_path_), &edges_vbo_);
		
		need_update = true;
	}

	void update_path_from_edge_set()
	{
		geodesic_path_.clear();
		if (selected_edges_set_ == nullptr)
			return;

		selected_edges_set_->foreach_cell([&] (Edge e) -> bool {
			geodesic_path_.push_back(e);
			return true;
		});
	}

	void update_path_from_vertex_set()
	{
		geodesic_path_.clear();
		if (selected_vertices_set_ == nullptr)
			return;
		
		std::vector<Vertex> points_list;
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
				ImGui::TextUnformatted("Selection set cell type");
				ImGui::BeginGroup();
				if (ImGui::RadioButton("Vertex##CellType", enum_selected_set_type == 0))
					enum_selected_set_type = 0;
				ImGui::SameLine();

				if (ImGui::RadioButton("Edge##CellType", enum_selected_set_type == 1))
					enum_selected_set_type = 1;
				ImGui::SameLine();
				ImGui::EndGroup();

				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

				if (enum_selected_set_type == 0)
					imgui_combo_cells_set(md, selected_vertices_set_, "Source vertices",
										  [&](CellsSet<MESH, Vertex>* cs) { selected_vertices_set_ = cs; });
				else if (enum_selected_set_type == 1)
					imgui_combo_cells_set(md, selected_edges_set_, "Source edges",
										  [&](CellsSet<MESH, Edge>* cs) { selected_edges_set_ = cs; });

				if (ImGui::Button("Compute path")) {
					if (enum_selected_set_type == 0)
						update_path_from_vertex_set();
					else if (enum_selected_set_type == 1)
						update_path_from_edge_set();
					intr = std::make_shared<geometry::IntrinsicTriangulation>(*selected_mesh_, selected_vertex_position_);
					update_vbo();
				}

				if (ImGui::Button("Compute geodesic")) {
					if (enum_selected_set_type == 0)
						update_path_from_vertex_set();
					else if (enum_selected_set_type == 1)
						update_path_from_edge_set();
					intr = std::make_shared<geometry::IntrinsicTriangulation>(*selected_mesh_, selected_vertex_position_);
					geometry::geodesic_path(*intr, geodesic_path_, flip_out_iteration);
					update_vbo();
				}

				ImGui::InputInt("Flip out iterations", &flip_out_iteration);
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
	bool show_topology = false;
	int enum_selected_set_type = 0; // 0 -> vertex, 1 -> edges
	int flip_out_iteration = 100;
	MESH* selected_mesh_ = nullptr;
	std::shared_ptr<geometry::IntrinsicTriangulation> intr;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	CellsSet<MESH, Vertex>* selected_vertices_set_ = nullptr;
	CellsSet<MESH, Edge>* selected_edges_set_ = nullptr;

	std::list<Edge> geodesic_path_;
	rendering::VBO edges_vbo_;
	std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;

	MeshProvider<MESH>* mesh_provider_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_CONSTRICTION_H_
