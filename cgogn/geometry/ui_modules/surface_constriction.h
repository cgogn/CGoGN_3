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
		if (!geodesic_path_.empty()) 
		{
			param_edge_->bind(view->projection_matrix(), view->modelview_matrix());
			glDrawArrays(GL_LINES, 0, 2 * geodesic_path_.size());
			param_edge_->release();
		}
	}

	/**
	 * Show geodesic path
	*/
	void update_vbo()
	{
		std::vector<Vec3> geodesic_segments;
		geodesic_segments.reserve(2 * geodesic_path_.size());
		for (Edge e : geodesic_path_)
		{
			foreach_incident_vertex(intr->getMesh(), e, [&](Vertex v) -> bool {
				geodesic_segments.push_back(value<Vec3>(intr->getMesh(), selected_vertex_position_, v));
				return true;
			});
		}
		rendering::update_vbo(geodesic_segments, &edges_vbo_);
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
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			if (selected_vertex_position_ && !intrinsic_made)
			{
				if (ImGui::Button("Build intrinsic mesh"))
				{
					// warning ! only CMap2, not yet templated
					intr = std::make_shared<geometry::IntrinsicTriangulation>(*selected_mesh_, selected_vertex_position_);
					//mesh_provider_->register_mesh(intr->getMesh(), "intrinsic");	// visualization of the topology
					intrinsic_made = true;
				}
			}
			if (intrinsic_made)
			{
				if (ImGui::Button("Display intrinsic topology"))	// TODO remove
				{
					geodesic_path_.clear();
					foreach_cell(intr->getMesh(), [&](Edge e) -> bool {
						geodesic_path_.push_back(e);
						return true;
					});
					update_vbo();
					need_update = true;
				}
				if (ImGui::Button("Init random path"))	//TODO dijsktra / selection
				{
					// "brute force and probably good" path
					std::vector<Edge> tmp_vec;
					foreach_cell(intr->getMesh(), [&](Edge e) -> bool {
						tmp_vec.push_back(e);
						return true;
					});
					Dart a = tmp_vec[rand()%tmp_vec.size()].dart;
					geodesic_path_.clear();
					for (int i = 3; i > 0; --i)
					{
						geodesic_path_.push_back(Edge(a));
						a = phi<1, 2, 1>(intr->getMesh(), a);
					}
					update_vbo();
					need_update = true;
				}
				if (geodesic_path_.size() > 0 && ImGui::Button("Compute geodesic"))
				{
					geometry::geodesic_path(*intr, geodesic_path_);
					update_vbo();
					need_update = true;
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	bool intrinsic_made = false;	// only one intrinsic triangulation
	MESH* selected_mesh_ = nullptr;
	std::shared_ptr<geometry::IntrinsicTriangulation> intr;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;

	std::list<Edge> geodesic_path_;
	rendering::VBO edges_vbo_;
	std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;

	MeshProvider<MESH>* mesh_provider_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_CONSTRICTION_H_
