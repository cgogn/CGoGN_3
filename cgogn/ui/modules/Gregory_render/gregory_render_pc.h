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

#ifndef CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_
#define CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/drawer.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class GregoryRenderPC : public ViewModule
{

	static_assert(mesh_traits<MESH>::dimension >= 2,
				  "GregoryPC can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> renderer_ = nullptr;
		bool show_ = true;
		Parameters() = default;
	};

public:
	GregoryRenderPC(const App& app)
		: ViewModule(app, "GregoryPC (" + std::string{mesh_traits<MESH>::name} + ")"),
		  initialized_(false),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~GregoryRenderPC()
	{
	}

private:
	void init_pc(const MESH* m)
	{
		auto vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
		auto& dr = drawers_[m];
		dr = std::make_unique<cgogn::rendering::DisplayListDrawer>();
		dr->new_list();

		dr->ball_size(1.0);
 		dr->begin(GL_POINTS); // or GL_POINTS, GL_LINES, GL_TRIANGLES
 		dr->color3f(1.0,0.0,0.0);
		foreach_cell(*m,[&](Vertex v) -> bool
		{
			const Vec3& P = value<Vec3>(*m,vertex_position,v);
			dr->vertex3f(P[0],P[1],P[2]);
			return true;
		});
		dr->end();
		dr->line_width(2.0);
		dr->begin(GL_LINES); // or GL_POINTS, GL_LINES, GL_TRIANGLES
 		dr->color3f(0.0,1.0,0.0);
		foreach_cell(*m,[&](Edge e) -> bool
		{
			const Vec3& P = value<Vec3>(*m,vertex_position,Vertex(e.dart));
			dr->vertex3f(P[0],P[1],P[2]);
			const MESH& rm = *m;
			const Vec3& Q = value<Vec3>(*m, vertex_position, Vertex(phi1(rm, e.dart)));
			dr->vertex3f(Q[0],Q[1],Q[2]);
			return true;
		});
		dr->end();

		dr->end_list();
		drawers_[m] = std::move(dr);

		for (View* v : linked_views_)
			parameters_[v][m].renderer_ = dr->generate_renderer();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));

		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) {init_pc(m);});

		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &GregoryRenderPC<MESH>::init_pc));
	}

	void draw(View* view) override
	{
		for (auto&[m,p] : parameters_[view])
		{
			if (p.show_)
			{
				MeshData<MESH>* md = mesh_provider_->mesh_data(m);
				const rendering::GLMat4& proj_matrix = view->projection_matrix();
				const rendering::GLMat4& view_matrix = view->modelview_matrix();
				p.renderer_->draw(proj_matrix, view_matrix);
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });
		imgui_mesh_selector(mesh_provider_, selected_mesh_, [&](MESH* m) {
			selected_mesh_ = m;
			mesh_provider_->mesh_data(m)->outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];
			ImGui::Separator();
			need_update |= ImGui::Checkbox("Show", &p.show_);
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	bool initialized_;
	MeshProvider<MESH>* mesh_provider_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	View* selected_view_;
	const MESH* selected_mesh_;
	std::unordered_map<const MESH*, std::unique_ptr<cgogn::rendering::DisplayListDrawer>> drawers_;
	std::unordered_map<View*, std::unordered_map<const MESH*, Parameters>> parameters_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_VECTOR_H_
