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

#ifndef CGOGN_MODULE_VOLUME_RENDER_H_
#define CGOGN_MODULE_VOLUME_RENDER_H_

#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/volume_drawer.h>
#include <boost/synapse/connect.hpp>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class VolumeRender : public ViewModule
{
    static_assert(mesh_traits<MESH>::dimension == 3, "VolumeRender can only be used with meshes of dimension 3");

    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Volume = typename mesh_traits<MESH>::Volume;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;
    using GLColor = rendering::GLColor;

    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	struct Parameters
	{
		Parameters() :
			vertex_position_(nullptr),
			render_volumes_(true),
                        render_volumes_edges_(false),
                        explode_factor_(0.9f)
		{
                    vol_drw_ = std::make_unique<cgogn::rendering::VolumeDrawer>();
                    vol_rend_ = vol_drw_->generate_renderer();
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

                std::shared_ptr<Attribute<Vec3>> vertex_position_;
//                std::shared_ptr<Attribute<Vec3>> volume_color_;

                std::unique_ptr<cgogn::rendering::VolumeDrawer> vol_drw_;
                std::unique_ptr<cgogn::rendering::VolumeDrawer::Renderer> vol_rend_;

                bool render_volumes_;
                bool render_volumes_edges_;
                GLColor color_volumes_;
                GLColor color_edges_;
		float32 explode_factor_;

		bool auto_update_center_;
	};

public:

	VolumeRender(const App& app) :
		ViewModule(app, "VolumeRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		selected_view_(app.current_view()),
		selected_mesh_(nullptr)
	{}

	~VolumeRender()
	{}

private:

	void init_mesh(MESH* m)
	{
		app_.foreach_view([this, m] (View* v)
		{
			parameters_[v][m];
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m] (Attribute<Vec3>* attribute)
					{
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
                                                {
                                                    p.vol_drw_->update_edge(*m, p.vertex_position_.get());
                                                    p.vol_drw_->update_face(*m, p.vertex_position_.get());
                                                }
						v->request_update();
					}
				)
			);
		});
	}

public:

	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>* md = mesh_provider_->mesh_data(&m);
		p.vertex_position_ = vertex_position;
		p.vol_drw_->update_edge(m, p.vertex_position_.get());
		p.vol_drw_->update_face(m, p.vertex_position_.get());

		v.request_update();
	}


protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this] (MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
				mesh_provider_, this, &VolumeRender<MESH>::init_mesh
			)
		);
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>* md = mesh_provider_->mesh_data(m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

                        p.vol_rend_->set_explode_volume(p.explode_factor_);
                        if (p.render_volumes_)
			{
                            p.vol_rend_->set_face_color(p.color_volumes_);
                            glEnable(GL_POLYGON_OFFSET_FILL);
                            glPolygonOffset(1.0f, 2.0f);
                            p.vol_rend_->draw_faces(proj_matrix, view_matrix);
                            glDisable(GL_POLYGON_OFFSET_FILL);
			}

                        if (p.render_volumes_edges_)
			{
                            p.vol_rend_->set_edge_color(p.color_edges_);
                            p.vol_rend_->draw_edges(proj_matrix, view_matrix);
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
			app_.foreach_view([this] (View* v)
			{
				bool is_selected = v == selected_view_;
				if (ImGui::Selectable(v->name().c_str(), is_selected))
					selected_view_ = v;
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			});
			ImGui::EndCombo();
		}

		if (ImGui::ListBoxHeader("Mesh"))
		{
			mesh_provider_->foreach_mesh([this] (MESH* m, const std::string& name)
			{
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
					selected_mesh_ = m;
			});
			ImGui::ListBoxFooter();
		}

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&] (const std::shared_ptr<Attribute<Vec3>>& attribute)
				{
					bool is_selected = attribute == p.vertex_position_;
					if (ImGui::Selectable(attribute->name().c_str(), is_selected))
						set_vertex_position(*selected_view_, *selected_mesh_, attribute);
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				});
				ImGui::EndCombo();
			}
//			if (p.vertex_position_)
//			{
//			set_vertex_position(*selected_view_, *selected_mesh_, nullptr);
//			}


			ImGui::Separator();
                        need_update |= ImGui::Checkbox("Edges", &p.render_volumes_edges_);
                        need_update |= ImGui::Checkbox("Volumes", &p.render_volumes_);

                        if (p.render_volumes_)
			{
                            //color
                            need_update |= ImGui::ColorEdit3("color##volumes", p.color_volumes_.data(), ImGuiColorEditFlags_NoInputs);
			}

                        if (p.render_volumes_edges_)
			{
                            // color
                            need_update |= ImGui::ColorEdit3("color##edges", p.color_edges_.data(), ImGuiColorEditFlags_NoInputs);
                        }
                        ImGui::Separator();
                        need_update |= ImGui::SliderFloat("explode",&(p.explode_factor_),0.01,1.0);
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

#endif // CGOGN_MODULE_VOLUME_RENDER_H_
