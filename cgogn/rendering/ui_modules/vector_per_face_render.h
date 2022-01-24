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

#ifndef CGOGN_MODULE_VECTOR_PER_FACE_RENDER_H_
#define CGOGN_MODULE_VECTOR_PER_FACE_RENDER_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class VectorPerFaceRender : public ViewModule
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_position_vbo_(nullptr), face_vector_(nullptr),
			  face_vector_vbo_(nullptr), face_center_(nullptr), face_center_vbo_(nullptr), vector_scale_factor_(1.0),
			  vector_base_size_(1.0)
		{
			param_vector_per_vertex_ = rendering::ShaderVectorPerVertex::generate_param();
			param_vector_per_vertex_->color_ = {1.0f, 0.0f, 0.0f, 1.0f};
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vertex_position_vbo_;
		std::shared_ptr<Attribute<Vec3>> face_vector_;
		rendering::VBO* face_vector_vbo_;

		std::shared_ptr<Attribute<Vec3>> face_center_;
		rendering::VBO* face_center_vbo_;

		std::unique_ptr<rendering::ShaderVectorPerVertex::Param> param_vector_per_vertex_;

		float32 vector_scale_factor_;
		float32 vector_base_size_;
	};

public:
	VectorPerFaceRender(const App& app)
		: ViewModule(app, "VectorPerFaceRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~VectorPerFaceRender()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		for (View* v : linked_views_)
		{
			Parameters& p = parameters_[v][m];

			p.face_center_ = add_attribute<Vec3, Face>(*m, "__face_center");

			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);

			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, v, m]() {
					Parameters& p = parameters_[v][m];
					if (p.vertex_position_)
					{
						p.vector_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 2.0);
						update_face_center(*v, *m);
					}
					if (p.vector_base_size_ == 0.0)
						p.vector_base_size_ = 1.0;
					v->request_update();
				}));
			mesh_connections_[m].push_back(
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					m, [this, v, m](Attribute<Vec3>* attribute) {
						Parameters& p = parameters_[v][m];
						if (p.vertex_position_.get() == attribute)
						{
							p.vector_base_size_ =
								float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 2.0);
							update_face_center(*v, *m);
							if (p.vector_base_size_ == 0.0)
								p.vector_base_size_ = 1.0;
						}
						v->request_update();
					}));
		}
	}

public:
	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_position_ == vertex_position)
			return;

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.vector_base_size_ = float32(geometry::mean_edge_length(m, vertex_position.get()) / 2.0);
			if (p.vector_base_size_ == 0.0)
				p.vector_base_size_ = 1.0;
			update_face_center(v, m);
		}
		else
			p.vertex_position_vbo_ = nullptr;

		p.param_vector_per_vertex_->set_vbos({p.face_center_vbo_, p.face_vector_vbo_});

		v.request_update();
	}

	void set_face_vector(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& face_vector)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.face_vector_ == face_vector)
			return;

		p.face_vector_ = face_vector;
		if (p.face_vector_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			p.face_vector_vbo_ = md.update_vbo(face_vector.get(), true);
		}
		else
			p.face_vector_vbo_ = nullptr;

		p.param_vector_per_vertex_->set_vbos({p.face_center_vbo_, p.face_vector_vbo_});

		v.request_update();
	}

protected:
	void update_face_center(View& v, const MESH& m)
	{
		Parameters& p = parameters_[&v][&m];
		MeshData<MESH>& md = mesh_provider_->mesh_data(m);

		if (p.vertex_position_)
		{
			geometry::compute_centroid<Vec3, Face>(m, p.vertex_position_.get(), p.face_center_.get());
			p.face_center_vbo_ = md.update_vbo(p.face_center_.get(), true);
		}
	}

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &VectorPerFaceRender<MESH>::init_mesh));
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_[view])
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.param_vector_per_vertex_->attributes_initialized())
			{
				p.param_vector_per_vertex_->length_ = p.vector_base_size_ * p.vector_scale_factor_;
				p.param_vector_per_vertex_->bind(proj_matrix, view_matrix);
				md.draw(rendering::INDEX_FACES);
				p.param_vector_per_vertex_->release();
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Mesh", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_view_, *selected_mesh_, attribute);
												});

			imgui_combo_attribute<Face, Vec3>(*selected_mesh_, p.face_vector_, "Vector",
											  [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
												  set_face_vector(*selected_view_, *selected_mesh_, attribute);
											  });

			ImGui::Separator();
			ImGui::TextUnformatted("Vectors");
			need_update |=
				ImGui::ColorEdit3("color", p.param_vector_per_vertex_->color_.data(), ImGuiColorEditFlags_NoInputs);
			need_update |= ImGui::SliderFloat("length", &p.vector_scale_factor_, 0.1f, 5.0f);
			need_update |= ImGui::SliderFloat("width", &p.param_vector_per_vertex_->width_, 1.0f, 9.0f);
			need_update |= ImGui::SliderFloat("lighted", &p.param_vector_per_vertex_->lighted_, 0.0f, 1.0f);

			float64 remain = mesh_provider_->mesh_data(*selected_mesh_).outlined_until_ - App::frame_time_;
			if (remain > 0)
				need_update = true;

			if (need_update)
				for (View* v : linked_views_)
					v->request_update();
		}
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

#endif // CGOGN_MODULE_VECTOR_PER_FACE_RENDER_H_
