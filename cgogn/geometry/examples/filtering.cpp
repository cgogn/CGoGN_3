/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_filter.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/curvature.h>
#include <cgogn/geometry/algos/filtering.h>

#include <cgogn/modeling/algos/subdivision.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/module.h>

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/drawer.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/io/surface_import.h>

#include <chrono>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Mesh = cgogn::CMap2;

using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

template <typename T>
using AttributePtr = typename cgogn::mesh_traits<Mesh>::AttributePtr<T>;

using namespace cgogn::numerics;


class Filtering : public cgogn::ui::Module
{
public:

	Filtering();
	virtual ~Filtering();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Filtering);

	void init();

	virtual void resize_event(int32 frame_width, int32 frame_height) override;
	virtual void close_event() override;

	virtual void mouse_press_event(cgogn::ui::View* view, int32 button, float64 x, float64 y) override;
	virtual void mouse_release_event(cgogn::ui::View* view, int32 button, float64 x, float64 y) override;
	virtual void mouse_dbl_click_event(cgogn::ui::View* view, int32 button, float64 x, float64 y) override;
	virtual void mouse_move_event(cgogn::ui::View* view, float64 x, float64 y) override;
	virtual void mouse_wheel_event(cgogn::ui::View* view, float64 x, float64 y) override;
	virtual void key_press_event(cgogn::ui::View* view, int32 key_code) override;
	virtual void key_release_event(cgogn::ui::View* view, int32 key_code) override;

	virtual void draw(cgogn::ui::View* view) override;

	virtual void interface() override;

	void import(const std::string& filename);
	void filter_mesh();
	void subdivide_mesh();

private:

	void update_bb();

	Mesh mesh_;
	cgogn::CellFilter<Mesh> filtered_mesh_;

	AttributePtr<Vec3> vertex_position_;
	AttributePtr<Vec3> vertex_position2_;
	AttributePtr<Vec3> vertex_normal_;

	AttributePtr<double> edge_angle_;

	AttributePtr<double> vertex_kmin_;
	AttributePtr<double> vertex_kmax_;
	AttributePtr<Vec3> vertex_Kmin_;
	AttributePtr<Vec3> vertex_Kmax_;
	AttributePtr<Vec3> vertex_Knormal_;

	Vec3 bb_min_, bb_max_;

	std::unique_ptr<cgogn::rendering::MeshRender> render_;

	std::unique_ptr<cgogn::rendering::VBO> vbo_position_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_normal_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_Kmin_;

	std::unique_ptr<cgogn::rendering::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;
	std::unique_ptr<cgogn::rendering::ShaderVectorPerVertex::Param> param_normal_;
	std::unique_ptr<cgogn::rendering::ShaderVectorPerVertex::Param> param_Kmin_;
	std::unique_ptr<cgogn::rendering::ShaderPhong::Param> param_phong_;
	std::unique_ptr<cgogn::rendering::ShaderPointSprite::Param> param_point_sprite_;

	std::unique_ptr<cgogn::rendering::DisplayListDrawer> drawer_;
	std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> drawer_rend_;

	double mel_;

	bool phong_rendering_;
	bool vertices_rendering_;
	bool edge_rendering_;
	bool normal_rendering_;
	bool bb_rendering_;
};





Filtering::Filtering() :
	mesh_(),
	filtered_mesh_(mesh_),
	phong_rendering_(true),
	vertices_rendering_(false),
	edge_rendering_(false),
	normal_rendering_(false),
	bb_rendering_(true)
{
	filtered_mesh_.set_filter<Vertex>([&] (Vertex v) -> bool
	{
		return cgogn::value<Vec3>(mesh_, vertex_position_, v)[0] < 0.0f;
	});
	filtered_mesh_.set_filter<Edge>([&] (Edge e) -> bool
	{
		std::vector<Vertex> vertices = cgogn::incident_vertices(mesh_, e);
		auto v = std::find_if(vertices.begin(), vertices.end(), [&] (Vertex v) { return cgogn::value<Vec3>(mesh_, vertex_position_, v)[0] < 0.0f; });
		return v != vertices.end();
	});
	filtered_mesh_.set_filter<Face>([&] (Face f) -> bool
	{
		std::vector<Vertex> vertices = cgogn::incident_vertices(mesh_, f);
		auto v = std::find_if(vertices.begin(), vertices.end(), [&] (Vertex v) { return cgogn::value<Vec3>(mesh_, vertex_position_, v)[0] < 0.0f; });
		return v != vertices.end();
	});
}

Filtering::~Filtering()
{}

void Filtering::init()
{
	glClearColor(0.1f, 0.1f, 0.3f, 0.0f);

	drawer_ = std::make_unique<cgogn::rendering::DisplayListDrawer>();
	drawer_rend_= drawer_->generate_renderer();

	update_bb();

	vbo_position_ = std::make_unique<cgogn::rendering::VBO>();
	vbo_normal_ = std::make_unique<cgogn::rendering::VBO>();
	vbo_Kmin_ = std::make_unique<cgogn::rendering::VBO>();

	cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());
	cgogn::rendering::update_vbo(vertex_Kmin_.get(), vbo_Kmin_.get());

	render_ = std::make_unique<cgogn::rendering::MeshRender>();

	render_->init_primitives(mesh_, cgogn::rendering::POINTS);
	render_->init_primitives(mesh_, cgogn::rendering::LINES);
	render_->init_primitives(mesh_, cgogn::rendering::TRIANGLES);

	param_point_sprite_ = cgogn::rendering::ShaderPointSprite::generate_param();
	param_point_sprite_->set_vbos(vbo_position_.get());
	param_point_sprite_->color_ = cgogn::rendering::GLColor(1, 0, 0, 1);
	param_point_sprite_->size_ = mel_ / 6.0f;

	param_edge_ = cgogn::rendering::ShaderBoldLine::generate_param();
	param_edge_->set_vbos(vbo_position_.get());
	param_edge_->color_ = cgogn::rendering::GLColor(1, 1, 0, 1);
	param_edge_->width_= 2.5f;

	param_flat_ =  cgogn::rendering::ShaderFlat::generate_param();
	param_flat_->set_vbos(vbo_position_.get());
	param_flat_->front_color_ = cgogn::rendering::GLColor(0, 0.8f, 0, 1);
	param_flat_->back_color_ = cgogn::rendering::GLColor(0, 0, 0.8f, 1);
	param_flat_->ambiant_color_ = cgogn::rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

	param_normal_ = cgogn::rendering::ShaderVectorPerVertex::generate_param();
	param_normal_->set_vbos(vbo_position_.get(), vbo_normal_.get());
	param_normal_->color_ = cgogn::rendering::GLColor(0.8f, 0, 0.8f, 1);
	param_normal_->length_ = (bb_max_ - bb_min_).norm() / 50.0f;

	param_Kmin_ = cgogn::rendering::ShaderVectorPerVertex::generate_param();
	param_Kmin_->set_vbos(vbo_position_.get(), vbo_Kmin_.get());
	param_Kmin_->color_ = cgogn::rendering::GLColor(0.8f, 0.8f, 0.8f, 1);
	param_Kmin_->length_ = (bb_max_ - bb_min_).norm() / 50.0f;

	param_phong_ = cgogn::rendering::ShaderPhong::generate_param();
	param_phong_->set_vbos(vbo_position_.get(), vbo_normal_.get());
}

void Filtering::resize_event(int32 frame_width, int32 frame_height)
{}

void Filtering::close_event()
{
	render_.reset();
	vbo_position_.reset();
	vbo_normal_.reset();
	vbo_Kmin_.reset();
	drawer_.reset();
	drawer_rend_.reset();
	cgogn::rendering::ShaderProgram::clean_all();
}

void Filtering::mouse_press_event(cgogn::ui::View* view, int32 button, float64 x, float64 y)
{}

void Filtering::mouse_release_event(cgogn::ui::View* view, int32 button, float64 x, float64 y)
{}

void Filtering::mouse_dbl_click_event(cgogn::ui::View* view, int32 button, float64 x, float64 y)
{}

void Filtering::mouse_move_event(cgogn::ui::View* view, float64 x, float64 y)
{}

void Filtering::mouse_wheel_event(cgogn::ui::View* view, float64 x, float64 y)
{}

void Filtering::key_press_event(cgogn::ui::View* view, cgogn::int32 k)
{
	switch (k)
	{
		case int('A') : {
			filter_mesh();
			break;
		}
		case int('S') : {
			subdivide_mesh();
			break;
		}
		default:
			break;
	}
}

void Filtering::key_release_event(cgogn::ui::View* view, int32 key_code)
{}

void Filtering::draw(cgogn::ui::View* view)
{
	Vec3 diagonal = bb_max_ - bb_min_;
	view->set_scene_radius(diagonal.norm() / 2.0f);
	Vec3 center = (bb_max_ + bb_min_) / 2.0f;
	view->set_scene_center(center);

	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	const cgogn::rendering::GLMat4& proj_matrix = view->projection_matrix();
	const cgogn::rendering::GLMat4& view_matrix = view->modelview_matrix();

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 2.0f);

	if (phong_rendering_)
	{
		param_phong_->bind(proj_matrix, view_matrix);
		render_->draw(cgogn::rendering::TRIANGLES);
		param_phong_->release();
	}
	else
	{
		param_flat_->bind(proj_matrix, view_matrix);
		render_->draw(cgogn::rendering::TRIANGLES);
		param_flat_->release();
	}
	
	glDisable(GL_POLYGON_OFFSET_FILL);

	if (vertices_rendering_)
	{
		param_point_sprite_->bind(proj_matrix, view_matrix);
		render_->draw(cgogn::rendering::POINTS);
		param_point_sprite_->release();
	}

	if (edge_rendering_)
	{
		param_edge_->bind(proj_matrix, view_matrix);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		render_->draw(cgogn::rendering::LINES);
		glDisable(GL_BLEND);
		param_edge_->release();
	}

	if (normal_rendering_)
	{
		param_normal_->bind(proj_matrix, view_matrix);
		render_->draw(cgogn::rendering::POINTS);
		param_normal_->release();
	}

	param_Kmin_->bind(proj_matrix, view_matrix);
	render_->draw(cgogn::rendering::POINTS);
	param_Kmin_->release();

	if (bb_rendering_)
		drawer_rend_->draw(proj_matrix, view_matrix);
}

void Filtering::interface()
{
	bool need_update = false;

	ImGui::Begin("Control Window", nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0, 0});

	need_update |= ImGui::Checkbox("BB", &bb_rendering_);
	need_update |= ImGui::Checkbox("Phong/Flat", &phong_rendering_);
	need_update |= ImGui::Checkbox("Vertices", &vertices_rendering_);
	need_update |= ImGui::Checkbox("Normals", &normal_rendering_);
	need_update |= ImGui::Checkbox("Edges", &edge_rendering_);

	if (phong_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Phong parameters");
		need_update |= ImGui::ColorEdit3("front color##phong", param_phong_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
		ImGui::SameLine();
		need_update |= ImGui::ColorEdit3("back color##phong", param_phong_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
		need_update |= ImGui::SliderFloat("spec##phong", &(param_phong_->specular_coef_), 10.0f, 1000.0f);
		need_update |= ImGui::Checkbox("double side##phong", &(param_phong_->double_side_));
	}
	else
	{
		ImGui::Separator();
		ImGui::Text("Flat parameters");
		need_update |= ImGui::ColorEdit3("front color##flat", param_flat_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
		ImGui::SameLine();
		need_update |= ImGui::ColorEdit3("back color##flat", param_flat_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
		need_update |= ImGui::Checkbox("single side##flat", &(param_flat_->bf_culling_));
	}

	if (normal_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Normal parameters");
		need_update |= ImGui::ColorEdit3("color##norm", param_normal_->color_.data(), ImGuiColorEditFlags_NoInputs);
		need_update |= ImGui::SliderFloat("length##norm", &(param_normal_->length_), 0.01f, 0.5f);
	}

	if (edge_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Edge parameters");
		need_update |= ImGui::ColorEdit3("color##edge", param_edge_->color_.data());
		need_update |= ImGui::SliderFloat("width##edge", &(param_edge_->width_), 1.0f, 10.0f);
	}

	if (vertices_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Vertices parameters");
		need_update |= ImGui::ColorEdit3("color##vert", param_point_sprite_->color_.data());
		need_update |= ImGui::SliderFloat("size##vert", &(param_point_sprite_->size_), mel_ / 12, mel_ / 3);
	}

	ImGui::Separator();
	if (ImGui::Button("Subdivide"))
		subdivide_mesh();
	if (ImGui::Button("Filter"))
		filter_mesh();

//	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

	ImGui::End();

	if (need_update)
		for (cgogn::ui::View* v : linked_views_)
			v->request_update();
}

void Filtering::import(const std::string& filename)
{
	cgogn::io::import_OFF(mesh_, filename);

	vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(mesh_, "position");
	if (!vertex_position_)
	{
		std::cerr << "Viewer::import: Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_position2_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "position2");

	vertex_normal_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "normal");

	std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
	
	t1 = std::chrono::high_resolution_clock::now();
	cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);
    t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "compute_normal: " << duration << " us" << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	mel_ = cgogn::geometry::mean_edge_length(mesh_, vertex_position_);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "mean_edge_length: " << duration << " us" << std::endl;

	edge_angle_ = cgogn::add_attribute<double, Edge>(mesh_, "angle");

	t1 = std::chrono::high_resolution_clock::now();
	cgogn::geometry::compute_angle(mesh_, vertex_position_, edge_angle_);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "compute_angle: " << duration << " us" << std::endl;

	vertex_kmin_ = cgogn::add_attribute<double, Vertex>(mesh_, "kmin");
	vertex_kmax_ = cgogn::add_attribute<double, Vertex>(mesh_, "kmax");
	vertex_Kmin_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "Kmin");
	vertex_Kmax_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "Kmax");
	vertex_Knormal_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "Knormal");

	t1 = std::chrono::high_resolution_clock::now();
	cgogn::geometry::compute_curvature(
		mesh_,
		mel_ * 4.0,
		vertex_position_,
		vertex_normal_,
		edge_angle_,
		vertex_kmax_,
		vertex_kmin_,
		vertex_Kmax_,
		vertex_Kmin_,
		vertex_Knormal_
	);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "compute_curvature: " << duration << " us" << std::endl;

	std::cout << std::endl;
}

void Filtering::filter_mesh()
{
	for (auto it = vertex_position2_->begin(), end = vertex_position2_->end(); it != end; ++it)
		*it = (*vertex_position_)[it.index()];
	cgogn::geometry::filter_average<Vec3>(filtered_mesh_, vertex_position_, vertex_position2_);
	vertex_position_->swap(vertex_position2_.get());

	std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
	
	t1 = std::chrono::high_resolution_clock::now();
	cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);
    t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "compute_normal: " << duration << " us" << std::endl;
	
	t1 = std::chrono::high_resolution_clock::now();
	mel_ = cgogn::geometry::mean_edge_length(mesh_, vertex_position_);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "mean_edge_length: " << duration << " us" << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	cgogn::geometry::compute_angle(mesh_, vertex_position_, edge_angle_);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "compute_angle: " << duration << " us" << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	cgogn::geometry::compute_curvature(
		mesh_,
		mel_ * 4.0,
		vertex_position_,
		vertex_normal_,
		edge_angle_,
		vertex_kmax_,
		vertex_kmin_,
		vertex_Kmax_,
		vertex_Kmin_,
		vertex_Knormal_
	);
    t2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "compute_curvature: " << duration << " us" << std::endl;

	std::cout << std::endl;

	cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());
	cgogn::rendering::update_vbo(vertex_Kmin_.get(), vbo_Kmin_.get());

	update_bb();
	Vec3 diagonal = bb_max_ - bb_min_;

	param_point_sprite_->size_ = mel_ / 6.0f;
	param_normal_->length_ = diagonal.norm() / 50.0f;
	param_Kmin_->length_ = diagonal.norm() / 50.0f;

	for (cgogn::ui::View* v : linked_views_)
		v->request_update();
}

void Filtering::subdivide_mesh()
{
	// cgogn::CellCache<Mesh> cached_map(mesh_);
	// cached_map.build<Vertex>();
	// cached_map.build<Edge>();
	// cached_map.build<Face>();

	// cgogn::CellFilter<cgogn::CellCache<Mesh>> filtered_map(cached_map);
	// filtered_map.set_filter<Edge>([&] (Edge e) -> bool
	// {
	// 	std::vector<Vertex> vertices = cgogn::incident_vertices(cached_map, e);
	// 	auto v = std::find_if(vertices.begin(), vertices.end(), [&] (Vertex v) { return cgogn::value<Vec3>(cached_map, vertex_position_, v)[0] < 0.0f; });
	// 	return v != vertices.end();
	// });
	// filtered_map.set_filter<Face>([&] (Face f) -> bool
	// {
	// 	std::vector<Vertex> vertices = cgogn::incident_vertices(cached_map, f);
	// 	auto v = std::find_if(vertices.begin(), vertices.end(), [&] (Vertex v) { return cgogn::value<Vec3>(cached_map, vertex_position_, v)[0] < 0.0f; });
	// 	return v != vertices.end();
	// });

	cgogn::modeling::subdivide(filtered_mesh_, vertex_position_);
	std::cout << "nbv: " << cgogn::nb_cells<Vertex>(mesh_) << std::endl;

	render_->init_primitives(mesh_, cgogn::rendering::POINTS);
	render_->init_primitives(mesh_, cgogn::rendering::LINES);
	render_->init_primitives(mesh_, cgogn::rendering::TRIANGLES);

	cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);
	mel_ = cgogn::geometry::mean_edge_length(mesh_, vertex_position_);
	cgogn::geometry::compute_angle(mesh_, vertex_position_, edge_angle_);
	cgogn::geometry::compute_curvature(
		mesh_,
		mel_ * 4.0,
		vertex_position_,
		vertex_normal_,
		edge_angle_,
		vertex_kmax_,
		vertex_kmin_,
		vertex_Kmax_,
		vertex_Kmin_,
		vertex_Knormal_
	);

	cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());
	cgogn::rendering::update_vbo(vertex_Kmin_.get(), vbo_Kmin_.get());
	
	update_bb();
	Vec3 diagonal = bb_max_ - bb_min_;

	param_point_sprite_->size_ = mel_ / 6.0f;
	param_normal_->length_ = diagonal.norm() / 50.0f;
	param_Kmin_->length_ = diagonal.norm() / 50.0f;

	for (cgogn::ui::View* v : linked_views_)
		v->request_update();
}

void Filtering::update_bb()
{
	for (cgogn::uint32 i = 0; i < 3; ++i)
	{
		bb_min_[i] = std::numeric_limits<cgogn::float64>::max();
		bb_max_[i] = std::numeric_limits<cgogn::float64>::lowest();
	}
	for (const Vec3& p : *vertex_position_)
	{
		for (cgogn::uint32 i = 0; i < 3; ++i)
		{
			if (p[i] < bb_min_[i])
				bb_min_[i] = p[i];
			if (p[i] > bb_max_[i])
				bb_max_[i] = p[i];
		}
	}

	drawer_->new_list();
	drawer_->line_width_aa(2.0);
	drawer_->begin(GL_LINE_LOOP);
		drawer_->color3f(1.0, 1.0, 1.0);
		drawer_->vertex3f(bb_min_[0], bb_min_[1], bb_min_[2]);
		drawer_->vertex3f(bb_max_[0], bb_min_[1], bb_min_[2]);
		drawer_->vertex3f(bb_max_[0], bb_max_[1], bb_min_[2]);
		drawer_->vertex3f(bb_min_[0], bb_max_[1], bb_min_[2]);
		drawer_->vertex3f(bb_min_[0], bb_max_[1], bb_max_[2]);
		drawer_->vertex3f(bb_max_[0], bb_max_[1], bb_max_[2]);
		drawer_->vertex3f(bb_max_[0], bb_min_[1], bb_max_[2]);
		drawer_->vertex3f(bb_min_[0], bb_min_[1], bb_max_[2]);
	drawer_->end();
	drawer_->begin(GL_LINES);
	drawer_->color3f(1.0, 1.0, 1.0);
		drawer_->vertex3f(bb_min_[0], bb_min_[1], bb_min_[2]);
		drawer_->vertex3f(bb_min_[0], bb_max_[1], bb_min_[2]);
		drawer_->vertex3f(bb_min_[0], bb_min_[1], bb_max_[2]);
		drawer_->vertex3f(bb_min_[0], bb_max_[1], bb_max_[2]);
		drawer_->vertex3f(bb_max_[0], bb_min_[1], bb_min_[2]);
		drawer_->vertex3f(bb_max_[0], bb_min_[1], bb_max_[2]);
		drawer_->vertex3f(bb_max_[0], bb_max_[1], bb_min_[2]);
		drawer_->vertex3f(bb_max_[0], bb_max_[1], bb_max_[2]);
	drawer_->end();
	drawer_->end_list();
}




int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	cgogn::thread_start(0, 0);

	cgogn::ui::App app;
	app.set_window_title("Filtering");

	Filtering f;
	f.import(filename);
	f.init();

	app.link_module(&f);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&f);

	cgogn::ui::View* v2 = app.add_view();
	v2->link_module(&f);

	// cgogn::ui::View* v3 = app.add_view();
	// v3->link_module(&f);

	// cgogn::ui::View* v4 = app.add_view();
	// v4->link_module(&f);

	return app.launch();
}
