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

#include <iostream>

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/normal.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/module.h>

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/drawer.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/rendering/shaders/shader_simple_color.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_frame2d.h>

#include <cgogn/io/surface_import.h>

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


class SimpleViewer : public cgogn::ui::Module
{
public:

	SimpleViewer();
	virtual ~SimpleViewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(SimpleViewer);

	void init();

	virtual void draw(cgogn::ui::View* view) override;
	virtual void interface() override;

	void import(const std::string& filename);

private:

	void update_bb();

	Mesh mesh_;
	AttributePtr<Vec3> vertex_position_;
	AttributePtr<Vec3> vertex_normal_;

	Vec3 bb_min_, bb_max_;

	std::unique_ptr<cgogn::rendering::MeshRender> render_;

	std::unique_ptr<cgogn::rendering::VBO> vbo_position_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_normal_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_color_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_sphere_size_;

	std::unique_ptr<cgogn::rendering::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;
	std::unique_ptr<cgogn::rendering::ShaderVectorPerVertex::Param> param_normal_;
	std::unique_ptr<cgogn::rendering::ShaderPhong::Param> param_phong_;
	std::unique_ptr<cgogn::rendering::ShaderPointSpriteColorSize::Param> param_point_sprite_;

	std::unique_ptr<cgogn::rendering::DisplayListDrawer> drawer_;
	std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> drawer_rend_;

	bool phong_rendering_;
	bool vertices_rendering_;
	bool edge_rendering_;
	bool normal_rendering_;
	bool bb_rendering_;
};

SimpleViewer::SimpleViewer() :
	mesh_(),
	render_(nullptr),
	vbo_position_(nullptr),
	vbo_normal_(nullptr),
	vbo_color_(nullptr),
	vbo_sphere_size_(nullptr),
	drawer_(nullptr),
	drawer_rend_(nullptr),
	phong_rendering_(true),
	vertices_rendering_(false),
	edge_rendering_(false),
	normal_rendering_(false),
	bb_rendering_(true)
{}

SimpleViewer::~SimpleViewer()
{}

void SimpleViewer::init()
{
	// drawer for simple old-school gl rendering
	drawer_ = std::make_unique<cgogn::rendering::DisplayListDrawer>();
	drawer_rend_= drawer_->generate_renderer();

	update_bb();
	Vec3 diagonal = bb_max_ - bb_min_;

	// create and fill VBO for positions
	vbo_position_ = std::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());

	// create and fill VBO for normals
	vbo_normal_ = std::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());

	// fill a color vbo with abs of normals
	vbo_color_ = std::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_color_.get(), [] (const Vec3& n) -> cgogn::geometry::Vec3f
	{
		return { float(std::abs(n[0])), float(std::abs(n[1])), float(std::abs(n[2])) };
	});

	// fill a sphere size vbo
	vbo_sphere_size_ = std::make_unique<cgogn::rendering::VBO>(1);
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_sphere_size_.get(), [&] (const Vec3& n) -> float
	{
		return diagonal.norm() / 1000.0 * (1 + 2 * std::abs(n[2]));
	});

	// map rendering object (primitive creation & sending to GPU)
	render_ = std::make_unique<cgogn::rendering::MeshRender>();

	render_->init_primitives(mesh_, cgogn::rendering::POINTS);
	render_->init_primitives(mesh_, cgogn::rendering::LINES);
	render_->init_primitives(mesh_, cgogn::rendering::TRIANGLES);

	// generation of one parameter set (for this shader) : vbo + uniforms
	param_point_sprite_ = cgogn::rendering::ShaderPointSpriteColorSize::generate_param();
	param_point_sprite_->set_vbos(vbo_position_.get(), vbo_color_.get(), vbo_sphere_size_.get());
//	param_point_sprite_->size_ = 0.01f;

	param_edge_ = cgogn::rendering::ShaderBoldLine::generate_param();
	param_edge_->set_vbos(vbo_position_.get());
	param_edge_->color_ =  cgogn::rendering::GLColor(1,1,1,1);
	param_edge_->width_= 2.5f;

	param_flat_ = cgogn::rendering::ShaderFlat::generate_param();
	param_flat_->set_vbos(vbo_position_.get());
	param_flat_->front_color_ =  cgogn::rendering::GLColor(0,0.8f,0,1);
	param_flat_->back_color_ =  cgogn::rendering::GLColor(0,0,0.8f,1);
	param_flat_->ambiant_color_ =  cgogn::rendering::GLColor(0.1f,0.1f,0.1f,1);

	param_normal_ = cgogn::rendering::ShaderVectorPerVertex::generate_param();
	param_normal_->set_vbos(vbo_position_.get(), vbo_normal_.get());
	param_normal_->color_ =  cgogn::rendering::GLColor(0.8f,0.8f,0.8f,1);
	param_normal_->length_ = diagonal.norm()/50;

	param_phong_ = cgogn::rendering::ShaderPhong::generate_param();
	param_phong_->front_color_ =  cgogn::rendering::GLColor(0,0.8f,0,1);
	param_phong_->back_color_ =  cgogn::rendering::GLColor(0,0,0.8f,1);
	param_phong_->ambiant_color_ =  cgogn::rendering::GLColor(0.1f,0.1f,0.1f,1);
	param_phong_->set_vbos(vbo_position_.get(), vbo_normal_.get());//, vbo_color_.get());
}

void SimpleViewer::draw(cgogn::ui::View* view)
{
	Vec3 diagonal = bb_max_ - bb_min_;
	view->set_scene_radius(diagonal.norm() / 2.0f);
	Vec3 center = (bb_max_ + bb_min_) / 2.0f;
	view->set_scene_center(center);

	glEnable(GL_DEPTH_TEST);
	glClearColor(0.25f, 0.25f, 0.29f, 1);
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

	if (bb_rendering_)
		drawer_rend_->draw(proj_matrix, view_matrix);
}

void SimpleViewer::interface()
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
		need_update |= ImGui::SliderFloat("Width##edge", &(param_edge_->width_), 1.0f, 10.0f);
	}

	// ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

	ImGui::End();

	if (need_update)
		for (cgogn::ui::View* v : linked_views_)
			v->request_update();
}

void SimpleViewer::import(const std::string& filename)
{
	cgogn::io::import_OFF(mesh_, filename);

	vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(mesh_, "position");
	if (!vertex_position_)
	{
		std::cerr << "Viewer::import: Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_normal_ = cgogn::get_attribute<Vec3, Vertex>(mesh_, "normal");
	if (!vertex_normal_)
	{
		vertex_normal_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "normal");
		cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);
	}
}

void SimpleViewer::update_bb()
{
	for (cgogn::uint32 i = 0; i < 3; ++i)
	{
		bb_min_[i] = std::numeric_limits<cgogn::float32>::max();
		bb_max_[i] = std::numeric_limits<cgogn::float32>::lowest();
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

	std::string filename2 = std::string(DEFAULT_MESH_PATH) + std::string("off/horse.off");

	cgogn::thread_start(0, 0);

	cgogn::ui::App app;
	app.set_window_title("Simple viewer");

	SimpleViewer sv;
	sv.import(filename);
	sv.init();

	app.link_module(&sv);

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&sv);
	
	return app.launch();
}
