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

#include <cgogn/modeling/algos/decimation/decimation.h>
#include <cgogn/modeling/algos/subdivision.h>

#include <cgogn/rendering_pureGL/imgui_viewer.h>
#include <cgogn/rendering_pureGL/mesh_render.h>
#include <cgogn/rendering_pureGL/drawer.h>
#include <cgogn/rendering_pureGL/vbo_update.h>

#include <cgogn/rendering_pureGL/shaders/shader_flat.h>
#include <cgogn/rendering_pureGL/shaders/shader_phong.h>
#include <cgogn/rendering_pureGL/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering_pureGL/shaders/shader_bold_line.h>
#include <cgogn/rendering_pureGL/shaders/shader_point_sprite.h>

#include <cgogn/io/surface_import.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Mesh = cgogn::CMap2;

using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

using Vec3 = cgogn::geometry::Vec3;

template <typename T>
using AttributePtr = typename cgogn::mesh_traits<Mesh>::AttributePtr<T>;

class App;

class Viewer : public cgogn::rendering_pgl::ImGUIViewer
{
	friend class App;

public:

	Viewer();
	virtual ~Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	void import(const std::string& filename);
	void update_bb();
	void draw() override;
	void init() override;
	void key_press_event(int k) override;
	void close_event() override;

private:

	Mesh mesh_;
	cgogn::CellFilter<Mesh> filtered_mesh_;

	AttributePtr<Vec3> vertex_position_;
	AttributePtr<Vec3> vertex_normal_;

	Vec3 bb_min_, bb_max_;

	std::unique_ptr<cgogn::rendering_pgl::MeshRender> render_;

	std::unique_ptr<cgogn::rendering_pgl::VBO> vbo_position_;
	std::unique_ptr<cgogn::rendering_pgl::VBO> vbo_normal_;

	std::unique_ptr<cgogn::rendering_pgl::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderFlat::Param> param_flat_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderVectorPerVertex::Param> param_normal_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderPhong::Param> param_phong_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderPointSprite::Param> param_point_sprite_;

	std::unique_ptr<cgogn::rendering_pgl::DisplayListDrawer> drawer_;
	std::unique_ptr<cgogn::rendering_pgl::DisplayListDrawer::Renderer> drawer_rend_;

	double mel_;

	bool phong_rendering_;
	bool vertices_rendering_;
	bool edge_rendering_;
	bool normal_rendering_;
	bool bb_rendering_;
};

class App: public cgogn::rendering_pgl::ImGUIApp
{
	int current_view_;

public:

	App() : current_view_(0) {}
	Viewer* view() { return static_cast<Viewer*>(viewers_[current_view_]); }
	bool interface() override;
	void key_press_event(int k) override;
};

/*****************************************************************************/
/*                          App IMPLEMENTATION                               */
/*****************************************************************************/

bool App::interface()
{
	ImGui::SetCurrentContext(context_);
	ImGui::GetIO().FontGlobalScale = interface_scaling_;

	bool inr = false;

	ImGui::Begin("Control Window", nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0, 0});

	inr |= ImGui::Checkbox("BB", &view()->bb_rendering_);
	inr |= ImGui::Checkbox("Phong", &view()->phong_rendering_);
	inr |= ImGui::Checkbox("Vertices", &view()->vertices_rendering_);
	inr |= ImGui::Checkbox("Normals", &view()->normal_rendering_);
	inr |= ImGui::Checkbox("Edges", &view()->edge_rendering_);

	if (view()->phong_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Phong parameters");
		inr |= ImGui::ColorEdit3("front color##phong", view()->param_phong_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
		ImGui::SameLine();
		inr |= ImGui::ColorEdit3("back color##phong", view()->param_phong_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
		inr |= ImGui::SliderFloat("spec##phong", &(view()->param_phong_->specular_coef_), 10.0f, 1000.0f);
		inr |= ImGui::Checkbox("double side##phong", &(view()->param_phong_->double_side_));
	}
	else
	{
		ImGui::Separator();
		ImGui::Text("Flat parameters");
		inr |= ImGui::ColorEdit3("front color##flat", view()->param_flat_->front_color_.data(), ImGuiColorEditFlags_NoInputs);
		ImGui::SameLine();
		inr |= ImGui::ColorEdit3("back color##flat", view()->param_flat_->back_color_.data(), ImGuiColorEditFlags_NoInputs);
		inr |= ImGui::Checkbox("single side##flat", &(view()->param_flat_->bf_culling_));
	}

	if (view()->normal_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Normal parameters");
		inr |= ImGui::ColorEdit3("color##norm", view()->param_normal_->color_.data(), ImGuiColorEditFlags_NoInputs);
		inr |= ImGui::SliderFloat("length##norm", &(view()->param_normal_->length_), 0.01f, 0.5f);
	}

	if (view()->edge_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Edge parameters");
		inr |= ImGui::ColorEdit3("color##edge", view()->param_edge_->color_.data());
		inr |= ImGui::SliderFloat("Width##edge", &(view()->param_edge_->width_), 1.0f, 10.0f);
	}

	if (view()->vertices_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Vertices parameters");
		inr |= ImGui::ColorEdit3("color##vert", view()->param_point_sprite_->color_.data());
		inr |= ImGui::SliderFloat("Size##vert", &(view()->param_point_sprite_->size_), view()->mel_ / 12, view()->mel_ / 3);
	}

//	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

	ImGui::End();

	if (inr)
		view()->request_update();
	
	return inr;
}

void App::key_press_event(int k)
{
	switch(k)
	{
		case int('I'):
			if (shift_pressed())
				interface_scaling_ += 0.1f;
			else
				interface_scaling_ -= 0.1f;
			break;
		case int(' '):
			show_imgui_ = !show_imgui_;
			break;
		default:
			break;
	}
	request_interface_update();
	ImGUIApp::key_press_event(k);
}

/*****************************************************************************/
/*                       Viewer IMPLEMENTATION                               */
/*****************************************************************************/

Viewer::Viewer() :
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

Viewer::~Viewer()
{}

void Viewer::close_event()
{
	render_.reset();
	vbo_position_.reset();
	vbo_normal_.reset();
	cgogn::rendering_pgl::ShaderProgram::clean_all();
}

void Viewer::import(const std::string& filename)
{
	cgogn::io::import_OFF(mesh_, filename);

	vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(mesh_, "position");
	if (!vertex_position_)
	{
		std::cerr << "Viewer::import: Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_normal_ = cgogn::add_attribute<Vec3, Vertex>(mesh_, "normal");
	cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);

	mel_ = cgogn::geometry::mean_edge_length(mesh_, vertex_position_);
}

void Viewer::update_bb()
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

void Viewer::draw()
{
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	cgogn::rendering_pgl::GLMat4 proj = get_projection_matrix();
	cgogn::rendering_pgl::GLMat4 view = get_modelview_matrix();

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 2.0f);

	if (phong_rendering_)
	{
		param_phong_->bind(proj, view);
		render_->draw(cgogn::rendering_pgl::TRIANGLES);
		param_phong_->release();
	}
	else
	{
		param_flat_->bind(proj, view);
		render_->draw(cgogn::rendering_pgl::TRIANGLES);
		param_flat_->release();
	}

	glDisable(GL_POLYGON_OFFSET_FILL);

	if (vertices_rendering_)
	{
		param_point_sprite_->bind(proj, view);
		render_->draw(cgogn::rendering_pgl::POINTS);
		param_point_sprite_->release();
	}

	if (edge_rendering_)
	{
		param_edge_->bind(proj, view);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		render_->draw(cgogn::rendering_pgl::LINES);
		glDisable(GL_BLEND);
		param_edge_->release();
	}

	if (normal_rendering_)
	{
		param_normal_->bind(proj, view);
		render_->draw(cgogn::rendering_pgl::POINTS);
		param_normal_->release();
	}

	if (bb_rendering_)
	{
		drawer_rend_->draw(proj, view);
	}
}

void Viewer::init()
{ 
	glClearColor(0.1f, 0.1f, 0.3f, 0.0f);

	drawer_ = cgogn::make_unique<cgogn::rendering_pgl::DisplayListDrawer>();
	drawer_rend_= drawer_->generate_renderer();

	update_bb();

	Vec3 diagonal = bb_max_ - bb_min_;
	set_scene_radius(diagonal.norm() / 2.0f);
	Vec3 center = (bb_max_ + bb_min_) / 2.0f;
	set_scene_center(center);

	vbo_position_ = cgogn::make_unique<cgogn::rendering_pgl::VBO>(3);
	vbo_normal_ = cgogn::make_unique<cgogn::rendering_pgl::VBO>(3);

	cgogn::rendering_pgl::update_vbo(vertex_position_.get(), vbo_position_.get());
	cgogn::rendering_pgl::update_vbo(vertex_normal_.get(), vbo_normal_.get());

	render_ = cgogn::make_unique<cgogn::rendering_pgl::MeshRender>();

	render_->init_primitives(mesh_, cgogn::rendering_pgl::POINTS);
	render_->init_primitives(mesh_, cgogn::rendering_pgl::LINES);
	render_->init_primitives(mesh_, cgogn::rendering_pgl::TRIANGLES);

	param_point_sprite_ = cgogn::rendering_pgl::ShaderPointSprite::generate_param();
	param_point_sprite_->set_vbos(vbo_position_.get());
	param_point_sprite_->color_ = cgogn::rendering_pgl::GLColor(1, 0, 0, 1);
	param_point_sprite_->size_ = mel_ / 6.0f;

	param_edge_ = cgogn::rendering_pgl::ShaderBoldLine::generate_param();
	param_edge_->set_vbos(vbo_position_.get());
	param_edge_->color_ = cgogn::rendering_pgl::GLColor(1, 1, 0, 1);
	param_edge_->width_= 2.5f;

	param_flat_ =  cgogn::rendering_pgl::ShaderFlat::generate_param();
	param_flat_->set_vbos(vbo_position_.get());
	param_flat_->front_color_ = cgogn::rendering_pgl::GLColor(0, 0.8f, 0, 1);
	param_flat_->back_color_ = cgogn::rendering_pgl::GLColor(0, 0, 0.8f, 1);
	param_flat_->ambiant_color_ = cgogn::rendering_pgl::GLColor(0.1f, 0.1f, 0.1f, 1);

	param_normal_ = cgogn::rendering_pgl::ShaderVectorPerVertex::generate_param();
	param_normal_->set_vbos(vbo_position_.get(), vbo_normal_.get());
	param_normal_->color_ = cgogn::rendering_pgl::GLColor(0.8f, 0, 0.8f, 1);
	param_normal_->length_ = (bb_max_ - bb_min_).norm() / 50.0f;

	param_phong_ = cgogn::rendering_pgl::ShaderPhong::generate_param();
	param_phong_->set_vbos(vbo_position_.get(), vbo_normal_.get());
}

void Viewer::key_press_event(cgogn::int32 k)
{
	switch (k)
	{
		case int('D') : {
			cgogn::modeling::decimate(mesh_, vertex_position_, cgogn::uint32(0.1 * cgogn::nb_cells<Vertex>(mesh_)));
			// cgogn::modeling::decimate(filtered_mesh_, vertex_position_, cgogn::uint32(0.1 * cgogn::nb_cells<Vertex>(filtered_mesh_)));
			std::cout << "nbv: " << cgogn::nb_cells<Vertex>(mesh_) << std::endl;
			cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);
			mel_ = cgogn::geometry::mean_edge_length(mesh_, vertex_position_);
			param_point_sprite_->size_ = mel_ / 6.0f;
			cgogn::rendering_pgl::update_vbo(vertex_position_.get(), vbo_position_.get());
			cgogn::rendering_pgl::update_vbo(vertex_normal_.get(), vbo_normal_.get());
			render_->init_primitives(mesh_, cgogn::rendering_pgl::POINTS);
			render_->init_primitives(mesh_, cgogn::rendering_pgl::LINES);
			render_->init_primitives(mesh_, cgogn::rendering_pgl::TRIANGLES);
			update_bb();
			Vec3 diagonal = bb_max_ - bb_min_;
			set_scene_radius(diagonal.norm() / 2.0f);
			break;
		}
		case int('S') : {
			cgogn::modeling::subdivide(mesh_, vertex_position_);
			// cgogn::modeling::subdivide(filtered_mesh_, vertex_position_);
			std::cout << "nbv: " << cgogn::nb_cells<Vertex>(mesh_) << std::endl;
			cgogn::geometry::compute_normal(mesh_, vertex_position_, vertex_normal_);
			mel_ = cgogn::geometry::mean_edge_length(mesh_, vertex_position_);
			param_point_sprite_->size_ = mel_ / 6.0f;
			cgogn::rendering_pgl::update_vbo(vertex_position_.get(), vbo_position_.get());
			cgogn::rendering_pgl::update_vbo(vertex_normal_.get(), vbo_normal_.get());
			render_->init_primitives(mesh_, cgogn::rendering_pgl::POINTS);
			render_->init_primitives(mesh_, cgogn::rendering_pgl::LINES);
			render_->init_primitives(mesh_, cgogn::rendering_pgl::TRIANGLES);
			update_bb();
			Vec3 diagonal = bb_max_ - bb_min_;
			set_scene_radius(diagonal.norm() / 2.0f);
			break;
		}
		default:
			break;
	}

	request_update();
}

int main(int argc, char** argv)
{
	std::string filename;
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	else
		filename = std::string(argv[1]);

	App app;
	app.set_window_title("Subdivision");

	gl3wInit();
	
	Viewer view;
	view.import(filename);
	
	app.add_view(&view);
	
	return app.launch();
}
