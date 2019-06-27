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

#include <cgogn/rendering_pureGL/imgui_viewer.h>
#include <GLFW/glfw3.h>
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

#include <cgogn/rendering_pureGL//mesh_render.h>
#include <cgogn/rendering_pureGL/vbo_update.h>

#include <cgogn/rendering_pureGL/shaders/shader_flat.h>
#include <cgogn/rendering_pureGL/shaders/shader_phong.h>
#include <cgogn/rendering_pureGL/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering_pureGL/shaders/shader_bold_line.h>
#include <cgogn/rendering_pureGL/shaders/shader_point_sprite.h>

#include <cgogn/io/surface_import.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Map2 = cgogn::CMap2;
using Vertex = Map2::Vertex;
using Edge = Map2::Edge;
using Face = Map2::Face;

using Vec3 = cgogn::geometry::Vec3;

template <typename T>
using AttributePtr = typename cgogn::mesh_traits<Map2>::AttributePtr<T>;

namespace GL = ::cgogn::rendering_pgl;

class App;

class Viewer : public GL::ImGUIViewer
{
	friend class App;

public:

	Viewer();
	virtual ~Viewer();

	Viewer(const Viewer&) = delete;
	Viewer& operator=(const Viewer&) = delete;

	void import(const std::string& surface_mesh);
	void update_bb();

	void draw() override;
	void init() override;
	void key_press_event(int k) override;
	void close_event() override;

private:

	Map2 map_;
	cgogn::CellFilter<Map2> filtered_map_;

	AttributePtr<Vec3> vertex_position_;
	AttributePtr<Vec3> vertex_normal_;

	Vec3 bb_min_, bb_max_;

	std::unique_ptr<GL::MeshRender> render_;

	std::unique_ptr<GL::VBO> vbo_position_;
	std::unique_ptr<GL::VBO> vbo_normal_;

	std::unique_ptr<cgogn::rendering_pgl::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderFlat::Param> param_flat_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderVectorPerVertex::Param> param_normal_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderPhong::Param> param_phong_;
	std::unique_ptr<cgogn::rendering_pgl::ShaderPointSprite::Param> param_point_sprite_;

	double mel_;

	bool phong_rendering_;
	bool flat_rendering_;
	bool vertices_rendering_;
	bool edge_rendering_;
	bool normal_rendering_;
};

class App: public GL::ImGUIApp
{
	int current_view_;

public:

	App():	current_view_(0) {}
	Viewer* view() { return static_cast<Viewer*>(viewers_[current_view_]); }
	bool interface() override;
	void key_press_event(int k) override;
};


////////////////////
// IMPLEMENTATION //
////////////////////

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
	ask_interface_update();
	ImGUIApp::key_press_event(k);
}

bool App::interface()
{
	std::cout << "draw interface" << std::endl;

	ImGui::SetCurrentContext(context_);
	ImGui::GetIO().FontGlobalScale = interface_scaling_;

	bool inr = false;

	ImGui::Begin("Control Window",nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0,0});

	inr |= ImGui::Checkbox("Phong/Flat", &view()->phong_rendering_);
	inr |= ImGui::Checkbox("Vertices", &view()->vertices_rendering_);
	inr |= ImGui::Checkbox("Normals", &view()->normal_rendering_);
	inr |= ImGui::Checkbox("Edges", &view()->edge_rendering_);
	inr |= ImGui::Checkbox("Vertice", &view()->vertices_rendering_);

	if (view()->phong_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Phong parameters");
		inr |= ImGui::ColorEdit3("front color##phong",view()->param_phong_->front_color_.data(),ImGuiColorEditFlags_NoInputs);
		ImGui::SameLine();
		inr |= ImGui::ColorEdit3("back color##phong",view()->param_phong_->back_color_.data(),ImGuiColorEditFlags_NoInputs);
		inr |= ImGui::SliderFloat("spec##phong", &(view()->param_phong_->specular_coef_), 10.0f, 1000.0f);
		inr |= ImGui::Checkbox("double side##phong", &(view()->param_phong_->double_side_));
	}
	else
	{
		ImGui::Separator();
		ImGui::Text("Flat parameters");
		inr |= ImGui::ColorEdit3("front color##flat",view()->param_flat_->front_color_.data(),ImGuiColorEditFlags_NoInputs);
		ImGui::SameLine();
		inr |= ImGui::ColorEdit3("back color##flat",view()->param_flat_->back_color_.data(),ImGuiColorEditFlags_NoInputs);
		inr |= ImGui::Checkbox("single side##flat", &(view()->param_flat_->bf_culling_));
	}
	if (view()->normal_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Normal parameters");
		inr |= ImGui::ColorEdit3("color##norm",view()->param_normal_->color_.data(),ImGuiColorEditFlags_NoInputs);
		inr |= ImGui::SliderFloat("length##norm", &(view()->param_normal_->length_), 0.01f, 0.5f);
	}

	if (view()->edge_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Edge parameters");
		inr |= ImGui::ColorEdit3("color##edge",view()->param_edge_->color_.data());
		inr |= ImGui::SliderFloat("Width##edge", &(view()->param_edge_->width_), 1.0f, 10.0f);
	}

	if (view()->vertices_rendering_)
	{
		ImGui::Separator();
		ImGui::Text("Vertices parameters");
		inr |= ImGui::ColorEdit3("color##vert",view()->param_point_sprite_->color_.data());
		inr |= ImGui::SliderFloat("Size##vert", &(view()->param_point_sprite_->size_), view()->mel_/12, view()->mel_/3);
	}

	//	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

	ImGui::End();

	if (inr)
	{
		view()->ask_update();
	}
	return inr;
}


Viewer::Viewer() :
	map_(),
	filtered_map_(map_),
	phong_rendering_(true),
	flat_rendering_(false),
	vertices_rendering_(false),
	edge_rendering_(false),
	normal_rendering_(false)
{
	filtered_map_.set_filter<Vertex>([&] (Vertex v) -> bool
	{
		return cgogn::value<Vec3>(map_, vertex_position_, v)[0] < 0.0f;
	});
	filtered_map_.set_filter<Edge>([&] (Edge e) -> bool
	{
		std::vector<Vertex> vertices = cgogn::incident_vertices(map_, e);
		auto v = std::find_if(vertices.begin(), vertices.end(), [&] (Vertex v) { return cgogn::value<Vec3>(map_, vertex_position_, v)[0] < 0.0f; });
		return v != vertices.end();
	});
	filtered_map_.set_filter<Face>([&] (Face f) -> bool
	{
		std::vector<Vertex> vertices = cgogn::incident_vertices(map_, f);
		auto v = std::find_if(vertices.begin(), vertices.end(), [&] (Vertex v) { return cgogn::value<Vec3>(map_, vertex_position_, v)[0] < 0.0f; });
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

void Viewer::import(const std::string& surface_mesh)
{
	cgogn::io::import_OFF(map_, surface_mesh);

	vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(map_, "position");
	if (!vertex_position_)
	{
		std::cerr << "Viewer::import: Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_normal_ = cgogn::add_attribute<Vec3, Vertex>(map_, "normal");
	cgogn::geometry::compute_normal(map_, vertex_position_, vertex_normal_);

	update_bb();

	Vec3 diagonal = bb_max_ - bb_min_;
	set_scene_radius(diagonal.norm() / 2.0f);
	Vec3 center = (bb_max_ + bb_min_) / 2.0f;
	set_scene_center(center);
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
}

void Viewer::key_press_event(cgogn::int32 k)
{
	switch(k)
	{
		case int('D'): {
			cgogn::modeling::decimate(filtered_map_, vertex_position_, cgogn::uint32(0.1 * cgogn::nb_cells<Vertex>(filtered_map_)));
			std::cout << "nbv: " << cgogn::nb_cells<Vertex>(map_) << std::endl;
			cgogn::geometry::compute_normal(map_, vertex_position_, vertex_normal_);
			double mel = cgogn::geometry::mean_edge_length(map_, vertex_position_);
			param_point_sprite_->size_ = mel / 6.0f;
			GL::update_vbo(vertex_position_.get(), vbo_position_.get());
			GL::update_vbo(vertex_normal_.get(), vbo_normal_.get());
			render_->init_primitives(map_, GL::POINTS);
			render_->init_primitives(map_, GL::LINES);
			render_->init_primitives(map_, GL::TRIANGLES);
			update_bb();
			Vec3 diagonal = bb_max_ - bb_min_;
			set_scene_radius(diagonal.norm() / 2.0f);
			break;
		}
		case int('S'): {
			cgogn::modeling::subdivide(filtered_map_, vertex_position_);
			std::cout << "nbv: " << cgogn::nb_cells<Vertex>(map_) << std::endl;
			cgogn::geometry::compute_normal(map_, vertex_position_, vertex_normal_);
			mel_ = cgogn::geometry::mean_edge_length(map_, vertex_position_);
			param_point_sprite_->size_ = mel_ / 6.0f;
			GL::update_vbo(vertex_position_.get(), vbo_position_.get());
			GL::update_vbo(vertex_normal_.get(), vbo_normal_.get());
			render_->init_primitives(map_, GL::POINTS);
			render_->init_primitives(map_, GL::LINES);
			render_->init_primitives(map_, GL::TRIANGLES);
			update_bb();
			Vec3 diagonal = bb_max_ - bb_min_;
			set_scene_radius(diagonal.norm() / 2.0f);
			break;
		}

		default:
			break;
	}
	ask_update();
}

void Viewer::draw()
{
	glEnable(GL_DEPTH_TEST);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	GL::GLMat4 proj = get_projection_matrix();
	GL::GLMat4 view = get_modelview_matrix();

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 2.0f);
	if (flat_rendering_)
	{
		param_flat_->bind(proj, view);
		render_->draw(GL::TRIANGLES);
		param_flat_->release();
	}

	if (phong_rendering_)
	{
		param_phong_->bind(proj, view);
		render_->draw(GL::TRIANGLES);
		param_phong_->release();
	}
	glDisable(GL_POLYGON_OFFSET_FILL);

	if (vertices_rendering_)
	{
		param_point_sprite_->bind(proj, view);
		render_->draw(GL::POINTS);
		param_point_sprite_->release();
	}

	if (edge_rendering_)
	{
		param_edge_->bind(proj, view);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		render_->draw(GL::LINES);
		glDisable(GL_BLEND);
		param_edge_->release();
	}

	if (normal_rendering_)
	{
		param_normal_->bind(proj, view);
		render_->draw(GL::POINTS);
		param_normal_->release();
	}
}

void Viewer::init()
{
	glClearColor(0.1f,0.1f,0.3f,0.0f);

	vbo_position_ = cgogn::make_unique<GL::VBO>();
	vbo_normal_ = cgogn::make_unique<GL::VBO>();

	GL::update_vbo(vertex_position_.get(), vbo_position_.get());
	GL::update_vbo(vertex_normal_.get(), vbo_normal_.get());

	render_ = cgogn::make_unique<GL::MeshRender>();

	render_->init_primitives(map_, GL::POINTS);
	render_->init_primitives(map_, GL::LINES);
	render_->init_primitives(map_, GL::TRIANGLES);

	param_point_sprite_ = GL::ShaderPointSprite::generate_param();
	param_point_sprite_->set_vbos(vbo_position_.get());
	param_point_sprite_->color_ = GL::GLColor(1, 0, 0,1);
	double mel = cgogn::geometry::mean_edge_length(map_, vertex_position_);
	param_point_sprite_->size_ = mel / 6.0f;

	param_edge_ = GL::ShaderBoldLine::generate_param();
	param_edge_->set_vbos(vbo_position_.get());
	param_edge_->color_ = GL::GLColor(1,1,0,1);
	param_edge_->width_= 2.5f;

	param_flat_ =  GL::ShaderFlat::generate_param();
	param_flat_->set_vbos(vbo_position_.get());
	param_flat_->front_color_ = GL::GLColor(0,0.8f,0,1);
	param_flat_->back_color_ = GL::GLColor(0,0,0.8f,1);
	param_flat_->ambiant_color_ = GL::GLColor(0.1f,0.1f,0.1f,1);

	param_normal_ = GL::ShaderVectorPerVertex::generate_param();
	param_normal_->set_vbos(vbo_position_.get(), vbo_normal_.get());
	param_normal_->color_ = GL::GLColor(0.8f,0,0.8f,1);
	param_normal_->length_ = (bb_max_ - bb_min_).norm() / 50.0f;

	param_phong_ = GL::ShaderPhong::generate_param();
	param_phong_->set_vbos(vbo_position_.get(), vbo_normal_.get());
}

int main(int argc, char** argv)
{
	std::string surface_mesh;
	if (argc < 2)
	{
		surface_mesh = std::string(DEFAULT_MESH_PATH) + std::string("off/socket.off");
	}
	else
		surface_mesh = std::string(argv[1]);

	App app;
	app.set_window_title("Subdivision");
	gl3wInit();
	Viewer view;
	view.import(surface_mesh);
	app.add_view(&view);
	return app.launch();
}
