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
#include <iostream>

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/io/surface_import.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/rendering_pureGL/map_render.h>
#include <cgogn/rendering_pureGL/shaders/shader_simple_color.h>
#include <cgogn/rendering_pureGL/shaders/shader_flat.h>
#include <cgogn/rendering_pureGL/shaders/shader_phong.h>
#include <cgogn/rendering_pureGL/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering_pureGL/shaders/shader_bold_line.h>
#include <cgogn/rendering_pureGL/shaders/shader_point_sprite.h>
#include <cgogn/rendering_pureGL/shaders/shader_frame2d.h>
#include <cgogn/rendering_pureGL/drawer.h>
#include <cgogn/rendering_pureGL/vbo_update.h>
#include <cgogn/rendering_pureGL/fbo.h>
#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Map2 = cgogn::CMap2;
using Vec3 = Eigen::Vector3f;

template <typename T>
using AttributePtr = typename cgogn::mesh_traits<Map2>::AttributePtr<T>;


using namespace cgogn;
namespace GL = ::cgogn::rendering_pgl;

class App;

class Viewer : public GL::ImGUIViewer
{
	friend class App;
public:
	Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	void import(const std::string& surface_mesh);

	Map2 map_;
//	cgogn::geometry::AABB<Vec3> bb_;
	std::unique_ptr<GL::MeshRender> render_;
	AttributePtr<Vec3> vertex_position_;
	AttributePtr<Vec3> vertex_normal_;

	std::unique_ptr<GL::VBO> vbo_pos_;
	std::unique_ptr<GL::VBO> vbo_norm_;
	std::unique_ptr<GL::VBO> vbo_color_;
	std::unique_ptr<GL::VBO> vbo_sphere_sz_;

	std::unique_ptr<GL::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<GL::ShaderFlat::Param> param_flat_;
	std::unique_ptr<GL::ShaderVectorPerVertex::Param> param_normal_;
	std::unique_ptr<GL::ShaderPhong::Param> param_phong_;
	std::unique_ptr<GL::ShaderPointSpriteColorSize::Param> param_point_sprite_;

	std::unique_ptr<GL::ShaderFrame2d::Param> param_frame_;
	std::unique_ptr<GL::DisplayListDrawer> drawer_;
	std::unique_ptr<GL::DisplayListDrawer::Renderer> drawer_rend_;

	std::unique_ptr<GL::ShaderFSTexture::Param> param_fst_;
//	std::unique_ptr<GL::FBO> fbo_;
//	std::unique_ptr<GL::Texture2D> tex_;

	Vec3 bb_min_, bb_max_;


	bool phong_rendering_;
	bool vertices_rendering_;
	bool edge_rendering_;
	bool normal_rendering_;
	bool bb_rendering_;
public:
	void update_bb();
	void draw() override;
	void init() override;
	void key_press_event(int k) override;
	void close_event() override;

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


bool App::interface()
{
	ImGui::SetCurrentContext(context_);
	ImGui::GetIO().FontGlobalScale = interface_scaling_;

	bool inr = false;

	ImGui::Begin("Control Window",nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0,0});

	inr |= ImGui::RadioButton("Left Viewer", &current_view_, 0);
	ImGui::SameLine();
	inr |= ImGui::RadioButton("Right Viewer", &current_view_, 1);
	inr |= ImGui::Checkbox("BB", &view()->bb_rendering_);
	inr |= ImGui::Checkbox("Phong/Flat", &view()->phong_rendering_);
	inr |= ImGui::Checkbox("Vertices", &view()->vertices_rendering_);
	inr |= ImGui::Checkbox("Normals", &view()->normal_rendering_);
	inr |= ImGui::Checkbox("Edges", &view()->edge_rendering_);
	inr |= ImGui::Checkbox("BB", &view()->bb_rendering_);

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

//	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

	ImGui::End();

	if (inr)
	{
		view()->ask_update();
	}
	return inr;
}

void App::key_press_event(int32 k)
{
	switch(k)
	{
		case int('S'):
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
	ImGUIApp::key_press_event(k);
}

Viewer::Viewer() :
	map_(),
	render_(nullptr),
	vbo_pos_(nullptr),
	vbo_norm_(nullptr),
	vbo_color_(nullptr),
	vbo_sphere_sz_(nullptr),
	drawer_(nullptr),
	drawer_rend_(nullptr),
	phong_rendering_(true),
	vertices_rendering_(false),
	edge_rendering_(false),
	normal_rendering_(false),
	bb_rendering_(true)
{}

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


void Viewer::close_event()
{
	render_.reset();
	vbo_pos_.reset();
	vbo_norm_.reset();
	vbo_color_.reset();
	vbo_sphere_sz_.reset();
	drawer_.reset();
	drawer_rend_.reset();
	fbo_.reset();
	GL::ShaderProgram::clean_all();
}

void Viewer::key_press_event(int k)
{
	switch(k)
	{
		case int('P'):
			phong_rendering_ = true;
			break;
		case int('F'):
			phong_rendering_ = false;
			break;
		case int('N'):
			normal_rendering_ = !normal_rendering_;
			break;
		case int('E'):
			edge_rendering_ = !edge_rendering_;
			break;
		case int('V'):
			vertices_rendering_ = !vertices_rendering_;
			break;
		case int('B'):
			bb_rendering_ = !bb_rendering_;
			break;
		case int('R'):
			cam_.reset();
			break;
		case int('C'):
			cam_.center_scene();
			break;
		default:
			break;
	}
	ImGUIViewer::key_press_event(k);
}


void Viewer::import(const std::string& surface_mesh)
{
	using Vertex = typename mesh_traits<Map2>::Vertex;

	cgogn::io::import_OFF<Vec3>(map_,surface_mesh);

	vertex_position_ = get_attribute<Vec3, Vertex>(map_,"position");

	vertex_normal_ = get_attribute<Vec3, Vertex>(map_,"normal");

	if (vertex_normal_ == nullptr)
	{
		vertex_normal_ = add_attribute<Vec3, Vertex>(map_,"normal");
		cgogn::geometry::compute_normal<Vec3>(map_, vertex_position_, vertex_normal_);
	}
}


void Viewer::init()
{
	update_bb();
	Vec3 diagonal = bb_max_ - bb_min_;
	set_scene_radius(diagonal.norm() / 2);
	Vec3 center = (bb_max_ + bb_min_ )/2 ;
	set_scene_center(center);

	// create and fill VBO for positions
	vbo_pos_ = cgogn::make_unique<GL::VBO>(3);
	GL::update_vbo(	vertex_position_.get(), vbo_pos_.get());

	// create and fill VBO for normals
	vbo_norm_ = cgogn::make_unique<GL::VBO>(3);
	GL::update_vbo(vertex_normal_.get(), vbo_norm_.get());

	// fill a color vbo with abs of normals
	vbo_color_ = cgogn::make_unique<GL::VBO>(3);
	GL::update_vbo(vertex_normal_.get(), vbo_color_.get(), [] (const Vec3& n) -> std::array<float,3>
	{
		return {float(std::abs(n[0])), float(std::abs(n[1])), float(std::abs(n[2]))};
	});
	// fill a sphere size vbo
	vbo_sphere_sz_ = cgogn::make_unique<GL::VBO>(1);
	GL::update_vbo(vertex_normal_.get(), vbo_sphere_sz_.get(), [&] (const Vec3& n) -> float
	{
		return diagonal.norm()/1000.0*(1 + 2*std::abs(n[2]));
	});

	// map rendering object (primitive creation & sending to GPU)
	render_ = cgogn::make_unique<GL::MeshRender>();
	render_->init_primitives(map_, GL::POINTS);
	render_->init_primitives(map_, GL::LINES);
	render_->init_primitives(map_, GL::TRIANGLES);

	// generation of one parameter set (for this shader) : vbo + uniforms
	param_point_sprite_ = GL::ShaderPointSpriteColorSize::generate_param();
	param_point_sprite_->set_vbos(vbo_pos_.get(), vbo_color_.get(), vbo_sphere_sz_.get());
//	param_point_sprite_->size_=0.01f;

	param_edge_ = GL::ShaderBoldLine::generate_param();
	param_edge_->set_vbos(vbo_pos_.get());
	param_edge_->color_ =  GL::GLColor(1,1,1,1);
	param_edge_->width_= 2.5f;

	param_flat_ = GL::ShaderFlat::generate_param();
	param_flat_->set_vbos(vbo_pos_.get());
	param_flat_->front_color_ =  GL::GLColor(0,0.8f,0,1);
	param_flat_->back_color_ =  GL::GLColor(0,0,0.8f,1);
	param_flat_->ambiant_color_ =  GL::GLColor(0.1f,0.1f,0.1f,1);

	param_normal_ = GL::ShaderVectorPerVertex::generate_param();
	param_normal_->set_vbos(vbo_pos_.get(), vbo_norm_.get());
	param_normal_->color_ =  GL::GLColor(0.8f,0.8f,0.8f,1);
	param_normal_->length_ = diagonal.norm()/50;

	param_phong_ = GL::ShaderPhong::generate_param();
	param_phong_->front_color_ =  GL::GLColor(0,0.8f,0,1);
	param_phong_->back_color_ =  GL::GLColor(0,0,0.8f,1);
	param_phong_->ambiant_color_ =  GL::GLColor(0.1f,0.1f,0.1f,1);
	param_phong_->set_vbos(vbo_pos_.get(), vbo_norm_.get());//, vbo_color_.get());

	param_frame_ = GL::ShaderFrame2d::generate_param();
	param_frame_->sz_ = 9.0f;

	// drawer for simple old-school g1 rendering
	drawer_ = cgogn::make_unique<GL::DisplayListDrawer>();
	drawer_rend_= drawer_->generate_renderer();
	drawer_->new_list();
	drawer_->line_width_aa(2.0);
	drawer_->begin(GL_LINE_LOOP);
		drawer_->color3f(1.0,1.0,1.0);
		drawer_->vertex3f(bb_min_[0],bb_min_[1],bb_min_[2]);
		drawer_->vertex3f(bb_max_[0],bb_min_[1],bb_min_[2]);
		drawer_->vertex3f(bb_max_[0],bb_max_[1],bb_min_[2]);
		drawer_->vertex3f(bb_min_[0],bb_max_[1],bb_min_[2]);
		drawer_->vertex3f(bb_min_[0],bb_max_[1],bb_max_[2]);
		drawer_->vertex3f(bb_max_[0],bb_max_[1],bb_max_[2]);
		drawer_->vertex3f(bb_max_[0],bb_min_[1],bb_max_[2]);
		drawer_->vertex3f(bb_min_[0],bb_min_[1],bb_max_[2]);
	drawer_->end();
	drawer_->begin(GL_LINES);
	drawer_->color3f(1.0,1.0,1.0);
		drawer_->vertex3f(bb_min_[0],bb_min_[1],bb_min_[2]);
		drawer_->vertex3f(bb_min_[0],bb_max_[1],bb_min_[2]);
		drawer_->vertex3f(bb_min_[0],bb_min_[1],bb_max_[2]);
		drawer_->vertex3f(bb_min_[0],bb_max_[1],bb_max_[2]);
		drawer_->vertex3f(bb_max_[0],bb_min_[1],bb_min_[2]);
		drawer_->vertex3f(bb_max_[0],bb_min_[1],bb_max_[2]);
		drawer_->vertex3f(bb_max_[0],bb_max_[1],bb_min_[2]);
		drawer_->vertex3f(bb_max_[0],bb_max_[1],bb_max_[2]);
	drawer_->end();
	drawer_->end_list();
}


void Viewer::draw()
{
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.25f,0.25f,0.29f,1);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	GL::GLMat4 proj = get_projection_matrix();
	GL::GLMat4 view = get_modelview_matrix();

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 4.0f);

	if (phong_rendering_)
	{
		param_phong_->bind(proj,view);
		render_->draw(GL::TRIANGLES);
		param_phong_->release();
	}
	else
	{
		param_flat_->bind(proj,view);
		render_->draw(GL::TRIANGLES);
		param_flat_->release();
	}

	glDisable(GL_POLYGON_OFFSET_FILL);

	if (vertices_rendering_)
	{
		param_point_sprite_->bind(proj,view);
		render_->draw(GL::POINTS);
		param_point_sprite_->release();
	}

	if (edge_rendering_)
	{
		param_edge_->bind(proj,view);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		render_->draw(GL::LINES);
		glDisable(GL_BLEND);
		param_edge_->release();
	}

	if (normal_rendering_)
	{
		param_normal_->bind(proj,view);
		render_->draw(GL::POINTS);
		param_normal_->release();
	}

	if (bb_rendering_)
	{
		drawer_rend_->draw(proj,view);
	}

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

	std::string surface_mesh2 = std::string(DEFAULT_MESH_PATH) + std::string("off/horse.off");



	App app;
	gl3wInit();
	Viewer view;
	view.import(surface_mesh);
	app.add_view(&view);
	Viewer view2;
	view2.import(surface_mesh2);
	app.add_view(&view2);
	return app.launch();
}

