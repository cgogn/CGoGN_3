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
#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/io/surface_import.h>
#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/rendering_pureGL/drawer.h>
#include <cgogn/rendering_pureGL/map_render.h>
#include <cgogn/rendering_pureGL/topo_drawer.h>
#include <cgogn/rendering_pureGL/shaders/shader_flat.h>
#include <cgogn/rendering_pureGL/shaders/shader_simple_color.h>
#include <cgogn/rendering_pureGL/vbo_update.h>


#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Map2 = cgogn::CMap2;
using Vertex = Map2::Vertex;

using Vec3 = Eigen::Vector3d;
//using Vec3 = cgogn::geometry::Vec_T<std::array<float64,3>>;

using namespace cgogn;

template <typename T>
using VertexAttribute = Map2::VertexAttribute<T>;

namespace GL= cgogn::rendering_pgl;

class App;

class Viewer : public GL::ImGUIViewer
{
	friend class App;
public:

	using MapRender = GL::MapRender;
	using TopoDrawer = GL::TopoDrawer;

	Viewer( GL::ImGUIViewer* share = nullptr);
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	void draw() override;
	void init() override;
	void key_press_event(int k) override;
	void mouse_press_event(int32 button, float64 x, float64 y) override;
	void close_event() override;
	void import(const std::string& surface_mesh);
	virtual ~Viewer() override;

private:
	Map2 map_;
	VertexAttribute<Vec3> vertex_position_;
	cgogn::geometry::AABB<Vec3> bb_;
	std::unique_ptr<MapRender> render_;
	std::unique_ptr<GL::VBO> vbo_pos_;
	std::unique_ptr<GL::ShaderFlat::Param> param_flat_;
	std::unique_ptr<TopoDrawer> topo_drawer_;
	std::unique_ptr<TopoDrawer::Renderer> topo_drawer_rend_;
	bool flat_rendering_;
	bool topo_drawing_;
};

class App: public GL::ImGUIApp
{
	int current_view_;
public:
	App():	current_view_(0) {}
	GL::GLColor clear_color;
	Viewer* view(int i) { return static_cast<Viewer*>(viewers_[i]); }
	void interface() override;
	void key_press_event(int k) override;
};



void App::interface()
{
	ImGui::SetCurrentContext(context_);
//	imgui_make_context_current();
	ImGui::GetIO().FontGlobalScale = interface_scaling_;

	ImGui::Begin("Control Window",nullptr, ImGuiWindowFlags_NoSavedSettings);
	ImGui::SetWindowSize({0,0});
	ImGui::SliderInt("View", &current_view_,0,viewers_.size()-1);
	ImGui::Checkbox("Flat", &view(current_view_)->flat_rendering_);
	ImGui::Checkbox("Topo", &view(current_view_)->topo_drawing_);
	ImGui::Text("Flat parameters");
	ImGui::ColorEdit3("front color##flat",view(current_view_)->param_flat_->front_color_.data(),ImGuiColorEditFlags_NoInputs);
	ImGui::SameLine();
	ImGui::ColorEdit3("back color##flat",view(current_view_)->param_flat_->back_color_.data(),ImGuiColorEditFlags_NoInputs);
	ImGui::Checkbox("single side##flat", &(view(current_view_)->param_flat_->bf_culling_));
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::End();
}

void App::key_press_event(int32 k)
{
	switch(k)
	{
		case int('S'):
			if (focused_->shift_pressed())
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

void Viewer::import(const std::string& surface_mesh)
{
	cgogn::io::import_surface<Vec3>(map_, surface_mesh);

	vertex_position_ = map_.template get_attribute<Vec3, Map2::Vertex>("position");
	if (!vertex_position_.is_valid())
	{
		cgogn_log_error("Viewer::import") << "Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	if (!map_.check_map_integrity())
	{
		cgogn_log_error("Viewer::import") << "Integrity of map not respected. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	cgogn::geometry::compute_AABB(vertex_position_, bb_);
	set_scene_center(cgogn::geometry::center(bb_));
	set_scene_radius(cgogn::geometry::diagonal(bb_).norm()/2.0);
}

Viewer::~Viewer()
{}

void Viewer::close_event()
{
	render_.reset();
	vbo_pos_.reset();
	topo_drawer_.reset();
	topo_drawer_rend_.reset();
	GL::ShaderProgram::clean_all();
}

Viewer::Viewer(GL::ImGUIViewer* share) :
	GL::ImGUIViewer(share),
	map_(),
	vertex_position_(),
	bb_(),
	render_(nullptr),
	vbo_pos_(nullptr),
	topo_drawer_(nullptr),
	topo_drawer_rend_(nullptr),
	flat_rendering_(true),
	topo_drawing_(true)
{}

void Viewer::key_press_event(int k)
{
	switch (k)
	{
		case int('F'):
			flat_rendering_ = !flat_rendering_;
			break;
		case int('T'):
			topo_drawing_ = !topo_drawing_;
			break;
		case int('E'):
			cgogn::io::export_surface(map_, cgogn::io::ExportOptions::create().filename("/tmp/pipo.vtp").position_attribute(Map2::Vertex::ORBIT, "position"));
			break;
		default:
			break;
	}
}

void Viewer::draw()
{
//glClearColor(0.5f,0.5f,0.9f,0);
//	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glViewport(vp_x_,vp_y_, vp_w_, vp_h_);
	glEnable(GL_DEPTH_TEST);

	GL::GLMat4 proj = get_projection_matrix();
	GL::GLMat4 view = get_modelview_matrix();

	if (flat_rendering_)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0f, 1.0f);
		param_flat_->bind(proj,view);
		render_->draw(GL::TRIANGLES);
		param_flat_->release();
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	if (topo_drawing_)
	{
		topo_drawer_rend_->draw(proj,view);
	}
}

void Viewer::init()
{
	glClearColor(0.1f, 0.1f, 0.3f, 0.0f);

	vbo_pos_ = cgogn::make_unique<GL::VBO>(3);
	GL::update_vbo(vertex_position_, vbo_pos_.get());

	render_ = cgogn::make_unique<GL::MapRender>();
	render_->init_primitives(map_, GL::TRIANGLES);

	param_flat_ = GL::ShaderFlat::generate_param();
	param_flat_->set_vbos(vbo_pos_.get());
	param_flat_->front_color_ = GL::GLColor(0,0.7f,0,1);
	param_flat_->back_color_ = GL::GLColor(0,0,0.7f,1);
	param_flat_->ambiant_color_ = GL::GLColor(0.1f,0.1f,0.1f,1);

	topo_drawer_ = cgogn::make_unique<GL::TopoDrawer>();
	topo_drawer_rend_ = topo_drawer_->generate_renderer();
	topo_drawer_->update(map_,vertex_position_);
}

void Viewer::mouse_press_event(int32 button, float64 x, float64 y)
{
	ImGUIViewer::mouse_press_event(button,x,y);
}

int main(int argc, char** argv)
{
	std::string surface_mesh;
	if (argc < 2)
	{
		cgogn_log_info("viewer_topo")<< "USAGE: " << argv[0] << " [filename]";
		surface_mesh = std::string(DEFAULT_MESH_PATH) + std::string("off/aneurysm_3D.off");
		cgogn_log_info("viewer_topo") << "Using default mesh \"" << surface_mesh << "\".";
	}
	else
		surface_mesh = std::string(argv[1]);

	// Instantiate the viewer.
//	Viewer view;
//	gl3wInit();
//	view.import(surface_mesh);
//	view.set_window_title("SimpleViewerIMGUI");
//	view.launch();

	App app;
	app.set_window_title("SimpleViewerIMGUI");
	app.set_window_size(1024,512);
	gl3wInit();
	Viewer view;
	view.import(surface_mesh);
	app.add_view(&view);
	Viewer view2;
	view2.import(surface_mesh);
	app.add_view(&view2);
	Viewer view3;
	view3.import(surface_mesh);
	app.add_view(&view3);
	Viewer view4;
	view4.import(surface_mesh);
	app.add_view(&view4);

	app.launch();

	// Run main loop.
	return 0;
}
