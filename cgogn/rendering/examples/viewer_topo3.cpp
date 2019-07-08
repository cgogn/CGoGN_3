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
#include <cgogn/rendering/imgui_viewer.h>
#include<GLFW/glfw3.h>

//#include <cgogn/core/utils/logger.h>
//#include <cgogn/core/cmap/cmap3_tetra.h>
//#include <cgogn/core/cmap/cmap3_hexa.h>
#include <cgogn/core/cmap/cmap3.h>
#include <cgogn/io/map_import.h>
#include <cgogn/rendering/map_render.h>
#include <cgogn/rendering/drawer.h>
#include <cgogn/rendering/volume_drawer.h>
#include <cgogn/rendering/topo_drawer.h>
#include <cgogn/rendering/frame_manipulator.h>
#include <cgogn/rendering/vbo_update.h>
#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/modeling/tiling/hexa_volume.h>


#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using namespace cgogn::numerics;
//using Map3 = cgogn::CMap3Tetra;
//using Map3 = cgogn::CMap3Hexa;
using Map3 = cgogn::CMap3;
using Vec3 = Eigen::Vector3d;
//using Vec3 = cgogn::geometry::Vec_T<std::array<float64,3>>;

using namespace cgogn;

template <typename T>
using VertexAttribute = Map3::VertexAttribute<T>;

namespace GL= cgogn::rendering;


class App;

class Viewer : public GL::ImGUIViewer
{
	friend class App;
public:

	using TopoDrawer = GL::TopoDrawer;
	using VolumeDrawer = GL::VolumeDrawer;
	using DisplayListDrawer = GL::DisplayListDrawer;

	Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	void draw() override;
	void init() override;
	void key_press_event(int k) override;
	void key_release_event(int k) override;
	void mouse_press_event(int32 button, float64 x, float64 y) override;
	void mouse_release_event(int32 button, float64 x, float64 y) override;
	void mouse_move_event(float64 dx, float64 dy) override;
	void close_event() override;
	void import(const std::string& surface_mesh);
	virtual ~Viewer() override;

private:

	void rayClick(float64 x, float64 y, GL::GLVec3d& P, GL::GLVec3d& Q);

	void plane_clip_from_frame();

	Map3 map_;
	VertexAttribute<Vec3> vertex_position_;

	cgogn::geometry::AABB<Vec3> bb_;

	std::unique_ptr<GL::VBO> vbo_pos_;

	std::unique_ptr<TopoDrawer> topo_drawer_;
	std::unique_ptr<TopoDrawer::Renderer> topo_drawer_rend_;

	std::unique_ptr<VolumeDrawer> volume_drawer_;
	std::unique_ptr<VolumeDrawer::Renderer> volume_drawer_rend_;

	std::unique_ptr<DisplayListDrawer> drawer_;
	std::unique_ptr<DisplayListDrawer::Renderer> drawer_rend_;

	std::unique_ptr<GL::FrameManipulator> frame_manip_;

	bool vol_rendering_;
	bool edge_rendering_;
	bool topo_drawering_;

	float32 expl_;

	GL::GLVec4 plane_clipping1_;
	float32 plane_thickness_;
	bool thick_plane_mode_;

	GL::GLMat4 proj_matrix_;
	GL::GLMat4 view_matrix_;

};


class App: public GL::ImGUIApp
{
public:
	App() {}
//	inline Viewer* view(int i) { return static_cast<Viewer*>(viewers_[i]); }
//	void interface() override;
};

//
// IMPLEMENTATION
//

void Viewer::rayClick(float64 x, float64 y, GL::GLVec3d &P, GL::GLVec3d &Q)
{
	float64 xgl = x/width() *2.0 - 1.0;
	float64 ygl = (1.0-y/height()) *2.0 - 1.0;

	GL::GLMat4d im =(proj_matrix_.cast<double>()*view_matrix_.cast<double>()).inverse();
	GL::GLVec4d P4 = im*GL::GLVec4d(xgl,ygl,0.0,1.0);
	GL::GLVec4d Q4 = im*GL::GLVec4d(xgl,ygl,1.0,1.0);
	P = GL::GLVec3d(P4.x()/P4.w(),P4.y()/P4.w(),P4.z()/P4.w());
	Q = GL::GLVec3d(Q4.x()/Q4.w(),Q4.y()/Q4.w(),Q4.z()/Q4.w());
}


bool menger(uint32 l, uint32 i, uint32 j, uint32 k)
{
	if (l == 1)
		return true;
	l /= 3;
	if ((i / l == 1) + (j / l == 1) + (k / l == 1) >= 2)
		return false;
	return menger(l, i % l, j % l, k % l);
}

void Viewer::import(const std::string& volumeMesh)
{
	if (volumeMesh=="")
	{
		vertex_position_ = map_.template add_attribute<Vec3, Map3::Vertex>("position");

		const uint32 N = 10;
		cgogn::modeling::TilingHexa tile1(map_, 2*N,2*N,2*N);
		tile1.embedded_grid3D([&](uint32 i, uint32 j, uint32 k)
		{
			i -= N;
			j -= N;
			k -= N;
			auto r = i * i + j * j + k * k;
			return ((r > N*N/4) && (r < N*N)) ||
					(r < 2*2);
		});
		tile1.update_positions([&](cgogn::CMap3::Vertex v, double r, double s, double t)
		{
			vertex_position_[v] = Vec3(r-1.2 ,s ,t*1.5);
		});

		const uint32 NB = 27;
		cgogn::modeling::TilingHexa tile2(map_,NB,NB,NB);
		tile2.embedded_grid3D( [&](uint32 i, uint32 j, uint32 k)
		{
			return menger(NB, i, j, k);
		});
		tile2.update_positions([&](cgogn::CMap3::Vertex v, double r, double s, double t) {
			vertex_position_[v] = Vec3(std::pow(1+2*r,2)/9,std::pow(1+2*s,2)/9,std::pow(1+2*t,2)/9);
		});

		cgogn::modeling::TilingHexa tile3(map_,3,3,50);
		tile3.grid_topo( [](uint32,uint32,uint32){return true;});
		tile3.self_sew_z();
		tile3.close_grid();
		tile3.embed();
		tile3.update_positions([&](cgogn::CMap3::Vertex v, double r, double s, double t) {
			vertex_position_[v] = Vec3((3+r)/2, 0.5+(1+s)/4*std::cos(t*6.283),0.5+(1+s)/4*std::sin(t*6.283));
		});
	}
	else
	{
		cgogn::io::import_volume<Vec3>(map_, volumeMesh);
		vertex_position_ = map_.template get_attribute<Vec3, Map3::Vertex>("position");
		if (!vertex_position_.is_valid())
		{
			cgogn_log_error("Viewer::import") << "Missing attribute position. Aborting.";
			std::exit(EXIT_FAILURE);
		}
	}

	if (!map_.check_map_integrity())
	{
		cgogn_log_error("Viewer::import") << "Integrity of map not respected. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	cgogn::geometry::compute_AABB(vertex_position_, bb_);

	set_scene_radius(cgogn::geometry::diagonal(bb_).norm()/2.0);
	Vec3 center = cgogn::geometry::center(bb_);
	set_scene_center(center);

	map_.check_map_integrity();
}

Viewer::~Viewer()
{}

void Viewer::close_event()
{
	vbo_pos_.reset();
	topo_drawer_.reset();
	topo_drawer_rend_.reset();
	volume_drawer_.reset();
	volume_drawer_rend_.reset();
	drawer_.reset();
	drawer_rend_.reset();
	GL::ShaderProgram::clean_all();
}

Viewer::Viewer() :
	map_(),
	vertex_position_(),
	bb_(),
	vbo_pos_(nullptr),
	topo_drawer_(nullptr),
	topo_drawer_rend_(nullptr),
	volume_drawer_(nullptr),
	volume_drawer_rend_(nullptr),
	drawer_(nullptr),
	drawer_rend_(nullptr),
	vol_rendering_(true),
	edge_rendering_(true),
	topo_drawering_(true),
	expl_(0.8f),
	plane_clipping1_(0,0,0,0),
	thick_plane_mode_(false)
{}

void Viewer::key_press_event(int key)
{
//	if ((ev->modifiers() & Qt::ShiftModifier) && (ev->modifiers() & Qt::ControlModifier))
//		setCursor(Qt::CrossCursor);

	switch (key)
	{
		case int32('V'):
			vol_rendering_ = !vol_rendering_;
			break;
		case int32('E'):
			edge_rendering_ = !edge_rendering_;
			break;
		case int32('T'):
			topo_drawering_ = !topo_drawering_;
			break;
		case GLFW_KEY_KP_MULTIPLY:
			expl_ = 1.0f;
			volume_drawer_rend_->set_explode_volume(expl_);
			topo_drawer_->set_explode_volume(expl_);
			topo_drawer_->update(map_,vertex_position_);
			break;
		case GLFW_KEY_KP_ADD:
			expl_ += 0.05f;
			volume_drawer_rend_->set_explode_volume(expl_);
			topo_drawer_->set_explode_volume(expl_);
			topo_drawer_->update(map_,vertex_position_);
			break;
		case GLFW_KEY_KP_SUBTRACT:
			expl_ -= 0.05f;
			volume_drawer_rend_->set_explode_volume(expl_);
			topo_drawer_->set_explode_volume(expl_);
			topo_drawer_->update(map_,vertex_position_);
			break;
		case int32('X'):
			frame_manip_->rotate(GL::FrameManipulator::Xr, 0.1507f);
			break;
		case int32('P'):
			if (control_pressed())
			{
				thick_plane_mode_ = !thick_plane_mode_;
			}
			else if (thick_plane_mode_)
			{
				if (shift_pressed())
					plane_thickness_ += cgogn::geometry::diagonal(bb_).norm()/200;
				else
				{
					if (plane_thickness_>= cgogn::geometry::diagonal(bb_).norm()/200)
						plane_thickness_ -= cgogn::geometry::diagonal(bb_).norm()/200;
				}
			}
			if (thick_plane_mode_)
			{
				volume_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
				topo_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
			}
			else
			{
				volume_drawer_rend_->set_clipping_plane(plane_clipping1_);
				topo_drawer_rend_->set_clipping_plane(plane_clipping1_);
			}
			break;
		default:
			break;
	}
	// enable QGLViewer keys

}

void Viewer::key_release_event(int32 k)
{
//	QOGLViewer::keyReleaseEvent(ev);
//	unsetCursor();
}

void Viewer::mouse_press_event(int32 button, float64 x, float64 y)
{
	GL::GLVec3d P;
	GL::GLVec3d Q;
	rayClick(x,y, P, Q);
	Vec3 A(P[0], P[1], P[2]);
	Vec3 B(Q[0], Q[1], Q[2]);


	if (control_pressed() && !shift_pressed())
		frame_manip_->pick(x,y,P,Q);

	if (shift_pressed() && !control_pressed())
	{
		drawer_->new_list();
		std::vector<Map3::Volume> selected;
		cgogn::geometry::picking(map_, vertex_position_, A, B, selected);
		cgogn_log_info("Viewer") << "Selected volumes: " << selected.size();
		if (!selected.empty())
		{
			drawer_->line_width(2.0);
			drawer_->begin(GL_LINES);
			// closest vol in red
			drawer_->color3f(1.0, 0.0, 0.0);
			GL::add_to_drawer(map_, selected[0], vertex_position_, drawer_.get());
			// others in yellow
			drawer_->color3f(1.0, 1.0, 0.0);
			for (uint32 i = 1u; i < selected.size(); ++i)
				GL::add_to_drawer(map_, selected[i], vertex_position_, drawer_.get());
			drawer_->end();
		}
		drawer_->line_width(4.0);
		drawer_->begin(GL_LINES);
		drawer_->color3f(1.0, 0.0, 1.0);
		drawer_->vertex3fv(A);
		drawer_->vertex3fv(B);
		drawer_->end();

		drawer_->end_list();
		ask_update();
	}

	if (shift_pressed() && control_pressed())
	{
		cgogn::Dart da;
		if (thick_plane_mode_)
			da = topo_drawer_->pick(A,B,plane_clipping1_, plane_thickness_);
		else
			da = topo_drawer_->pick(A,B,plane_clipping1_);
		if (!da.is_nil())
		{
			topo_drawer_->update_color(da, Vec3(1.0,0.0,0.0));
		}
		ask_update();
	}

	if (shift_pressed() && !control_pressed())
		GL::ImGUIViewer::mouse_press_event(button,x,y);
}

void Viewer::mouse_release_event(int32 button, float64 x, float64 y)
{
	if (control_pressed())
		frame_manip_->release();

	GL::ImGUIViewer::mouse_release_event(button,x,y);
	ask_update();
}

void Viewer::mouse_move_event(float64 dx, float64 dy)
{
	if (control_pressed())
	{
		bool local_manip = mouse_buttons_ & 2;
		frame_manip_->drag(local_manip, last_mouse_x_, last_mouse_y_);

		plane_clip_from_frame();

		if (thick_plane_mode_)
		{
			volume_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
			topo_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
		}
		else
		{
			volume_drawer_rend_->set_clipping_plane(plane_clipping1_);
			topo_drawer_rend_->set_clipping_plane(plane_clipping1_);
		}
	}

	GL::ImGUIViewer::mouse_move_event(dx,dy);
	ask_update();
}

void Viewer::draw()
{
	glViewport(vp_x_,vp_y_, vp_w_, vp_h_);
	glEnable(GL_DEPTH_TEST);

	proj_matrix_ = get_projection_matrix();
	view_matrix_ = get_modelview_matrix();

	if (vol_rendering_)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0f, 1.0f);
		volume_drawer_rend_->draw_faces(proj_matrix_,view_matrix_);
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	if (edge_rendering_)
		volume_drawer_rend_->draw_edges(proj_matrix_,view_matrix_);

	if (topo_drawering_)
		topo_drawer_rend_->draw(proj_matrix_,view_matrix_);

	drawer_rend_->draw(proj_matrix_, view_matrix_);

	frame_manip_->draw(true,true,proj_matrix_, view_matrix_);
}

void Viewer::init()
{
	glClearColor(0.1f,0.1f,0.3f,0.0f);

	vbo_pos_ = cgogn::make_unique<GL::VBO>(3);
	GL::update_vbo(vertex_position_, vbo_pos_.get());

	topo_drawer_ =  cgogn::make_unique<GL::TopoDrawer>();
	topo_drawer_rend_ = topo_drawer_->generate_renderer();
	topo_drawer_->set_explode_volume(expl_);
	topo_drawer_->update(map_,vertex_position_);

	volume_drawer_ = cgogn::make_unique<GL::VolumeDrawer>();
	volume_drawer_->update_face(map_,vertex_position_);
	volume_drawer_->update_edge(map_,vertex_position_);

	volume_drawer_rend_ = volume_drawer_->generate_renderer();
	volume_drawer_rend_->set_explode_volume(expl_);

	drawer_ = cgogn::make_unique<GL::DisplayListDrawer>();
	drawer_rend_ = drawer_->generate_renderer();


	frame_manip_ = cgogn::make_unique<GL::FrameManipulator>();
	frame_manip_->set_size(cgogn::geometry::diagonal(bb_).norm()/4);
	frame_manip_->set_position(bb_.max());
	frame_manip_->z_plane_param(GL::GLColor(0.8f,0.8f,0.8f,1.0f),-1.5f,-1.5f, 2.0f);

	plane_thickness_ = cgogn::geometry::diagonal(bb_).norm()/20;

	plane_clip_from_frame();

	if (thick_plane_mode_)
	{
		volume_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
		topo_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
	}
	else
	{
		volume_drawer_rend_->set_clipping_plane(plane_clipping1_);
		topo_drawer_rend_->set_clipping_plane(plane_clipping1_);
	}
}

void Viewer::plane_clip_from_frame()
{
	Vec3 position;
	Vec3 axis_z;
	frame_manip_->get_position(position);
	frame_manip_->get_axis(GL::FrameManipulator::Zt,axis_z);
	float32 d = -(position.dot(axis_z));
	plane_clipping1_ = GL::GLVec4(axis_z[0],axis_z[1],axis_z[2],d);
}

int main(int argc, char** argv)
{
	std::string volume_mesh;
	if (argc < 2)
	{
		cgogn_log_debug("viewer_topo3") << "USAGE: " << argv[0] << " [filename]";
		volume_mesh = "";//std::string(DEFAULT_MESH_PATH) + std::string("vtk/nine_hexas.vtu");
		cgogn_log_debug("viewer_topo3") << "Using procedural.";
	}
	else
		volume_mesh = std::string(argv[1]);

	App app;
	app.set_window_title("viewer_topo3_IMGUI");
	app.set_window_size(800,800);
	gl3wInit();
	Viewer view;
	view.import(volume_mesh);
	app.add_view(&view);
	app.launch();

	return 0;
}
