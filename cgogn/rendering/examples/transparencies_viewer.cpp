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

#include <QApplication>
#include <QMatrix4x4>
#include <QKeyEvent>
#include <chrono>

#include <QOGLViewer/qoglviewer.h>

#include <cgogn/core/cmap/cmap2.h>

#include <cgogn/io/map_import.h>

#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/geometry/algos/normal.h>

#include <cgogn/rendering/map_render.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/transparency_drawer.h>

#include <cgogn/rendering/drawer.h>


using namespace cgogn::numerics;

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Map2 = cgogn::CMap2;

using Vec3 = Eigen::Vector3d;

template <typename T>
using VertexAttribute = Map2::VertexAttribute<T>;



class ViewerTransparency : public QOGLViewer
{
public:

	ViewerTransparency();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ViewerTransparency);

	virtual void draw();
	virtual void init();

	void import(const std::string& surface_mesh);
	virtual ~ViewerTransparency();
	virtual void closeEvent(QCloseEvent *e);
	virtual void resizeGL(int width, int height);

private:

	Map2 map_;
	VertexAttribute<Vec3> vertex_position_;
	VertexAttribute<Vec3> vertex_normal_;
	cgogn::geometry::AABB<Vec3> bb_;
	std::unique_ptr<cgogn::rendering::MapRender> render_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_pos_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_norm_;

	std::unique_ptr<cgogn::rendering::SurfaceTransparencyDrawer> transp_drawer_;

	std::unique_ptr<cgogn::rendering::ShaderFlatTransp::Param> tr_flat_param_;
	std::unique_ptr<cgogn::rendering::ShaderPhongTransp::Param> tr_phong_param_;

	std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;

	std::unique_ptr<cgogn::rendering::DisplayListDrawer> drawer_;
	std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> drawer_rend_;


	std::chrono::time_point<std::chrono::system_clock> start_fps_;
	int nb_fps_;

};


//
// IMPLEMENTATION
//


void ViewerTransparency::import(const std::string& surface_mesh)
{
	cgogn::io::import_surface<Vec3>(map_, surface_mesh);

	vertex_position_ = map_.get_attribute<Vec3, Map2::Vertex::ORBIT>("position");
	if (!vertex_position_.is_valid())
	{
		cgogn_log_error("ViewerTransparency::import") << "Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_normal_ = map_.template get_attribute<Vec3, Map2::Vertex>("normal");
	if (!vertex_normal_.is_valid())
	{
		vertex_normal_ = map_.template add_attribute<Vec3, Map2::Vertex>("normal");
		cgogn::geometry::compute_normal(map_, vertex_position_, vertex_normal_);
	}


	cgogn::geometry::compute_AABB(vertex_position_, bb_);
	setSceneRadius(cgogn::geometry::diagonal(bb_).norm());
	Vec3 center = cgogn::geometry::center(bb_);
	setSceneCenter(qoglviewer::Vec(center[0], center[1], center[2]));
	showEntireScene();
}


ViewerTransparency::~ViewerTransparency()
{}


void ViewerTransparency::closeEvent(QCloseEvent*)
{
	render_.reset();
	vbo_pos_.reset();
	vbo_norm_.reset();
	transp_drawer_.reset();
	cgogn::rendering::ShaderProgram::clean_all();
}

ViewerTransparency::ViewerTransparency() :
	map_(),
	vertex_position_(),
	vertex_normal_(),
	bb_(),
	render_(nullptr),
	vbo_pos_(nullptr),
	vbo_norm_(nullptr),
	transp_drawer_(nullptr),
	tr_flat_param_(nullptr),
	tr_phong_param_(nullptr),
	param_flat_(nullptr)
{}


void ViewerTransparency::draw()
{
	QMatrix4x4 proj;
	QMatrix4x4 view;
	camera()->getProjectionMatrix(proj);
	camera()->getModelViewMatrix(view);

	// draw opaque first
	param_flat_->bind(proj,view);
	render_->draw(cgogn::rendering::TRIANGLES);
	param_flat_->release();

	drawer_rend_->draw(proj,view);

	// the the transparents objects.

	QMatrix4x4 tr1;
	tr1.translate(-0.25*cgogn::geometry::diagonal(bb_).norm(),0,0);
	QMatrix4x4 tr2;
	tr2.translate(0.25*cgogn::geometry::diagonal(bb_).norm(),0,0);

	transp_drawer_->draw( [&] ()
	{
		tr_flat_param_->set_alpha(150);
		tr_flat_param_->bind(proj,view*tr1);
		render_->draw(cgogn::rendering::TRIANGLES);
		tr_flat_param_->release();

		tr_phong_param_->set_alpha(170);
		tr_phong_param_->bind(proj,view*tr2);
		render_->draw(cgogn::rendering::TRIANGLES);
		tr_phong_param_->release();

		tr_flat_param_->set_alpha(100);
		tr_flat_param_->bind(proj,view*tr1*tr1);
		render_->draw(cgogn::rendering::TRIANGLES);
		tr_flat_param_->release();

		tr_phong_param_->set_alpha(90);
		tr_phong_param_->bind(proj,view*tr2*tr2);
		render_->draw(cgogn::rendering::TRIANGLES);
		tr_phong_param_->release();
	});


	nb_fps_++;
	std::chrono::duration<float64> elapsed_seconds = std::chrono::system_clock::now() - start_fps_;
	if (elapsed_seconds.count()>= 5)
	{
		cgogn_log_info("fps") << double(nb_fps_)/elapsed_seconds.count();
		nb_fps_ = 0;
		start_fps_ = std::chrono::system_clock::now();
	}

}

void ViewerTransparency::init()
{
	glClearColor(0.1f,0.1f,0.1f,0.0f);

	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());

	vbo_norm_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_normal_, vbo_norm_.get());

	render_ = cgogn::make_unique<cgogn::rendering::MapRender>();
	render_->init_primitives(map_, cgogn::rendering::TRIANGLES, &vertex_position_);
	render_->init_primitives(map_, cgogn::rendering::POINTS, &vertex_position_);


	param_flat_ = cgogn::rendering::ShaderFlat::generate_param();
	param_flat_->set_position_vbo(vbo_pos_.get());
	param_flat_->front_color_ = QColor(0,50,200);
	param_flat_->back_color_ = QColor(0,50,200);

	transp_drawer_ = cgogn::make_unique<cgogn::rendering::SurfaceTransparencyDrawer>();
//	transp_drawer_->set_max_nb_layers(16);

	tr_flat_param_ = cgogn::rendering::ShaderFlatTransp::generate_param();
	tr_flat_param_->set_position_vbo(vbo_pos_.get());
	tr_flat_param_->front_color_ = QColor(0,250,0,100);
	tr_flat_param_->back_color_ = QColor(0,250,0,100);

	tr_phong_param_ = cgogn::rendering::ShaderPhongTransp::generate_param();
	tr_phong_param_->set_position_vbo(vbo_pos_.get());
	tr_phong_param_->set_normal_vbo(vbo_norm_.get());
	tr_phong_param_->front_color_ = QColor(250,0,0,120);
	tr_phong_param_->back_color_ = QColor(250,0,0,120);


	drawer_ = cgogn::make_unique<cgogn::rendering::DisplayListDrawer>();
	drawer_rend_= drawer_->generate_renderer();
	drawer_->new_list();

	drawer_->ball_size(cgogn::geometry::diagonal(bb_).norm()/50);
	drawer_->begin(GL_POINTS);
		drawer_->color3f(1,1,0);
		Vec3 P = cgogn::geometry::center(bb_);
		drawer_->vertex3fv(P);
		P[0] -= 0.25*cgogn::geometry::diagonal(bb_).norm();
		drawer_->vertex3fv(P);
		P[0] -= 0.25*cgogn::geometry::diagonal(bb_).norm();
		drawer_->vertex3fv(P);
		P[0] += 0.75*cgogn::geometry::diagonal(bb_).norm();
		drawer_->vertex3fv(P);
		P[0] += 0.25*cgogn::geometry::diagonal(bb_).norm();
		drawer_->vertex3fv(P);
	drawer_->end();
	drawer_->end_list();


	start_fps_ = std::chrono::system_clock::now();
	nb_fps_ = 0;
}

void ViewerTransparency::resizeGL(int w ,int h)
{
	transp_drawer_->resize(this->devicePixelRatio()*w,this->devicePixelRatio()*h);
	QOGLViewer::resizeGL(w,h);
}


int main(int argc, char** argv)
{
	std::string surface_mesh;
	if (argc < 2)
	{
		cgogn_log_info("transparency_viewer") << "USAGE: " << argv[0] << " [filename]";
		surface_mesh = std::string(DEFAULT_MESH_PATH) + std::string("off/horse.off");
		cgogn_log_info("simple_viewer") << "Using default mesh \"" << surface_mesh << "\".";
	}
	else
		surface_mesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	ViewerTransparency viewer;
	viewer.setWindowTitle("transparencies_viewer");
	viewer.import(surface_mesh);
	viewer.show();

	// Run main loop.
	return application.exec();
}
