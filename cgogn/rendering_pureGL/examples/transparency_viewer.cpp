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
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/rendering/transparency_drawer.h>
#include <cgogn/rendering/wall_paper.h>


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

	virtual void keyPressEvent(QKeyEvent *);
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
	std::unique_ptr<cgogn::rendering::ShaderPointSprite::Param> param_point_sprite_;

	std::shared_ptr<cgogn::rendering::WallPaper> wp_;
	std::unique_ptr<cgogn::rendering::WallPaper::Renderer> wp_rend_;

	std::chrono::time_point<std::chrono::system_clock> start_fps_;
	int nb_fps_;

	bool draw_points_;
	int mesh_transparency_;
	bool lighted_;
	bool bfc_;
	bool phong_rendered_;
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
	setSceneRadius(cgogn::geometry::diagonal(bb_).norm()/2.0);
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
	param_point_sprite_.reset();
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
	param_point_sprite_(nullptr),
	wp_(nullptr),
	wp_rend_(nullptr),
	draw_points_(true),
	mesh_transparency_(130),
	lighted_(true),
	bfc_(false),
	phong_rendered_(false)
{}

void ViewerTransparency::keyPressEvent(QKeyEvent *ev)
{
	switch (ev->key())
	{
		case Qt::Key_P:
			draw_points_ = !draw_points_;
			break;
		case Qt::Key_Plus:
			std::cout <<"mesh_transparency_ " << mesh_transparency_<< std::endl;
			if (mesh_transparency_<254) mesh_transparency_++;
			transp_drawer_->set_front_color(QColor(0,250,0,mesh_transparency_));
			transp_drawer_->set_back_color(QColor(0,250,0,mesh_transparency_));
			break;
		case Qt::Key_Minus:
			std::cout <<"mesh_transparency_ " << mesh_transparency_<< std::endl;
			if (mesh_transparency_>0) mesh_transparency_--;
			transp_drawer_->set_front_color(QColor(0,250,0,mesh_transparency_));
			transp_drawer_->set_back_color(QColor(0,250,0,mesh_transparency_));
			break;
		case Qt::Key_L:
			lighted_ = !lighted_;
			transp_drawer_->set_lighted(lighted_);
			break;
		case Qt::Key_C:
			bfc_ = !bfc_;
			transp_drawer_->set_back_face_culling(bfc_);
			break;
		case Qt::Key_R:
			phong_rendered_ = !phong_rendered_;
			break;
		default:
			break;
	}
	// enable QGLViewer keys
	QOGLViewer::keyPressEvent(ev);
	//update drawing
	update();
}

void ViewerTransparency::draw()
{
	QMatrix4x4 proj;
	QMatrix4x4 view;
	camera()->getProjectionMatrix(proj);
	camera()->getModelViewMatrix(view);

	wp_rend_->draw();

	if (draw_points_)
	{
		param_point_sprite_->bind(proj,view);
		render_->draw(cgogn::rendering::POINTS);
		param_point_sprite_->release();
	}

	if (phong_rendered_)
		transp_drawer_->draw_phong(proj,view, [&] { render_->draw(cgogn::rendering::TRIANGLES); });
	else
		transp_drawer_->draw_flat(proj,view, [&] { render_->draw(cgogn::rendering::TRIANGLES); });

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

	param_point_sprite_ = cgogn::rendering::ShaderPointSprite::generate_param();
	param_point_sprite_->set_position_vbo(vbo_pos_.get());
	param_point_sprite_->size_ = cgogn::geometry::diagonal(bb_).norm()/1000;
	param_point_sprite_->color_ = QColor(200,200,0);

	transp_drawer_ = cgogn::make_unique<cgogn::rendering::SurfaceTransparencyDrawer>();
	transp_drawer_->set_position_vbo(vbo_pos_.get());
	transp_drawer_->set_normal_vbo(vbo_norm_.get());
	std::cout <<"mesh_transparency_ " << mesh_transparency_<< std::endl;
	transp_drawer_->set_front_color(QColor(0,250,0,mesh_transparency_));
	transp_drawer_->set_back_color(QColor(0,250,0,mesh_transparency_));

	wp_ = std::make_shared<cgogn::rendering::WallPaper>(QImage(QString(DEFAULT_MESH_PATH) + QString("../images/cgogn2.png")));
	wp_rend_ = wp_->generate_renderer();

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
		surface_mesh = std::string(DEFAULT_MESH_PATH) + std::string("off/aneurysm_3D.off");
		cgogn_log_info("simple_viewer") << "Using default mesh \"" << surface_mesh << "\".";
	}
	else
		surface_mesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	ViewerTransparency viewer;
	viewer.setWindowTitle("transparency_viewer");
	viewer.import(surface_mesh);
	viewer.show();

	// Run main loop.
	return application.exec();
}
