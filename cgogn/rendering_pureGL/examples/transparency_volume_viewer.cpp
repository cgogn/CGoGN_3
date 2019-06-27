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

#include <cgogn/core/utils/logger.h>
//#include <cgogn/core/cmap/cmap3_tetra.h>
//#include <cgogn/core/cmap/cmap3_hexa.h>
#include <cgogn/core/cmap/cmap3.h>
#include <cgogn/io/map_import.h>
#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/rendering/shaders/vbo.h>

#include <cgogn/rendering/transparency_drawer.h>
#include <cgogn/rendering/transparency_volume_drawer.h>
#include <cgogn/rendering/frame_manipulator.h>
#include <cgogn/rendering/wall_paper.h>


#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using namespace cgogn::numerics;
//using Map3 = cgogn::CMap3Tetra;
//using Map3 = cgogn::CMap3Hexa;
using Map3 = cgogn::CMap3;
using Vec3 = Eigen::Vector3d;
//using Vec3 = cgogn::geometry::Vec_T<std::array<float64,3>>;

template <typename T>
using VertexAttribute = Map3::VertexAttribute<T>;


class Viewer : public QOGLViewer
{
public:

	using VolumeDrawer = cgogn::rendering::VolumeTransparencyDrawer;

	Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	virtual void draw();
	virtual void init();
	virtual void keyPressEvent(QKeyEvent*);
	virtual void keyReleaseEvent(QKeyEvent*);
	virtual void mousePressEvent(QMouseEvent*);
	virtual void mouseReleaseEvent(QMouseEvent*);
	virtual void mouseMoveEvent(QMouseEvent*);
	virtual void resizeGL(int w ,int h);

	void import(const std::string& volumeMesh);
	virtual ~Viewer();
	virtual void closeEvent(QCloseEvent *e);

private:

	void rayClick(QMouseEvent* event, qoglviewer::Vec& P, qoglviewer::Vec& Q);

	void plane_clip_from_frame();

	Map3 map_;
	VertexAttribute<Vec3> vertex_position_;

	cgogn::geometry::AABB<Vec3> bb_;

	std::unique_ptr<cgogn::rendering::VBO> vbo_pos_;

	std::unique_ptr<cgogn::rendering::SurfaceTransparencyDrawer> transp_drawer_;
	std::unique_ptr<VolumeDrawer> volume_drawer_;
	std::unique_ptr<VolumeDrawer::Renderer> volume_drawer_rend_;

	std::unique_ptr<cgogn::rendering::FrameManipulator> frame_manip_;

	float32 expl_;

	QVector4D plane_clipping1_;
	float32 plane_thickness_;
	bool thick_plane_mode_;
	int mesh_transparency_;
	bool lighted_;
	bool bfc_;

	std::chrono::time_point<std::chrono::system_clock> start_fps_;
	int nb_fps_;

};

//
// IMPLEMENTATION
//

void Viewer::rayClick(QMouseEvent* event, qoglviewer::Vec& P, qoglviewer::Vec& Q)
{
	P = camera()->unprojectedCoordinatesOf(qoglviewer::Vec(event->x(), event->y(), 0.0));
	Q = camera()->unprojectedCoordinatesOf(qoglviewer::Vec(event->x(), event->y(), 1.0));
}

void Viewer::import(const std::string& volumeMesh)
{
	cgogn::io::import_volume<Vec3>(map_, volumeMesh);

	vertex_position_ = map_.template get_attribute<Vec3, Map3::Vertex>("position");
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

	setSceneRadius(cgogn::geometry::diagonal(bb_).norm()/2.0);
	Vec3 center = cgogn::geometry::center(bb_);
	setSceneCenter(qoglviewer::Vec(center[0], center[1], center[2]));
	showEntireScene();

	map_.check_map_integrity();
}

Viewer::~Viewer()
{}

void Viewer::closeEvent(QCloseEvent*)
{
	vbo_pos_.reset();
	volume_drawer_.reset();
	volume_drawer_rend_.reset();
	cgogn::rendering::ShaderProgram::clean_all();
}

Viewer::Viewer() :
	map_(),
	vertex_position_(),
	bb_(),
	vbo_pos_(nullptr),
	volume_drawer_(nullptr),
	volume_drawer_rend_(nullptr),
	frame_manip_(nullptr),
	expl_(0.8f),
	plane_clipping1_(0,0,0,0),
	thick_plane_mode_(false),
	mesh_transparency_(130),
	lighted_(true),
	bfc_(false)
{}

void Viewer::keyPressEvent(QKeyEvent *ev)
{
	if ((ev->modifiers() & Qt::ShiftModifier) && (ev->modifiers() & Qt::ControlModifier))
		setCursor(Qt::CrossCursor);

	switch (ev->key())
	{
		case Qt::Key_Plus:
			expl_ += 0.05f;
			volume_drawer_rend_->set_explode_volume(expl_);
			break;
		case Qt::Key_Minus:
			expl_ -= 0.05f;
			volume_drawer_rend_->set_explode_volume(expl_);
			break;
		case Qt::Key_P:
			if (ev->modifiers() & Qt::ControlModifier)
			{
				thick_plane_mode_ = !thick_plane_mode_;
			}
			else if (thick_plane_mode_)
			{
				if (ev->modifiers() & Qt::ShiftModifier)
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
			}
			else
			{
				volume_drawer_rend_->set_clipping_plane(plane_clipping1_);
			}
			break;


		case Qt::Key_A:
			mesh_transparency_ += 1;
			if (mesh_transparency_ > 254)
				mesh_transparency_ = 254;
			volume_drawer_rend_->set_color(QColor(0,250,0,mesh_transparency_));
			break;
		case Qt::Key_Z:
			mesh_transparency_ -= 1;
			if (mesh_transparency_ < 0)
				mesh_transparency_ = 0;
			volume_drawer_rend_->set_color(QColor(0,250,0,mesh_transparency_));
			break;
		case Qt::Key_L:
			lighted_ = !lighted_;
			volume_drawer_rend_->set_lighted(lighted_);
			break;
		case Qt::Key_C:
			bfc_ = !bfc_;
			volume_drawer_rend_->set_back_face_culling(bfc_);
			break;

//		case Qt::Key_O:
//			volume_drawer_rend_->set_ogl(this);
//			break;
		default:
			break;
	}
	// enable QGLViewer keys
	QOGLViewer::keyPressEvent(ev);
	//update drawing
	update();
}

void Viewer::keyReleaseEvent(QKeyEvent* ev)
{
	QOGLViewer::keyReleaseEvent(ev);
	unsetCursor();
}

void Viewer::mousePressEvent(QMouseEvent* event)
{
	qoglviewer::Vec P;
	qoglviewer::Vec Q;
	rayClick(event, P, Q);
	Vec3 A(P[0], P[1], P[2]);
	Vec3 B(Q[0], Q[1], Q[2]);


	if ((event->modifiers() & Qt::ControlModifier) && !(event->modifiers() & Qt::ShiftModifier))
		frame_manip_->pick(event->x(), event->y(),P,Q);

	QOGLViewer::mousePressEvent(event);
	update();
}

void Viewer::mouseReleaseEvent(QMouseEvent* event)
{
	if (event->modifiers() & Qt::ControlModifier)
		frame_manip_->release();

	QOGLViewer::mouseReleaseEvent(event);
	update();
}

void Viewer::mouseMoveEvent(QMouseEvent* event)
{
	if (event->modifiers() & Qt::ControlModifier)
	{
		bool local_manip = (event->buttons() & Qt::RightButton);
		frame_manip_->drag(local_manip, event->x(), event->y());

		plane_clip_from_frame();

		if (thick_plane_mode_)
		{
			volume_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
		}
		else
		{
			volume_drawer_rend_->set_clipping_plane(plane_clipping1_);
		}
	}

	QOGLViewer::mouseMoveEvent(event);
	update();
}

void Viewer::draw()
{
	QMatrix4x4 proj;
	QMatrix4x4 view;
	camera()->getProjectionMatrix(proj);
	camera()->getModelViewMatrix(view);

	frame_manip_->draw(true,true,proj, view); // draw opaque first

	transp_drawer_->draw([&]() -> void
	{
		volume_drawer_rend_->draw_faces(proj, view);
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


void Viewer::init()
{
	glClearColor(0.1f,0.1f,0.5f,0.0f);

	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());

	volume_drawer_ = cgogn::make_unique<VolumeDrawer>();
	volume_drawer_->update_face(map_,vertex_position_);

	volume_drawer_rend_ = volume_drawer_->generate_renderer();
	volume_drawer_rend_->set_explode_volume(expl_);
	volume_drawer_rend_->set_color(QColor(0,250,0,mesh_transparency_));

	transp_drawer_ = cgogn::make_unique<cgogn::rendering::SurfaceTransparencyDrawer>();
	transp_drawer_->set_max_nb_layers(16);

	frame_manip_ = cgogn::make_unique<cgogn::rendering::FrameManipulator>();
	frame_manip_->set_size(cgogn::geometry::diagonal(bb_).norm()/4);
	frame_manip_->set_position(bb_.max());
	frame_manip_->z_plane_param(QColor(200,200,200),-1.5f,-1.5f, 2.0f);

	plane_thickness_ = cgogn::geometry::diagonal(bb_).norm()/20;

	plane_clip_from_frame();

	if (thick_plane_mode_)
	{
		volume_drawer_rend_->set_thick_clipping_plane(plane_clipping1_,plane_thickness_);
	}
	else
	{
		volume_drawer_rend_->set_clipping_plane(plane_clipping1_);
	}

	start_fps_ = std::chrono::system_clock::now();
	nb_fps_ = 0;
}

void Viewer::resizeGL(int w ,int h)
{
	transp_drawer_->resize(this->devicePixelRatio()*w,this->devicePixelRatio()*h);
	QOGLViewer::resizeGL(w,h);
}


void Viewer::plane_clip_from_frame()
{
	Vec3 position;
	Vec3 axis_z;
	frame_manip_->get_position(position);
	frame_manip_->get_axis(cgogn::rendering::FrameManipulator::Zt,axis_z);
	float32 d = -(position.dot(axis_z));
	plane_clipping1_ = QVector4D(axis_z[0],axis_z[1],axis_z[2],d);
}

int main(int argc, char** argv)
{
	std::string volumeMesh;
	if (argc < 2)
	{
		cgogn_log_debug("transparency_volume_viewer") << "USAGE: " << argv[0] << " [filename]";
		volumeMesh = std::string(DEFAULT_MESH_PATH) + std::string("vtk/nine_hexas.vtu");
		cgogn_log_debug("viewer_topo3") << "Using default mesh \"" << volumeMesh << "\".";
	}
	else
		volumeMesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	Viewer viewer;
	viewer.setWindowTitle("viewer_topo3");
	viewer.import(volumeMesh);
	viewer.show();

	// Run main loop.
	return application.exec();
}
