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
#include <QMouseEvent>
#include <QVector3D>

#include <QOGLViewer/qoglviewer.h>
#include <QOGLViewer/vec.h>

#include <cgogn/core/cmap/cmap2.h>

#include <cgogn/io/map_import.h>

#include <cgogn/geometry/algos/bounding_box.h>

#include <cgogn/rendering/map_render.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/vbo.h>
#include <cgogn/rendering/drawer.h>

#include <cgogn/geometry/algos/picking.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using namespace cgogn::numerics;

using Map2 = cgogn::CMap2;
//using Vec3 = Eigen::Vector3d;
using Vec3 = cgogn::geometry::Vec_T<std::array<float64,3>>;

template <typename T>
using VertexAttribute = Map2::VertexAttribute<T>;

class Viewer : public QOGLViewer
{
public:

	Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	virtual void draw();
	virtual void init();
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void keyPressEvent(QKeyEvent *);

	void import(const std::string& surfaceMesh);

	virtual ~Viewer();

private:

	void rayClick(QMouseEvent* event, qoglviewer::Vec& P, qoglviewer::Vec& Q);

	QMatrix4x4 proj_;
	QMatrix4x4 view_;

	Map2 map_;

	VertexAttribute<Vec3> vertex_position_;

	cgogn::geometry::AABB<Vec3> bb_;

	std::unique_ptr<cgogn::rendering::MapRender> render_;


	std::unique_ptr<cgogn::rendering::VBO> vbo_pos_;

	std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;

	std::unique_ptr<cgogn::rendering::DisplayListDrawer> drawer_;
	std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> drawer_rend_;

	int32 cell_picking;
};

//
// IMPLEMENTATION
//

void Viewer::import(const std::string& surfaceMesh)
{
	cgogn::io::import_surface<Vec3>(map_, surfaceMesh);

	vertex_position_ = map_.template get_attribute<Vec3, Map2::Vertex>("position");

	cgogn::geometry::compute_AABB(vertex_position_, bb_);

	setSceneRadius(cgogn::geometry::diagonal(bb_).norm()/2.0);
	Vec3 center = cgogn::geometry::center(bb_);
	setSceneCenter(qoglviewer::Vec(center[0], center[1], center[2]));
	showEntireScene();
}

Viewer::~Viewer()
{}

Viewer::Viewer() :
	map_(),
	vertex_position_(),
	bb_(),
	render_(nullptr),
	vbo_pos_(nullptr),
	drawer_(nullptr),
	cell_picking(0)
{}

void Viewer::draw()
{
	camera()->getProjectionMatrix(proj_);
	camera()->getModelViewMatrix(view_);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 1.0f);

	param_flat_->bind(proj_, view_);
	render_->draw(cgogn::rendering::TRIANGLES);
	param_flat_->release();

	glDisable(GL_POLYGON_OFFSET_FILL);

	drawer_rend_->draw(proj_, view_);
}

void Viewer::init()
{
	glClearColor(0.1f, 0.1f, 0.3f, 0.0f);

	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());

	render_ = cgogn::make_unique<cgogn::rendering::MapRender>();
	render_->init_primitives(map_, cgogn::rendering::TRIANGLES, &vertex_position_);

	param_flat_ = cgogn::rendering::ShaderFlat::generate_param();

	param_flat_->set_position_vbo(vbo_pos_.get());
	param_flat_->front_color_ = QColor(0,200,0);
	param_flat_->back_color_ = QColor(200,0,0);
	param_flat_->ambiant_color_ = QColor(5,5,5);

	drawer_ = cgogn::make_unique<cgogn::rendering::DisplayListDrawer>();
	drawer_rend_ = drawer_->generate_renderer();
}

void Viewer::rayClick(QMouseEvent* event, qoglviewer::Vec& P, qoglviewer::Vec& Q)
{
	P = camera()->unprojectedCoordinatesOf(qoglviewer::Vec(event->x(), event->y(), 0.0));
	Q = camera()->unprojectedCoordinatesOf(qoglviewer::Vec(event->x(), event->y(), 1.0));
}

void Viewer::keyPressEvent(QKeyEvent *ev)
{
	switch (ev->key())
	{
		case Qt::Key_0:
			cell_picking = 0;
			break;
		case Qt::Key_1:
			cell_picking = 1;
			break;
		case Qt::Key_2:
			cell_picking = 2;
			break;
		case Qt::Key_3:
			cell_picking = 3;
			break;
	}
	QOGLViewer::keyPressEvent(ev);
}

void Viewer::mousePressEvent(QMouseEvent* event)
{
	if (event->modifiers() & Qt::ShiftModifier)
	{
		qoglviewer::Vec P;
		qoglviewer::Vec Q;
		rayClick(event,P,Q);

		Vec3 A(P[0],P[1],P[2]);
		Vec3 B(Q[0],Q[1],Q[2]);

		drawer_->new_list();
		switch(cell_picking)
		{
			case 0:
			{
				std::vector<Map2::Vertex> selected;
                cgogn::geometry::picking(map_,vertex_position_, A, B, selected);
				cgogn_log_info("picking_viewer") << "Selected vertices: "<< selected.size();
				if (!selected.empty())
				{
					drawer_->point_size_aa(4.0);
					drawer_->begin(GL_POINTS);
					// closest point in red
					drawer_->color3f(1.0,0.0,0.0);
					drawer_->vertex3fv(vertex_position_[selected[0]]);
					// others in yellow
					drawer_->color3f(1.0,1.0,0.0);
					for(uint32 i=1u;i<selected.size();++i)
						drawer_->vertex3fv(vertex_position_[selected[i]]);
					drawer_->end();
				}
			}
			break;
			case 1:
			{
				std::vector<Map2::Edge> selected;
				cgogn::geometry::picking(map_, vertex_position_, A, B, selected);
				cgogn_log_info("picking_viewer") << "Selected edges: "<< selected.size();
				if (!selected.empty())
				{
					drawer_->line_width(2.0);
					drawer_->begin(GL_LINES);
					// closest face in red
					drawer_->color3f(1.0, 0.0, 0.0);
					cgogn::rendering::add_to_drawer(map_, selected[0], vertex_position_, drawer_.get());
					// others in yellow
					drawer_->color3f(1.0, 1.0, 0.0);
					for(uint32 i = 1u; i < selected.size(); ++i)
						cgogn::rendering::add_to_drawer(map_, selected[i], vertex_position_, drawer_.get());
					drawer_->end();
				}
			}
			break;
			case 2:
			{
				std::vector<Map2::Face> selected;
				cgogn::geometry::picking(map_, vertex_position_, A, B, selected);
				cgogn_log_info("picking_viewer") << "Selected faces: "<< selected.size();
				if (!selected.empty())
				{
					drawer_->line_width(2.0);
					drawer_->begin(GL_LINES);
					// closest face in red
					drawer_->color3f(1.0, 0.0, 0.0);
					cgogn::rendering::add_to_drawer(map_, selected[0], vertex_position_, drawer_.get());
					// others in yellow
					drawer_->color3f(1.0, 1.0, 0.0);
					for(uint32 i = 1u; i < selected.size(); ++i)
						cgogn::rendering::add_to_drawer(map_, selected[i], vertex_position_, drawer_.get());
					drawer_->end();
				}
			}
			break;
			case 3:
			{
				std::vector<Map2::Volume> selected;
				cgogn::geometry::picking(map_, vertex_position_, A, B, selected);
				cgogn_log_info("picking_viewer") << "Selected volumes: "<< selected.size();
				if (!selected.empty())
				{
					drawer_->line_width(2.0);
					drawer_->begin(GL_LINES);
					// closest face in red
					drawer_->color3f(1.0, 0.0, 0.0);
					cgogn::rendering::add_to_drawer(map_, selected[0], vertex_position_, drawer_.get());
					// others in yellow
					drawer_->color3f(1.0, 1.0, 0.0);
					for(uint32 i = 1u; i < selected.size(); ++i)
						cgogn::rendering::add_to_drawer(map_, selected[i], vertex_position_, drawer_.get());
					drawer_->end();
				}
			}
			break;
		}
		drawer_->line_width(4.0);
		drawer_->begin(GL_LINES);
		drawer_->color3f(1.0, 0.0, 1.0);
		drawer_->vertex3fv(A);
		drawer_->vertex3fv(B);
		drawer_->end();
		drawer_->end_list();
	}

	QOGLViewer::mousePressEvent(event);
}

int main(int argc, char** argv)
{
	std::string surfaceMesh;
	if (argc < 2)
	{
		cgogn_log_info("picking_viewer") << "USAGE: " << argv[0] << " [filename]";
		surfaceMesh = std::string(DEFAULT_MESH_PATH) + std::string("off/aneurysm_3D.off");
		cgogn_log_info("picking_viewer") << "Using default mesh \"" << surfaceMesh << "\".";
	}
	else
		surfaceMesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	Viewer viewer;
	viewer.setWindowTitle("picking_viewer");
	viewer.import(surfaceMesh);
	viewer.show();
	viewer.resize(800, 600);

	// Run main loop.
	return application.exec();
}
