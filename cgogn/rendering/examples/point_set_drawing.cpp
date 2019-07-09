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

#include <QOGLViewer/qoglviewer.h>

#include <cgogn/core/cmap/cmap0.h>

#include <cgogn/io/map_import.h>

#include <cgogn/geometry/algos/bounding_box.h>

#include <cgogn/rendering/map_render.h>
#include <cgogn/rendering/shaders/vbo.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_round_point.h>

#include <cgogn/rendering/drawer.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using Map0 = cgogn::CMap0;


using Vec3 = Eigen::Vector3d;
//using Vec3 = cgogn::geometry::Vec_T<std::array<double,3>>;

template <typename T>
using VertexAttribute = Map0::VertexAttribute<T>;

class Viewer : public QOGLViewer
{
public:

	Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	virtual void draw();
	virtual void init();

	virtual void keyPressEvent(QKeyEvent *);
	void import(const std::string& surface_mesh);
	virtual ~Viewer();
	virtual void closeEvent(QCloseEvent *e);

private:

	Map0 map_;
	VertexAttribute<Vec3> vertex_position_;

	cgogn::geometry::AABB<Vec3> bb_;

	std::unique_ptr<cgogn::rendering::MapRender> render_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_pos_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_color_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_sphere_sz_;
	std::unique_ptr<cgogn::rendering::ShaderPointSpriteColorSize::Param> param_point_sprite_;

	std::unique_ptr<cgogn::rendering::DisplayListDrawer> drawer_;
	std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> drawer_rend_;




	bool vertices_rendering_;
	bool bb_rendering_;
};

//
// IMPLEMENTATION
//


void Viewer::import(const std::string& surface_mesh)
{
	cgogn::io::import_point_set<Vec3>(map_, surface_mesh);

	vertex_position_ = map_.template get_attribute<Vec3, Map0::Vertex>("position");
	if (!vertex_position_.is_valid())
	{
		cgogn_log_error("Viewer::import") << "Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	cgogn::geometry::compute_AABB(vertex_position_, bb_);
	setSceneRadius(cgogn::geometry::diagonal(bb_).norm()/2.0);
	Vec3 center = cgogn::geometry::center(bb_);
	setSceneCenter(qoglviewer::Vec(center[0], center[1], center[2]));
	showEntireScene();
}

Viewer::~Viewer()
{}

void Viewer::closeEvent(QCloseEvent*)
{
	render_.reset();
	vbo_pos_.reset();
	vbo_color_.reset();
	vbo_sphere_sz_.reset();
	param_point_sprite_.reset();
	drawer_.reset();
	drawer_rend_.reset();

	cgogn::rendering::ShaderProgram::clean_all();

}

Viewer::Viewer() :
	map_(),
	vertex_position_(),
	bb_(),
	render_(nullptr),
	vbo_pos_(nullptr),
	vbo_color_(nullptr),
	vbo_sphere_sz_(nullptr),
	param_point_sprite_(nullptr),
	drawer_(nullptr),
	drawer_rend_(nullptr),
	vertices_rendering_(true),
	bb_rendering_(true)
{}

void Viewer::keyPressEvent(QKeyEvent *ev)
{
	switch (ev->key())
	{
		case Qt::Key_V:
			vertices_rendering_ = !vertices_rendering_;
			break;
		case Qt::Key_B:
			bb_rendering_ = !bb_rendering_;
			break;
		default:
			break;
	}
	// enable QGLViewer keys
	QOGLViewer::keyPressEvent(ev);
	//update drawing
	update();
}

void Viewer::draw()
{
	QMatrix4x4 proj;
	QMatrix4x4 view;
	camera()->getProjectionMatrix(proj);
	camera()->getModelViewMatrix(view);

	if (vertices_rendering_)
	{
		param_point_sprite_->bind(proj,view);
		render_->draw(cgogn::rendering::POINTS);
		param_point_sprite_->release();
	}

	if (bb_rendering_)
		drawer_rend_->draw(proj,view);
}

void Viewer::init()
{
	glClearColor(0.1f,0.1f,0.3f,0.0f);

	// create and fill VBO for positions
	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());

	// fill a color vbo with abs of normals
	vbo_color_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_, vbo_color_.get(), [] (const Vec3& n) -> std::array<float,3>
	{
		return {float(1.0), float(1.0), float(0.0)};
	});

	// fill a sphere size vbo
	vbo_sphere_sz_ = cgogn::make_unique<cgogn::rendering::VBO>(1);
	cgogn::rendering::update_vbo(vertex_position_, vbo_sphere_sz_.get(), [&] (const Vec3& n) -> float
	{
		return cgogn::geometry::diagonal(bb_).norm()/1000.0;// * (1.0 + 2.0*std::abs(n[2]));
	});

	render_ = cgogn::make_unique<cgogn::rendering::MapRender>();
	render_->init_primitives(map_, cgogn::rendering::POINTS);
	param_point_sprite_ = cgogn::rendering::ShaderPointSpriteColorSize::generate_param();
	// set vbo param (see param::set_vbo signature)
	param_point_sprite_->set_all_vbos(vbo_pos_.get(), vbo_color_.get(), vbo_sphere_sz_.get());

	// drawer for simple old-school g1 rendering
	drawer_ = cgogn::make_unique<cgogn::rendering::DisplayListDrawer>();
	drawer_rend_= drawer_->generate_renderer();
	drawer_->new_list();
	drawer_->line_width_aa(2.0);
	drawer_->begin(GL_LINE_LOOP);
		drawer_->color3f(1.0,1.0,1.0);
		drawer_->vertex3f(bb_.min()[0],bb_.min()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.min()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.max()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.min()[0],bb_.max()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.min()[0],bb_.max()[1],bb_.max()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.max()[1],bb_.max()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.min()[1],bb_.max()[2]);
		drawer_->vertex3f(bb_.min()[0],bb_.min()[1],bb_.max()[2]);
	drawer_->end();
	drawer_->begin(GL_LINES);
	drawer_->color3f(1.0,1.0,1.0);
		drawer_->vertex3f(bb_.min()[0],bb_.min()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.min()[0],bb_.max()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.min()[0],bb_.min()[1],bb_.max()[2]);
		drawer_->vertex3f(bb_.min()[0],bb_.max()[1],bb_.max()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.min()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.min()[1],bb_.max()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.max()[1],bb_.min()[2]);
		drawer_->vertex3f(bb_.max()[0],bb_.max()[1],bb_.max()[2]);
	drawer_->end();
	drawer_->end_list();
}

int main(int argc, char** argv)
{
	std::string surface_mesh;
	if (argc < 2)
	{
		cgogn_log_info("simple_viewer") << "USAGE: " << argv[0] << " [filename]";
		return 1;
	}
	else
		surface_mesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	Viewer viewer;
	viewer.setWindowTitle("simple_viewer");
	viewer.import(surface_mesh);
	viewer.show();

	// Run main loop.
	return application.exec();
}
