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
#include <chrono>
#include <QApplication>
#include <QMatrix4x4>
#include <QKeyEvent>

#include <QOGLViewer/qoglviewer.h>

#include <cgogn/core/cmap/cmap2.h>
//#include <cgogn/core/cmap/cmap2_tri.h>
//#include <cgogn/core/cmap/cmap2_quad.h>

#include <cgogn/core/utils/masks.h>


#include <cgogn/io/map_import.h>

#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/geometry/algos/normal.h>

#include <cgogn/rendering/map_render.h>
#include <cgogn/rendering/shaders/shader_simple_color.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering/shaders/vbo.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/geometry/algos/ear_triangulation.h>

#include <cgogn/rendering/drawer.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)

using namespace cgogn::numerics;

using Map2 = cgogn::CMap2;


using Vec3 = Eigen::Vector3d;

template <typename T>
using VertexAttribute = Map2::VertexAttribute<T>;


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

	Map2 map_;
	VertexAttribute<Vec3> vertex_position_;
	VertexAttribute<Vec3> vertex_pos2_;


	cgogn::geometry::AABB<Vec3> bb_;

	std::unique_ptr<cgogn::rendering::MapRender> render_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_pos_;
	std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;

//	bool flat_rendering_;

	std::mutex mut_update_;
	std::condition_variable cv_update_;
	std::atomic_bool need_update_;

	std::future<void> future_;

};


//
// IMPLEMENTATION
//


void Viewer::import(const std::string& surface_mesh)
{
	cgogn::io::import_surface<Vec3>(map_, surface_mesh);

	vertex_position_ = map_.template get_attribute<Vec3, Map2::Vertex>("position");
	if (!vertex_position_.is_valid())
	{
		cgogn_log_error("Viewer::import") << "Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_pos2_ = map_.template add_attribute<Vec3, Map2::Vertex>("position2");

	map_.foreach_cell([&](Map2::Vertex v)
	{
		vertex_pos2_[v] = vertex_position_[v];
	});


	cgogn::external_thread_pool()->set_nb_workers(1);


	cgogn::geometry::compute_AABB(vertex_position_, bb_);
	setSceneRadius(cgogn::geometry::diagonal(bb_).norm() / 2.0);
	Vec3 center = cgogn::geometry::center(bb_);
	setSceneCenter(qoglviewer::Vec(center[0], center[1], center[2]));
	showEntireScene();
}

Viewer::~Viewer()
{}

void Viewer::closeEvent(QCloseEvent*)
{
	cv_update_.notify_all();
	if (future_.valid())
		future_.wait();
	render_.reset();
	vbo_pos_.reset();
	cgogn::rendering::ShaderProgram::clean_all();

}

Viewer::Viewer() :
	map_(),
	vertex_position_(),
	bb_(),
	render_(nullptr),
	vbo_pos_(nullptr),
	need_update_(false)
{

}

void Viewer::keyPressEvent(QKeyEvent *ev)
{
	switch (ev->key())
	{
	case Qt::Key_R:
		map_.foreach_cell([&](Map2::Vertex v)
		{
			vertex_position_[v] = vertex_pos2_[v];
		});
		cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());
		break;

	case Qt::Key_A:
	{
		if (future_.valid())
		{
			std::future_status status = future_.wait_for(std::chrono::nanoseconds::min());
			if (status == std::future_status::timeout)
				break;
		}
		cgogn_log_info("Asyncrone") << "let's go 20s";

		future_ = cgogn::launch_thread([&]() -> void
		{
			auto start = std::chrono::system_clock::now();
			decltype(start) end;

			VertexAttribute<Vec3> vertex_normal = map_.template get_attribute<Vec3, Map2::Vertex>("normal");
			if (!vertex_normal.is_valid())
				vertex_normal = map_.template add_attribute<Vec3, Map2::Vertex>("normal");
			do
			{
				cgogn::geometry::compute_normal(map_, vertex_position_, vertex_normal);
				for (int i = 0; i < 100; ++i)
				{
					map_.foreach_cell([&](Map2::Vertex v)
					{
						Vec3& P = vertex_position_[v];
						Vec3 N = vertex_normal[v];
						P += 0.0002*N;
					});
					// tag vbo need update
					need_update_ = true;
					// ask for update
					update();
					// wait vbo update finished (we can write in attribute)
					std::unique_lock<std::mutex> lk(mut_update_);
					cv_update_.wait(lk);
				}
				for (int i = 0; i < 100; ++i)
				{
					map_.foreach_cell([&](Map2::Vertex v)
					{
						Vec3& P = vertex_position_[v];
						Vec3 N = vertex_normal[v];
						P -= 0.0002*N;
					});
					// tag vbo need update
					need_update_ = true;
					// ask for update
					update();
					// lock & wait vbo update finished (we can write in attribute)
					std::unique_lock<std::mutex> lk(mut_update_);
					cv_update_.wait(lk);
				}
				end = std::chrono::system_clock::now();
			} while (std::chrono::duration<float>(end - start).count() < 20);
			cgogn_log_info("Asyncrone") << "finished";
		});
	}
	break;



	default:
		break;
	}
	// enable QGLViewer keys
	QOGLViewer::keyPressEvent(ev);

}

void Viewer::draw()
{
	QMatrix4x4 proj;
	QMatrix4x4 view;
	camera()->getProjectionMatrix(proj);
	camera()->getModelViewMatrix(view);

	// VBO need to be sync with attribute ?
	if (need_update_)
	{
		// lock
		std::unique_lock<std::mutex> lk(mut_update_);
		// update
		cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());
		// tag
		need_update_ = false;
		// ask compute thread to continue
		cv_update_.notify_all();
	}

	param_flat_->bind(proj, view);
	render_->draw(cgogn::rendering::TRIANGLES);
	param_flat_->release();

}

void Viewer::init()
{
	glClearColor(0.1f, 0.1f, 0.3f, 0.0f);

	// create and fill VBO for positions
	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	cgogn::rendering::update_vbo(vertex_position_, vbo_pos_.get());

	// create and fill VBO for normals
// map rendering object (primitive creation & sending to GPU)
	render_ = cgogn::make_unique<cgogn::rendering::MapRender>();
	render_->init_primitives(map_, cgogn::rendering::TRIANGLES, &vertex_position_);


	param_flat_ = cgogn::rendering::ShaderFlat::generate_param();
	param_flat_->set_position_vbo(vbo_pos_.get());
	param_flat_->front_color_ = QColor(0, 200, 0);
	param_flat_->back_color_ = QColor(0, 0, 200);
	param_flat_->ambiant_color_ = QColor(5, 5, 5);

}

int main(int argc, char** argv)
{
	std::string surface_mesh;
	if (argc < 2)
	{
		cgogn_log_info("simple_viewer") << "USAGE: " << argv[0] << " [filename]";
		surface_mesh = std::string(DEFAULT_MESH_PATH) + std::string("off/horse.off");
		cgogn_log_info("simple_viewer") << "Using default mesh \"" << surface_mesh << "\".";
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
