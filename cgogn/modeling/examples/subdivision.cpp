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

#include <Eigen/Dense>

#include <QOGLViewer/qoglviewer.h>

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

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/vbo.h>

#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/io/surface_import.h>

using Map2 = cgogn::CMap2;
using Vertex = Map2::Vertex;
using Edge = Map2::Edge;
using Face = Map2::Face;

using Vec3 = Eigen::Vector3f;
using Scalar = typename cgogn::geometry::vector_traits<Vec3>::Scalar;

template <typename T>
using AttributePtr = typename cgogn::mesh_traits<Map2>::AttributePtr<T>;

class Viewer : public QOGLViewer
{
public:

	Viewer();
	virtual ~Viewer();

	Viewer(const Viewer&) = delete;
	Viewer& operator=(const Viewer&) = delete;

	virtual void closeEvent(QCloseEvent *e);

	void import(const std::string& surface_mesh);
	void update_bb();

	virtual void keyPressEvent(QKeyEvent *);

	virtual void draw();
	virtual void init();

private:

	Map2 map_;
	cgogn::CellFilter<Map2> filtered_map_;

	AttributePtr<Vec3> vertex_position_;
	AttributePtr<Vec3> vertex_normal_;

	Vec3 bb_min_, bb_max_;

	std::unique_ptr<cgogn::rendering::MeshRender> render_;

	std::unique_ptr<cgogn::rendering::VBO> vbo_position_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_normal_;

	std::unique_ptr<cgogn::rendering::ShaderBoldLine::Param> param_edge_;
	std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;
	std::unique_ptr<cgogn::rendering::ShaderVectorPerVertex::Param> param_normal_;
	std::unique_ptr<cgogn::rendering::ShaderPhong::Param> param_phong_;
	std::unique_ptr<cgogn::rendering::ShaderPointSprite::Param> param_point_sprite_;

	bool phong_rendering_;
	bool flat_rendering_;
	bool vertices_rendering_;
	bool edge_rendering_;
	bool normal_rendering_;
};

////////////////////
// IMPLEMENTATION //
////////////////////

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

void Viewer::closeEvent(QCloseEvent*)
{
	render_.reset();
	vbo_position_.reset();
	vbo_normal_.reset();
	cgogn::rendering::ShaderProgram::clean_all();
}

void Viewer::import(const std::string& surface_mesh)
{
	cgogn::io::import_OFF<Vec3>(map_, surface_mesh);

	vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(map_, "position");
	if (!vertex_position_)
	{
		std::cerr << "Viewer::import: Missing attribute position. Aborting.";
		std::exit(EXIT_FAILURE);
	}

	vertex_normal_ = cgogn::add_attribute<Vec3, Vertex>(map_, "normal");
	cgogn::geometry::compute_normal<Vec3>(map_, vertex_position_, vertex_normal_);

	update_bb();

	Vec3 diagonal = bb_max_ - bb_min_;
	setSceneRadius(diagonal.norm() / 2.0f);
	Vec3 center = (bb_max_ + bb_min_) / 2.0f;
	setSceneCenter(qoglviewer::Vec(center[0], center[1], center[2]));
	showEntireScene();
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

void Viewer::keyPressEvent(QKeyEvent *ev)
{
	switch (ev->key())
	{
		case Qt::Key_P:
			phong_rendering_ = true;
			flat_rendering_ = false;
			break;
		case Qt::Key_F:
			flat_rendering_ = true;
			phong_rendering_ = false;
			break;
		case Qt::Key_N:
			normal_rendering_ = !normal_rendering_;
			break;
		case Qt::Key_E:
			edge_rendering_ = !edge_rendering_;
			break;
		case Qt::Key_V:
			vertices_rendering_ = !vertices_rendering_;
			break;
//		case Qt::Key_B:
//			bb_rendering_ = !bb_rendering_;
//			break;
		case Qt::Key_D: {
			cgogn::modeling::decimate<Vec3>(filtered_map_, vertex_position_, cgogn::uint32(0.1 * cgogn::nb_cells<Vertex>(filtered_map_)));
			std::cout << "nbv: " << cgogn::nb_cells<Vertex>(map_) << std::endl;
			cgogn::geometry::compute_normal<Vec3>(map_, vertex_position_, vertex_normal_);
			Scalar mel = cgogn::geometry::mean_edge_length<Vec3>(map_, vertex_position_);
			param_point_sprite_->size_ = mel / 6.0f;
			cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());
			cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());
			render_->init_primitives(map_, cgogn::rendering::POINTS);
			render_->init_primitives(map_, cgogn::rendering::LINES);
			render_->init_primitives(map_, cgogn::rendering::TRIANGLES);
			update_bb();
			Vec3 diagonal = bb_max_ - bb_min_;
			setSceneRadius(diagonal.norm() / 2.0f);
			break;
		}
		case Qt::Key_S: {
			cgogn::modeling::subdivide<Vec3>(map_, vertex_position_);
			std::cout << "nbv: " << cgogn::nb_cells<Vertex>(map_) << std::endl;
			cgogn::geometry::compute_normal<Vec3>(map_, vertex_position_, vertex_normal_);
			Scalar mel = cgogn::geometry::mean_edge_length<Vec3>(map_, vertex_position_);
			param_point_sprite_->size_ = mel / 6.0f;
			cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());
			cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());
			render_->init_primitives(map_, cgogn::rendering::POINTS);
			render_->init_primitives(map_, cgogn::rendering::LINES);
			render_->init_primitives(map_, cgogn::rendering::TRIANGLES);
			update_bb();
			Vec3 diagonal = bb_max_ - bb_min_;
			setSceneRadius(diagonal.norm() / 2.0f);
			break;
		}
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

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 2.0f);
	if (flat_rendering_)
	{
		param_flat_->bind(proj, view);
		render_->draw(cgogn::rendering::TRIANGLES);
		param_flat_->release();
	}

	if (phong_rendering_)
	{
		param_phong_->bind(proj, view);
		render_->draw(cgogn::rendering::TRIANGLES);
		param_phong_->release();
	}
	glDisable(GL_POLYGON_OFFSET_FILL);

	if (vertices_rendering_)
	{
		param_point_sprite_->bind(proj, view);
		render_->draw(cgogn::rendering::POINTS);
		param_point_sprite_->release();
	}

	if (edge_rendering_)
	{
		param_edge_->bind(proj, view);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		render_->draw(cgogn::rendering::LINES);
		glDisable(GL_BLEND);
		param_edge_->release();
	}

	if (normal_rendering_)
	{
		param_normal_->bind(proj, view);
		render_->draw(cgogn::rendering::POINTS);
		param_normal_->release();
	}
}

void Viewer::init()
{
	glClearColor(0.1f,0.1f,0.3f,0.0f);

	vbo_position_ = cgogn::make_unique<cgogn::rendering::VBO>();
	vbo_normal_ = cgogn::make_unique<cgogn::rendering::VBO>();

	cgogn::rendering::update_vbo(vertex_position_.get(), vbo_position_.get());
	cgogn::rendering::update_vbo(vertex_normal_.get(), vbo_normal_.get());

	render_ = cgogn::make_unique<cgogn::rendering::MeshRender>();

	render_->init_primitives(map_, cgogn::rendering::POINTS);
	render_->init_primitives(map_, cgogn::rendering::LINES);
	render_->init_primitives(map_, cgogn::rendering::TRIANGLES);

	param_point_sprite_ = cgogn::rendering::ShaderPointSprite::generate_param();
	param_point_sprite_->set_position_vbo(vbo_position_.get());
	param_point_sprite_->color_ = QColor(255, 0, 0);
	Scalar mel = cgogn::geometry::mean_edge_length<Vec3>(map_, vertex_position_);
	param_point_sprite_->size_ = mel / 6.0f;

	param_edge_ = cgogn::rendering::ShaderBoldLine::generate_param();
	param_edge_->set_position_vbo(vbo_position_.get());
	param_edge_->color_ = QColor(255,255,0);
	param_edge_->width_= 2.5f;

	param_flat_ =  cgogn::rendering::ShaderFlat::generate_param();
	param_flat_->set_position_vbo(vbo_position_.get());
	param_flat_->front_color_ = QColor(0,200,0);
	param_flat_->back_color_ = QColor(0,0,200);
	param_flat_->ambiant_color_ = QColor(5,5,5);

	param_normal_ = cgogn::rendering::ShaderVectorPerVertex::generate_param();
	param_normal_->set_all_vbos(vbo_position_.get(), vbo_normal_.get());
	param_normal_->color_ = QColor(200,0,200);
	param_normal_->length_ = (bb_max_ - bb_min_).norm() / 50.0f;

	param_phong_ = cgogn::rendering::ShaderPhong::generate_param();
	param_phong_->set_all_vbos(vbo_position_.get(), vbo_normal_.get());
}

int main(int argc, char** argv)
{
	std::string surface_mesh;
	if (argc < 2)
	{
		std::cout << "filtering: USAGE: " << argv[0] << " [filename]" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		surface_mesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	Viewer viewer;
	viewer.setWindowTitle("subdivision");
	viewer.import(surface_mesh);
	viewer.show();

	// Run main loop.
	return application.exec();
}
