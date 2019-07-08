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

#include <cgogn/core/cmap/cmap2.h>

#include <cgogn/io/map_import.h>

#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/geometry/algos/normal.h>

#include <cgogn/rendering/map_render.h>
//#include <cgogn/rendering/shaders/shader_simple_color.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/vbo.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_TEST_MESHES_PATH)


//USING_CGOGN_NUMERICS;
using namespace cgogn::numerics;

using Map2 = cgogn::CMap2;
using Vec3 = Eigen::Vector3d;
//using Vec3 = cgogn::geometry::Vec_T<std::array<float64,3>>;

template <typename T>
using VertexAttribute = Map2::VertexAttribute<T>;

template <typename T>
using FaceAttribute = Map2::FaceAttribute<T>;


class Viewer : public QOGLViewer
{
public:

	Viewer();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(Viewer);

	virtual void draw();
	virtual void init();

	virtual void keyPressEvent(QKeyEvent *);
	void import(const std::string& surfaceMesh);
	virtual ~Viewer();
	virtual void closeEvent(QCloseEvent *e);

private:

	Map2 map_;
	VertexAttribute<Vec3> vertex_position_;
	FaceAttribute<Vec3> face_normal_;

	cgogn::geometry::AABB<Vec3> bb_;

	std::unique_ptr<cgogn::rendering::VBO> vbo_pos_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_norm_;
	std::unique_ptr<cgogn::rendering::VBO> vbo_color_;

	std::unique_ptr<cgogn::rendering::ShaderFlatColor::Param> param_flat_;
	std::unique_ptr<cgogn::rendering::ShaderPhongColor::Param> param_phong_;

	bool phong_rendering_;
	bool flat_rendering_;
};


//
// IMPLEMENTATION
//


void Viewer::import(const std::string& surfaceMesh)
{
	cgogn::io::import_surface<Vec3>(map_, surfaceMesh);

	vertex_position_ = map_.template get_attribute<Vec3, Map2::Vertex>("position");
	face_normal_ = map_.template add_attribute<Vec3, Map2::Face>("normal");

	cgogn::geometry::compute_normal(map_, vertex_position_, face_normal_);
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
	vbo_pos_.reset();
	vbo_norm_.reset();
	vbo_color_.reset();
	cgogn::rendering::ShaderProgram::clean_all();
}

Viewer::Viewer() :
	map_(),
	vertex_position_(),
	face_normal_(),
	bb_(),
	vbo_pos_(nullptr),
	vbo_norm_(nullptr),
	vbo_color_(nullptr),
	phong_rendering_(true),
	flat_rendering_(false)
{}

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


	if (flat_rendering_)
	{
		param_flat_->bind(proj,view);
		glDrawArrays(GL_TRIANGLES,0,vbo_pos_->size());
		param_flat_->release();
	}

	if (phong_rendering_)
	{
		param_phong_->bind(proj,view);
		glDrawArrays(GL_TRIANGLES,0,vbo_pos_->size());
		param_phong_->release();
	}
}

void Viewer::init()
{
	glClearColor(0.1f,0.1f,0.3f,0.0f);


	vbo_pos_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	vbo_norm_ = cgogn::make_unique<cgogn::rendering::VBO>(3);
	vbo_color_ = cgogn::make_unique<cgogn::rendering::VBO>(3);

	// indices of vertices emb (f1_v1,f1_v2,f1_v3, f2_v1,f2_v2,f2_v3, f3_v1...)
	std::vector<uint32> ind_v;
	// indices of faces emb (f1,f1,f1, f2,f2,f2, f3,f3,f3...)
	std::vector<uint32> ind_f;

	// create indices ( need to be done only after topo modifications
	cgogn::rendering::create_indices_vertices_faces(map_,vertex_position_,ind_v,ind_f);

	// generate VBO: positions
	cgogn::rendering::generate_vbo(vertex_position_, ind_v, vbo_pos_.get(), [] (const Vec3& v) -> std::array<float32,3>
	{
		return {float32(v[0]), float32(v[1]), float32(v[2])};
	});

	// generate VBO: normals
	cgogn::rendering::generate_vbo(face_normal_, ind_f, vbo_norm_.get(), [] (const Vec3& n) -> std::array<float32,3>
	{
		return {float32(n[0]), float32(n[1]), float32(n[2])};
	});

	// generate VBO: colors (here abs of normals)
	cgogn::rendering::generate_vbo(face_normal_, ind_f, vbo_color_.get(), [] (const Vec3& n) -> std::array<float32,3>
	{
		return {float32(std::abs(n[0])), float32(std::abs(n[1])), float32(std::abs(n[2]))};
	});

	param_phong_ = cgogn::rendering::ShaderPhongColor::generate_param();
	param_phong_->set_all_vbos(vbo_pos_.get(), vbo_norm_.get(), vbo_color_.get());
	param_phong_->ambiant_color_ = QColor(5, 5, 5);
	param_phong_->double_side_ = true;
	param_phong_->specular_color_ = QColor(255, 255, 255);
	param_phong_->specular_coef_ = 100.0;

	param_flat_ = cgogn::rendering::ShaderFlatColor::generate_param();
	param_flat_->set_all_vbos(vbo_pos_.get(), vbo_color_.get());
	param_flat_->ambiant_color_ = QColor(5, 5, 5);
}

int main(int argc, char** argv)
{
	std::string surfaceMesh;
	if (argc < 2)
	{
		cgogn_log_info("viewer_per_face") << "USAGE: " << argv[0] << " [filename]";
		surfaceMesh = std::string(DEFAULT_MESH_PATH) + std::string("off/aneurysm_3D.off");
		cgogn_log_info("viewer_per_face") << "Using default mesh \"" << surfaceMesh << "\".";
	}
	else
		surfaceMesh = std::string(argv[1]);

	QApplication application(argc, argv);
	qoglviewer::init_ogl_context();

	// Instantiate the viewer.
	Viewer viewer;
	viewer.setWindowTitle("viewer_per_face");
	viewer.import(surfaceMesh);
	viewer.show();

	// Run main loop.
	return application.exec();
}
