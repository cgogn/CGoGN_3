/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
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

//import CGoGn
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>

#include <cgogn/io/utils.h>
#include <cgogn/io/point/point_import.h>
#include <cgogn/io/volume/volume_import.h>


//import CGAL
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/read_off_points.h>

#include <vector>
#include <random>
#include <chrono>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<uint32_t, K>;
using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
using Point = Delaunay::Point;
using Kernel =  CGAL::Simple_cartesian<double>;
using Cgal_Surface_mesh = CGAL::Surface_mesh<Kernel::Point_3> ;

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Mesh_Point = cgogn::CMap0;
template <typename T> 
using Attribute_Point = typename cgogn::mesh_traits<Mesh_Point>::Attribute<T>;
using Vertex_Point = typename cgogn::mesh_traits<Mesh_Point>::Vertex;

using Mesh_Surface = cgogn::CMap2;
template <typename T>
using Attribute_Surface = typename cgogn::mesh_traits<Mesh_Surface>::Attribute<T>;
using Vertex_Surface = typename cgogn::mesh_traits<Mesh_Surface>::Vertex;

using Mesh_Volumn = cgogn::CMap3;
template <typename T>
using Attribute_Volumn = typename cgogn::mesh_traits<Mesh_Volumn>::Attribute<T>;
using namespace cgogn::numerics;

using Vertex_Volumn = typename cgogn::mesh_traits<Mesh_Volumn>::Vertex;
using Volume = typename cgogn::mesh_traits<Mesh_Volumn>::Volume;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;



//#define PERF_TEST

void test_delaunay(Mesh_Point* mp, Mesh_Volumn* mv, Cgal_Surface_mesh csm)
{
	cgogn::io::VolumeImportData volume_data;
	cgogn::io::PointImportData point_data;

	std::vector<Kernel::Point_3> mesh_samples;
	CGAL::Polygon_mesh_processing::sample_triangle_mesh(csm, 
														std::back_inserter(mesh_samples),
														CGAL::parameters::use_grid_sampling(true).grid_spacing(5));

	//Add vertex into the volume data to construct the delauney tredrahedron
	for (auto s : mesh_samples)
	{
		volume_data.vertex_position_.emplace_back(s[0], s[1], s[2]);
		point_data.vertex_position_.emplace_back(s[0], s[1], s[2]);
	}
	import_point_data(*mp, point_data);
	uint32 nb_vertices = volume_data.vertex_position_.size();

	auto start_timer = std::chrono::high_resolution_clock::now();

	std::vector<unsigned> indices;
	indices.reserve(nb_vertices);
	for (unsigned int i = 0; i < nb_vertices; ++i) 
		indices.push_back(i);

	//Construct delauney tredrahedron using CGAL
	std::vector<Point>* points = reinterpret_cast<std::vector<Point>*>(&volume_data.vertex_position_);

	Delaunay tri(boost::make_zip_iterator(boost::make_tuple(points->begin(), indices.begin())),
				 boost::make_zip_iterator(boost::make_tuple(points->end(), indices.end())));

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
	std::cout << "CGAL delaunay for " << nb_vertices << " points in " << elapsed_seconds.count() << std::endl;

	start_timer = std::chrono::high_resolution_clock::now();

	uint32 nb_volumes = tri.number_of_finite_cells();
	volume_data.reserve(nb_vertices, nb_volumes);

	for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
	{
		volume_data.volumes_types_.push_back(cgogn::io::VolumeType::Tetra);
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(0)->info());
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(1)->info());
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(2)->info());
		volume_data.volumes_vertex_indices_.push_back(cit->vertex(3)->info());
	}

	import_volume_data(*mv, volume_data);
	

	end_timer = std::chrono::high_resolution_clock::now();
	elapsed_seconds = end_timer - start_timer;
	std::cout << "transfer to CGoGN in " << elapsed_seconds.count() << std::endl;
}

int main(int argc, char** argv)
{
	std::string filename;
	if (argc > 1)
		filename = std::string(argv[1]);
	
	//Load surface mesh using CGAL;
	Cgal_Surface_mesh surface_mesh;
	
	std::ifstream input(filename);
	if (!input || !(input >> surface_mesh))
	{
		std::cerr << "Error: input OFF file could not be read" << std::endl;
		return EXIT_FAILURE;
	}

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Delaunay volume viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh_Volumn> mpv(app);
	cgogn::ui::MeshProvider<Mesh_Surface> mps(app);
	cgogn::ui::MeshProvider<Mesh_Point> mpp(app);

	cgogn::ui::VolumeRender<Mesh_Volumn> vr(app);
	cgogn::ui::SurfaceRender<Mesh_Surface> sr(app);
	cgogn::ui::PointCloudRender<Mesh_Point> pcr(app);
	
	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mpv);
	v1->link_module(&vr);
	v1->link_module(&mps);
	v1->link_module(&sr);
	v1->link_module(&mpp);
	v1->link_module(&pcr);


 	Mesh_Volumn* mv = new Mesh_Volumn{};
	Mesh_Point* mpoint = new Mesh_Point{};
	test_delaunay(mpoint, mv, surface_mesh);
	mpv.register_mesh(mv, "delaunay");
	mpp.register_mesh(mpoint, "sample point");


	std::shared_ptr<Attribute_Volumn<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex_Volumn>(*mv, "position");
	mpv.set_mesh_bb_vertex_position(*mv, vertex_position);
	vr.set_vertex_position(*v1, *mv, vertex_position);

	std::shared_ptr<Attribute_Point<Vec3>> vertex_position_2 =
		cgogn::get_attribute<Vec3, Vertex_Point>(*mv, "position");
	pcr.set_vertex_position(*v1, *mpoint, vertex_position_2);

	//Load surface mesh using CGoGn
	Mesh_Surface* original_m = mps.load_surface_from_file(filename);
	if (!original_m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute_Surface<Vec3>> vertex_position_3 = cgogn::get_attribute<Vec3, Vertex_Surface>(*original_m, "position");
	std::shared_ptr<Attribute_Surface<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex_Surface>(*original_m, "normal");

	sr.set_vertex_position(*v1, *original_m, vertex_position_3);
	sr.set_vertex_normal(*v1, *original_m, vertex_normal);

	return app.launch();
}
