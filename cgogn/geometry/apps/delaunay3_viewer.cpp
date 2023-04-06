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
#include <cgogn/core/types/cmap/graph.h>
#include <cgogn/core/functions/attributes.h>


//import CGAL
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Object.h>
#include <CGAL/Side_of_triangle_mesh.h>




#include <vector>
#include <random>
#include <chrono>


//Kernel for construct delauney 
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<uint32_t, K>;
using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
using Regular = CGAL::Regular_triangulation_3<K, Tds, CGAL::Fast_location>;
using Point = K::Point_3;
using Weight = K::FT;
using Weight_Point = Regular::Weighted_point;
using Cgal_Surface_mesh = CGAL::Surface_mesh<Point>;
using Point_inside = CGAL::Side_of_triangle_mesh<Cgal_Surface_mesh, K>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Cgal_Surface_mesh>;
using Traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

//Kernel for test
// using Point = Delaunay::Point;
// using Kernel =  CGAL::Simple_cartesian<double>;
// using Point_3 = K::Point_3;
// using Cgal_Surface_mesh = CGAL::Surface_mesh<Point_3> ;
// using Point_inside = CGAL::Side_of_triangle_mesh<Cgal_Surface_mesh,Kernel>;
// using Primitive = CGAL::AABB_face_graph_triangle_primitive<Cgal_Surface_mesh>;
// using Traits =  CGAL::AABB_traits<Kernel, Primitive>;
// using Tree = CGAL::AABB_tree<Traits>;


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
using Graph = cgogn::Graph;


//#define PERF_TEST

bool pointInside(Tree &tree, Point& query)
{
	// Initialize the point-in-polyhedron tester
	Point_inside inside_tester(tree);

	// Determine the side and return true if inside!
	return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
}

void test_delaunay(Mesh_Point* mp, Mesh_Volumn * mv, Cgal_Surface_mesh csm)
{
	cgogn::io::VolumeImportData Delaunay_tri_3_data;
	cgogn::io::PointImportData Voronoi_Vertices_data;
	cgogn::io::VolumeImportData Medial_axis_data;

	std::vector<Point> Delaunay_tri_point;
	std::vector<Weight_Point> Power_point;
	// 1. Compute the Voronoi diagram of the sample points S.
	Tree tree(faces(csm).first, faces(csm).second, csm);
	tree.accelerate_distance_queries();
	// Compute the Bounding box of mesh
	std::array<Point, 8> obb_points;
	Point acc(0,0,0);
	CGAL::oriented_bounding_box(csm, obb_points, CGAL::parameters::use_convex_hull(true));
	for (size_t i = 0; i < obb_points.size(); i++)
	{
		acc += K::Vector_3(obb_points[i].x(), obb_points[i].y(), obb_points[i].z());
	}
	std::array<double, 3> center{acc.x() / 8, acc.y() / 8, acc.z() / 8};
	//Create a large box surrounding object so that the Voronoi vertices are bounded
	double offset = 2.0;
	std::array<double, 3> offset_array = {offset, offset, offset};
	std::array<std::array<double, 3>, 8> cube_corners = {
		{{center[0] - offset_array[0], center[1] - offset_array[1],
		  center[2] - offset_array[2]},
		 {center[0] - offset_array[0], center[1] - offset_array[1],
		  center[2] + offset_array[2]},
		 {center[0] - offset_array[0], center[1] + offset_array[1],
		  center[2] - offset_array[2]},
		 {center[0] - offset_array[0], center[1] + offset_array[1],
		  center[2] + offset_array[2]},
		 {center[0] + offset_array[0], center[1] - offset_array[1],
		  center[2] - offset_array[2]},
		 {center[0] + offset_array[0], center[1] - offset_array[1],
		  center[2] + offset_array[2]},
		 {center[0] + offset_array[0], center[1] + offset_array[1],
		  center[2] - offset_array[2]},
		 {center[0] + offset_array[0], center[1] + offset_array[1],
		  center[2] + offset_array[2]}}};

	// Sampling the mesh surface
	std::vector<Point> mesh_samples;
	CGAL::Polygon_mesh_processing::sample_triangle_mesh(
		csm, std::back_inserter(mesh_samples),
		// CGAL::parameters::use_monte_carlo_sampling(true).number_of_points_per_area_unit(0.1));
		CGAL::parameters::use_grid_sampling(true).grid_spacing(1));

	// 	Add bounding box vertices in the sample points set
	for (auto p : cube_corners)
  	{
		Delaunay_tri_point.emplace_back(p[0], p[1], p[2]);
  	}

	//Add sampled vertices into the volume data to construct the delauney tredrahedron
	for (auto s : mesh_samples)
	{
		
		Delaunay_tri_point.emplace_back(s[0], s[1], s[2]);
	}

	uint32 nb_vertices = Delaunay_tri_point.size();

	auto start_timer = std::chrono::high_resolution_clock::now();

	std::vector<unsigned> indices;
	indices.reserve(nb_vertices);
	for (unsigned int i = 0; i < nb_vertices; ++i) 
		indices.push_back(i);

	//Construct delauney tredrahedron using CGAL
	
	Delaunay tri(boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.begin(), indices.begin())),
				 boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.end(), indices.end())));

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
	std::cout << "CGAL delaunay for " << nb_vertices << " points in " << elapsed_seconds.count() << std::endl;
	
	start_timer = std::chrono::high_resolution_clock::now();
	
	uint32 nb_volumes = tri.number_of_finite_cells();
	Delaunay_tri_3_data.reserve(nb_vertices, nb_volumes);
	Voronoi_Vertices_data.reserve(nb_volumes);

	for (auto p : Delaunay_tri_point)
	{
		Delaunay_tri_3_data.vertex_position_.push_back({p[0], p[1], p[2]});
	}

	for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
	{
		Delaunay_tri_3_data.volumes_types_.push_back(cgogn::io::VolumeType::Tetra);
		Delaunay_tri_3_data.volumes_vertex_indices_.push_back(cit->vertex(0)->info());
		Delaunay_tri_3_data.volumes_vertex_indices_.push_back(cit->vertex(1)->info());
		Delaunay_tri_3_data.volumes_vertex_indices_.push_back(cit->vertex(2)->info());
		Delaunay_tri_3_data.volumes_vertex_indices_.push_back(cit->vertex(3)->info());
	}
	//Find inside poles 
 	for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
 	{
 		double distance = 0.0f;
 		Delaunay::Vertex_handle v = vit;
 		Point farthest_vertex;
 		std::vector<Delaunay::Facet> facets;
 		tri.finite_incident_facets(vit, std::back_inserter(facets));
 		if (facets.size() != 0)
 		{
 			for (auto f = facets.begin(); f != facets.end(); ++f)
 			{
 				K::Object_3 o = tri.dual(*f);
 				if (const K::Segment_3* s = CGAL::object_cast<K::Segment_3>(&o))
 				{
 					auto start = s->start(), source = s->source();
					if (f == facets.begin() && pointInside(tree, start))
					{
 						auto dis = CGAL::squared_distance(start, v->point());
 						if (dis > distance)
 						{
 							distance = dis;
 							farthest_vertex = start;
 						}
 					}
					if (pointInside(tree, source))
					{
						auto dis = CGAL::squared_distance(source, v->point());
						if (dis > distance)
						{
							distance = dis;
							farthest_vertex = source;
						}
					}
	
 				}
 			}
 			Voronoi_Vertices_data.vertex_position_.push_back(
 				{farthest_vertex[0], farthest_vertex[1], farthest_vertex[2]});
			Power_point.push_back(Weight_Point(farthest_vertex, distance));
 		}
 	}
	//Build weighted delaunay tredrahedron

	cgogn::io::VolumeImportData Power_shape_data;
	Regular reg;
	reg.insert(Power_point.begin(), Power_point.end());
	for (auto p : Power_point)
	{
		Power_shape_data.vertex_position_.push_back({p.x(), p.y(), p.z()});
	}

	for (auto cit = reg.finite_cells_begin(); cit != reg.finite_cells_end(); ++cit)
	{
		Power_shape_data.volumes_types_.push_back(cgogn::io::VolumeType::Tetra);
		Power_shape_data.volumes_vertex_indices_.push_back(cit->vertex(0)->info());
		Power_shape_data.volumes_vertex_indices_.push_back(cit->vertex(1)->info());
		Power_shape_data.volumes_vertex_indices_.push_back(cit->vertex(2)->info());
		Power_shape_data.volumes_vertex_indices_.push_back(cit->vertex(3)->info());
	}

	import_volume_data(*mv, Power_shape_data);
	import_point_data(*mp, Voronoi_Vertices_data);
	
	// 2. For each sample point, compute its poles .

	
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
	Mesh_Point* mp = new Mesh_Point{};

	test_delaunay(mp, mv, surface_mesh);

 	mpv.register_mesh(mv, "delaunay tredrehedra");
 	mpp.register_mesh(mp, "power vertices");


	std::shared_ptr<Attribute_Volumn<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex_Volumn>(*mv, "position");
	mpv.set_mesh_bb_vertex_position(*mv, vertex_position);
	vr.set_vertex_position(*v1, *mv, vertex_position);
	
	std::shared_ptr<Attribute_Point<Vec3>> vertex_position_2 =
		cgogn::get_attribute<Vec3, Vertex_Point>(*mp, "position");
	mpp.set_mesh_bb_vertex_position(*mp, vertex_position_2);
	pcr.set_vertex_position(*v1, *mp, vertex_position_2);

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
