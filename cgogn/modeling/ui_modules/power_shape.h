/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_POWER_SHAPE_H_
#define CGOGN_MODULE_POWER_SHAPE_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn\modeling\algos\decimation\SQEM_helper.h>
#include <cgogn/io/point/point_import.h>
#include <cgogn/geometry/types/slab_quadric.h>
#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>
// import CGAL
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/double.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Bbox_3.h>
#include <Highs.h>
#include <iomanip>
namespace cgogn
{

namespace ui
{

template <typename POINT, typename SURFACE, typename NONMANIFOLD>
class PowerShape : public Module
{
	static_assert(mesh_traits<SURFACE>::dimension >= 2, "PowerShape can only be used with meshes of dimension >= 2");
	// Kernel for construct Delaunay
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point = K::Point_3;
	using Weight_Point = K::Weighted_point_3;

	class VertexInfo
	{
	public:
		Point inside_pole;
		Point outside_pole;
		double inside_pole_distance = 0.0;
		double outside_pole_distance = 0.0;
	};

	class DelaunayCellInfo
	{
	public:
		uint32 id = -1;
		bool inside = false;
		Point centroid;
		double radius2;//radius square
	};
	class RegularVertexInfo
	{
	public:
		bool inside = false;
		uint32 id = -1;
	};

	using Vb = CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>;
	// Delaunay
	using Cb = CGAL::Triangulation_cell_base_with_info_3<DelaunayCellInfo,K>;
	using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
	using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
	// Regular
	using Vb0 = CGAL::Regular_triangulation_vertex_base_3<K>;
	using RVb = CGAL::Triangulation_vertex_base_with_info_3<RegularVertexInfo, K, Vb0>;
	using RCb = CGAL::Regular_triangulation_cell_base_3<K>;
	using RTds = CGAL::Triangulation_data_structure_3<RVb, RCb>;
	using Regular = CGAL::Regular_triangulation_3<K, RTds, CGAL::Fast_location>;
	
	using Cgal_Surface_mesh = CGAL::Surface_mesh<Point>;
	using Point_inside = CGAL::Side_of_triangle_mesh<Cgal_Surface_mesh, K>;
	using Primitive = CGAL::AABB_face_graph_triangle_primitive<Cgal_Surface_mesh>;
	using Traits = CGAL::AABB_traits<K, Primitive>;
	using Tree = CGAL::AABB_tree<Traits>;

	template <typename T>
	using NonManifoldAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using PointVertex = typename mesh_traits<POINT>::Vertex;
	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NonManifoldFace = typename mesh_traits<NONMANIFOLD>::Face;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	using Vec4 = geometry::Vec4;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	PowerShape(const App& app) : Module(app, "PowerShape"), selected_surface_mesh_(nullptr)
	{
	}
	~PowerShape()
	{
	}

private:
	struct point_hash
	{
		std::size_t operator()(const Point& p) const
		{
			return ((std::hash<double>()(p.x()) ^ (std::hash<double>()(p.y()) << 1)) >> 1) ^
				   (std::hash<double>()(p.z()) << 1);
		}
	};

	struct edge_hash
	{
		std::size_t operator()(const std::pair<uint32, uint32>& edge) const
		{
			return std::hash<uint32>()(edge.first) + std::hash<uint32>()(edge.second);
		}
	};
	struct edge_equal
	{
		bool operator()(const std::pair<uint32, uint32>& edge1, const std::pair<uint32, uint32>& edge2) const
		{
			return ((edge1.first == edge2.first && edge1.second == edge2.second) ||
				(edge1.first == edge2.second && edge1.second == edge2.first));
		}
	};
	bool inside_sphere(const Vec3& point, const Vec3& center, double radius)
	{
		return (point - center).norm() <= radius;
	}
	bool pointInside(Tree& tree, Point& query)
	{
		// Initialize the point-in-polyhedron tester
		Point_inside inside_tester(tree);

		// Determine the side and return true if inside!
		return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
	}

	void load_model_in_cgal(SURFACE& surface, Cgal_Surface_mesh& csm)
	{
		std::string filename = surface_provider_->mesh_filename(surface);
		if (!filename.empty())
		{
			std::ifstream input(filename);
			if (!input || !(input >> csm))
			{
				std::cerr << "Error: input file could not be read" << std::endl;
				return;
			}
		}
	}

public:
/*
 	void detect_boundary_cells(NONMANIFOLD& nm)
 	{
 		parallel_foreach_cell(nm, [&](NonManifoldVertex v) {
 			auto ie = incident_edges(nm, v);
 			auto iface = incident_faces(nm, v);
 			set_boundary(nm, v, ie.size() == 1 && iface.size() == 0);
 			return true;
 			});
 		parallel_foreach_cell(nm, [&](NonManifoldEdge e) {
 			auto iface = incident_faces(nm, e);
 			set_boundary(nm, e, iface.size() == 1);
 
 			return true;
 			});
 	}
*/

	std::array<std::array<double, 3>, 8> compute_big_box(SURFACE &surface, Cgal_Surface_mesh& csm)
	{
		std::array<Point, 8> obb_points;
		Point acc(0, 0, 0);
		CGAL::oriented_bounding_box(csm, obb_points, CGAL::parameters::use_convex_hull(true));
		for (size_t i = 0; i < obb_points.size(); i++)
		{
			acc += K::Vector_3(obb_points[i].x(), obb_points[i].y(), obb_points[i].z());
		}
		std::array<double, 3> center{acc.x() / 8, acc.y() / 8, acc.z() / 8};
		// Create a large box surrounding object so that the Voronoi vertices are bounded
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		double offset = (md.bb_max_ - md.bb_min_).norm() * 20.0;
		std::array<double, 3> offset_array = {offset, offset, offset};
		std::array<std::array<double, 3>, 8> cube_corners = {
			{{center[0] - offset_array[0], center[1] - offset_array[1], center[2] - offset_array[2]},
			 {center[0] - offset_array[0], center[1] - offset_array[1], center[2] + offset_array[2]},
			 {center[0] - offset_array[0], center[1] + offset_array[1], center[2] - offset_array[2]},
			 {center[0] - offset_array[0], center[1] + offset_array[1], center[2] + offset_array[2]},
			 {center[0] + offset_array[0], center[1] - offset_array[1], center[2] - offset_array[2]},
			 {center[0] + offset_array[0], center[1] - offset_array[1], center[2] + offset_array[2]},
			 {center[0] + offset_array[0], center[1] + offset_array[1], center[2] - offset_array[2]},
			 {center[0] + offset_array[0], center[1] + offset_array[1], center[2] + offset_array[2]}}};
		return cube_corners;
	}

	void plot_surface_samples(std::vector<Point>& mesh_samples)
	{
		cgogn::io::PointImportData samples;
		uint32 nb_vertices = mesh_samples.size();
		samples.reserve(nb_vertices);
		for (auto& s : mesh_samples)
		{
			samples.vertex_position_.emplace_back(s[0], s[1], s[2]);
		}
		cgogn::io::import_point_data(*surface_sample_, samples);
		auto position = get_attribute<Vec3, PointVertex>(*surface_sample_, "position");
		if (position)
			point_provider_->set_mesh_bb_vertex_position(*surface_sample_, position);
	}

	Delaunay compute_delaunay_tredrahedron(SURFACE& surface, Cgal_Surface_mesh& csm, Tree& tree)
	{
	
		std::vector<Point> Delaunay_tri_point;
		
		// Sampling the mesh surface
		std::vector<Point> mesh_samples;
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(
			csm, std::back_inserter(mesh_samples),
			//CGAL::parameters::use_monte_carlo_sampling(true).number_of_points_per_area_unit(50));
			CGAL::parameters::use_grid_sampling(true).grid_spacing(1));

 		/*std::array<std::array<double, 3>, 8> cube_corners = compute_big_box(surface, csm);
 		// 	Add bounding box vertices in the sample points set
 		for (auto& p : cube_corners)
 		{
 			Delaunay_tri_point.emplace_back(p[0], p[1], p[2]);
 		}*/

		// Add sampled vertices into the volume data to construct the delauney tredrahedron
		for (auto& s : mesh_samples)
		{
			Delaunay_tri_point.emplace_back(s[0], s[1], s[2]);
		}
		
		plot_surface_samples(mesh_samples);		

		//auto start_timer = std::chrono::high_resolution_clock::now();
		
		// Construct delauney tredrahedron using CGAL
		int vertex_index = 0;
		int inside_cell_index = 0;
		int general_cell_index = 0;
		Delaunay tri;
		for (Point p : Delaunay_tri_point)
		{
			Delaunay::Vertex_handle cell= tri.insert(p);
			
		}
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			cit->info().id = -1;
			cit->info().centroid = CGAL::circumcenter(tri.tetrahedron(cit));
			cit->info().radius2 = CGAL::squared_distance(cit->info().centroid, cit->vertex(0)->point());
			
			if (pointInside(tree, cit->info().centroid))
			{
				cit->info().inside = true;
				cit->info().id = inside_cell_index;
				inside_cell_index++;
			}
			else
			{
				cit->info().inside = false;
			}
		}
		return tri;
	}

	void compute_initial_non_manifold(Delaunay& tri, Tree& tree, string name)
	{
		std::vector<double> sphere_radius;
		std::vector<Point> sphere_center;
		cgogn::io::IncidenceGraphImportData Initial_non_manifold;
		std::unordered_map<std::pair<uint32, uint32>, size_t, edge_hash, edge_equal> edge_indices;
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(name +
														  std::to_string(nonmanifold_provider_->number_of_meshes()));
		uint32 vertex_count = 0, edge_count = 0;
		// Add vertices
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			if (cit->info().inside)
			{
				Point centroid = cit->info().centroid;
				double radius = std::sqrt(cit->info().radius2);
				Initial_non_manifold.vertex_position_.emplace_back(centroid[0], centroid[1], centroid[2]);
				sphere_radius.push_back(radius);
				sphere_center.push_back(centroid);
			
			}
		}
		// Add edges
		for (auto fit = tri.finite_facets_begin(); fit != tri.finite_facets_end(); ++fit)
		{
			if (fit->first->info().inside && tri.mirror_facet(*fit).first->info().inside){
				uint32 v1_ind = fit->first->info().id;
				uint32 v2_ind = tri.mirror_facet(*fit).first->info().id;
				Initial_non_manifold.edges_vertex_indices_.push_back(v1_ind);
				Initial_non_manifold.edges_vertex_indices_.push_back(v2_ind);
				edge_indices.insert({{v1_ind, v2_ind}, edge_indices.size()});
			}
		}
		bool all_finite_inside;
		
		std::vector<Delaunay::Cell_handle> incells;
		for (auto eit = tri.finite_edges_begin(); eit != tri.finite_edges_end(); ++eit)
		{
			all_finite_inside = true;
			incells.clear();
 			Delaunay::Cell_circulator& cc = tri.incident_cells(*eit);
 			do
 			{
 				if (tri.is_infinite(cc))
 				{
 					all_finite_inside = false;
 					break;
 				}
 				else if (cc->info().inside == false)
 				{
 					all_finite_inside = false;
 					break;
 				}
				incells.push_back(cc);
				cc++;
			} while (cc != tri.incident_cells(*eit));
			if (!all_finite_inside)
				continue;
			for (size_t k = 2; k < incells.size() - 1; k++)
			{
				uint32 ev1 = incells[0]->info().id;
				uint32 ev2 = incells[k]->info().id;
				// Check if the edge is already added
				if (edge_indices.find({ev1, ev2}) == edge_indices.end())
				{
					Initial_non_manifold.edges_vertex_indices_.push_back(ev1);
					Initial_non_manifold.edges_vertex_indices_.push_back(ev2);
					edge_indices.insert({{ev1, ev2}, edge_indices.size()});
				}
			}
			for (size_t k = 1; k < incells.size() - 1; k++)
			{
				uint32 v1 = incells[0]->info().id;
				uint32 v2 = incells[k]->info().id;
				uint32 v3 = incells[k+1]->info().id;
				uint32 e1, e2, e3;
				e1 = edge_indices[{v1, v2}];
				e2 = edge_indices[{v2, v3}];
				e3 = edge_indices[{v3, v1}];
				
				if (e1 == e2 || e2 == e3 || e1 == e3)
				{
					std::cout << "e1: " << e1 << ", e2: " << e2 << ", e3: " << e3 << std::endl;
					std::cout << "v1: " << v1 << ", v2: " << v2 << ", v3: " << v3 << std::endl;
				}
				Initial_non_manifold.faces_nb_edges_.push_back(3);
				Initial_non_manifold.faces_edge_indices_.push_back(e1);
				Initial_non_manifold.faces_edge_indices_.push_back(e2);
				Initial_non_manifold.faces_edge_indices_.push_back(e3);
			}
		}
		
		uint32 Initial_non_manifold_nb_vertices = Initial_non_manifold.vertex_position_.size();
		uint32 Initial_non_manifold_nb_edges = Initial_non_manifold.edges_vertex_indices_.size() / 2;
		uint32 Initial_non_manifold_nb_faces = Initial_non_manifold.faces_nb_edges_.size();
		Initial_non_manifold.set_parameter(Initial_non_manifold_nb_vertices, Initial_non_manifold_nb_edges,
										   Initial_non_manifold_nb_faces);

		import_incidence_graph_data(*mv, Initial_non_manifold);

		/*foreach_cell(*mv, [&](NonManifoldFace f) {
			
			auto ifv = incident_vertices(*mv, f);
			std::set<NonManifoldVertex> sv(ifv.begin(), ifv.end());
			if (sv.size() != 3)
			{
				std::cout << "face: " << index_of(*mv, f) << " is not a triangle" << std::endl;
				for (NonManifoldVertex ifaceiiv : ifv)
				{
					std::cout << " vertex: " << index_of(*mv, ifaceiiv) << ", ";
				}
				std::cout << std::endl;
			}
			std::cout << "-------------------------" << std::endl;
			return true;
			});*/

		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
		auto sphere_info = add_attribute<Vec4, NonManifoldVertex>(*mv, "sphere_info");
		for (size_t idx = 0; idx < sphere_center.size(); idx++)
		{
			(*sphere_raidus)[idx] = sphere_radius[idx];
			(*sphere_info)[idx] =
				Vec4(sphere_center[idx].x(), sphere_center[idx].y(), sphere_center[idx].z(), sphere_radius[idx]);
		}
		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}

	Regular compute_regular_tredrahedron(Tree& tree, Delaunay& tri)
	{
		Regular power_shape;
		uint32 inside_pole_index = 0;
		uint32 count_p = 0;
		bool outside = false;
		std::vector<Delaunay::Cell_handle> vertex_incident_cells;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			outside = false;
			vertex_incident_cells.clear();
			tri.finite_incident_cells(vit, std::back_inserter(vertex_incident_cells));
			for (auto cell : vertex_incident_cells)
			{
				
				if (cell->info().inside)
				{
					if (cell->info().radius2 > vit->info().inside_pole_distance)
					{
						vit->info().inside_pole_distance = cell->info().radius2;
						vit->info().inside_pole = cell->info().centroid;
					}
				}
				else
				{
					
					if (cell->info().radius2 > vit->info().outside_pole_distance)
					{
						vit->info().outside_pole_distance = cell->info().radius2;
						vit->info().outside_pole = cell->info().centroid;
						outside = true;
					}
				}
			}
			Weight_Point wp1(vit->info().inside_pole, vit->info().inside_pole_distance);
		    power_shape.insert(wp1);
			Weight_Point wp2(vit->info().outside_pole, vit->info().outside_pole_distance);
			power_shape.insert(wp2);
			
		}
		int count = 0;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			Point p = vit->point().point();
			if (pointInside(tree, p))
			{
				vit->info().id = count;
				vit->info().inside = true;
				count++;
			}

		}
		return power_shape;
	}

	void constrcut_power_shape_non_manifold(Regular& power_shape,string name)
	{
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(name +
														  std::to_string(nonmanifold_provider_->number_of_meshes()));
		cgogn::io::IncidenceGraphImportData Inner_Power_shape_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		uint32 edge_count = 0;
		std::vector<Weight_Point> power_point;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			if (vit->info().inside)
			{
				power_point.push_back(vit->point());
				Inner_Power_shape_data.vertex_position_.push_back(
					{vit->point().x(), vit->point().y(), vit->point().z()});
			}
		}
	
		for (auto eit = power_shape.finite_edges_begin(); eit != power_shape.finite_edges_end(); ++eit)
		{
			Regular::Vertex_handle v1 = eit->first->vertex(eit->second);
			Regular::Vertex_handle v2 = eit->first->vertex(eit->third);
			if (v1->info().inside && v2->info().inside)
			{
				// Add edge
				uint32 v1_ind = v1->info().id;
				uint32 v2_ind = v2->info().id;
				edge_indices.insert({{v1_ind, v2_ind}, edge_count});
				Inner_Power_shape_data.edges_vertex_indices_.push_back(v1_ind);
				Inner_Power_shape_data.edges_vertex_indices_.push_back(v2_ind);
				edge_count++;
			}
		}

		for (auto fit = power_shape.finite_facets_begin(); fit != power_shape.finite_facets_end(); ++fit)
		{
			Regular::Vertex_handle v[3];
			int count = 0;
			bool inside = true;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != fit->second)
				{
					v[count]= fit->first->vertex(idx);
					inside &= v[count]->info().inside;
					count++;
				}
			}
			if (inside)
			{
				Inner_Power_shape_data.faces_nb_edges_.push_back(3);
				Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v[0]->info().id, v[1]->info().id}]);
				Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v[1]->info().id, v[2]->info().id}]);
				Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v[2]->info().id, v[0]->info().id}]);
				
			}
		}
		uint32 inner_power_nb_vertices = Inner_Power_shape_data.vertex_position_.size();
		uint32 inner_power_nb_edges = Inner_Power_shape_data.edges_vertex_indices_.size() / 2;
		uint32 inner_power_nb_faces = Inner_Power_shape_data.faces_nb_edges_.size();

		Inner_Power_shape_data.set_parameter(inner_power_nb_vertices, inner_power_nb_edges, inner_power_nb_faces);

		import_incidence_graph_data(*mv, Inner_Power_shape_data);
		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
		auto sphere_info = add_attribute<Vec4, NonManifoldVertex>(*mv, "sphere_info");
		for (uint32 idx = 0u; idx < power_point.size(); ++idx)
		{
			(*sphere_raidus)[idx] = std::sqrt(power_point[idx].weight());
			(*sphere_info)[idx] =
				Vec4(power_point[idx].x(), power_point[idx].y(), power_point[idx].z(), power_point[idx].weight());
			
		}

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}

	void compute_power_shape(SURFACE& surface)
	{
		surface_sample_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "surface_samples");
		
		Cgal_Surface_mesh csm;
		load_model_in_cgal(surface, csm);
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Delaunay tri = compute_delaunay_tredrahedron(surface, csm, tree);
		Regular reg = compute_regular_tredrahedron(tree, tri);
		constrcut_power_shape_non_manifold( reg,surface_provider_->mesh_name(surface) +"_inner_power_shape");
	}

	void compute_original_power_diagram(SURFACE& surface)
	{
		surface_sample_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "surface_samples");
		Cgal_Surface_mesh csm;
		load_model_in_cgal(surface, csm);
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Delaunay tri = compute_delaunay_tredrahedron(surface, csm, tree);
		compute_initial_non_manifold(tri, tree, surface_provider_->mesh_name(surface)+ "_inner_voronoi_diagram");
	}



	void compute_stability_ratio_edge(NONMANIFOLD& nm, NonManifoldEdge e,
									  std::shared_ptr<NonManifoldAttribute<double>>& stability_ratio,
									  std::shared_ptr<NonManifoldAttribute<Vec3>>& stability_color,
									  std::shared_ptr<NonManifoldAttribute<Vec4>>& sphere_info)
	{
		auto iv = incident_vertices(nm, e);
		NonManifoldVertex v1 = iv[0];
		NonManifoldVertex v2 = iv[1];
		const Vec3& v1_p = value<Vec4>(nm, sphere_info, v1).head<3>();
		const Vec3& v2_p = value<Vec4>(nm, sphere_info, v2).head<3>();
		const double& r1 = value<Vec4>(nm, sphere_info, v1).w();
		const double& r2 = value<Vec4>(nm, sphere_info, v2).w();
		const double center_dist = (v1_p - v2_p).norm();
		double dis = std::max(0.0, (center_dist - std::abs(r1 - r2)));
		if (center_dist == 0.0)
		{
			(*stability_ratio)[e.index_] = 0.0;
			(*stability_color)[e.index_] = Vec3(0, 0, 0.5);
			return;
		}
		double stability = dis / center_dist;
		value<double>(nm, stability_ratio,e) = stability;
		value<Vec3>(nm, stability_color, e) =
			(stability <= 0.5) ? Vec3(0, stability, (0.5 - stability)) : Vec3(stability - 0.5, (1 - stability), 0);
	}

	void compute_stability_ratio(NONMANIFOLD& nm)
	{
		auto stability_ratio = add_attribute<double, NonManifoldEdge>(nm, "stability_ratio");
		auto stability_color = add_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		auto sphere_info =
			get_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info"); // {center, radius} = {x, y, z, r}
		parallel_foreach_cell(nm, [&](NonManifoldEdge e) -> bool {
			compute_stability_ratio_edge(nm, e, stability_ratio, stability_color, sphere_info);
			return true;
		});
	}

	 void collapse_non_manifold_using_QMat(NONMANIFOLD& nm, uint32 number_vertices_remain, float k)
	{
		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;
		using Slab_Quadric = geometry::Slab_Quadric;
		
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(nm, "sphere_radius");
		auto stability_ratio = get_attribute<double, NonManifoldEdge>(nm, "stability_ratio");
		auto stability_color = get_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		auto position = get_attribute<Vec3, NonManifoldVertex>(nm, "position");
		auto sphere_info = get_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info");
		auto fixed_vertex = add_attribute<bool, NonManifoldVertex>(nm, "fixed_vertex");
		foreach_cell(nm, [&](NonManifoldVertex v) -> bool {
			value<bool>(nm, fixed_vertex, v) = false;
			return true;
		});

		QMatHelper helper(k, nm, position, sphere_info, stability_color, stability_ratio, sphere_radius,fixed_vertex); 

		helper.initial_slab_mesh();
		helper.initial_boundary_mesh();
		helper.initial_collapse_queue();
		helper.simplify(number_vertices_remain, true);
		
		remove_attribute<NonManifoldVertex>(nm, fixed_vertex);
		nonmanifold_provider_->emit_connectivity_changed(nm);
		nonmanifold_provider_->emit_attribute_changed(nm, position.get());
		nonmanifold_provider_->emit_attribute_changed(nm, sphere_radius.get());
		
	}


	void coverage_axis_PD(SURFACE& surface, NONMANIFOLD& mv, HighsSolution& solution, double dilation_factor)
	{
		Cgal_Surface_mesh csm;
		load_model_in_cgal(surface, csm);
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		auto inner_position = get_attribute<Vec3, NonManifoldVertex>(mv, "position");
		auto sample_position = get_attribute<Vec3, Vertex>(surface, "position");
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(mv, "sphere_radius");
		std::vector<Weight_Point> power_point;
		Regular power_shape;
		foreach_cell(mv, [&](NonManifoldVertex nv) {
			if (solution.col_value[index_of(mv, nv)] > 1e-5)
			{
				Vec3 pos = value<Vec3>(mv, inner_position, nv);
				power_shape.insert(Weight_Point(Point(pos[0], pos[1], pos[2]),
												(value<double>(mv, sphere_radius, nv) + dilation_factor) *
													(value<double>(mv, sphere_radius, nv) + dilation_factor)));
			
				std::cout << "Vertex: " << index_of(mv, nv) << " position is : " << pos.x() << " " << pos.y() << " "
						  << pos.z()
						  << std::endl;
			}
			return true;
		});

		foreach_cell(surface, [&](Vertex v) {
			Vec3 pos = value<Vec3>(surface, sample_position, v);
		
			power_shape.insert(
				Weight_Point(Point(pos[0], pos[1], pos[2]), (dilation_factor + 0.001) * (dilation_factor + 0.001)));
			
			return true;
		});
		int count = 0;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			Point p = vit->point().point();
			if (pointInside(tree, p))
			{
				vit->info().id = count;
				vit->info().inside = true;
				count++;
				
			}
		}
		constrcut_power_shape_non_manifold(power_shape,surface_provider_->mesh_name(surface)+ "_coverage_axis");
	}

	void coverage_axis_collapse(NONMANIFOLD& nm, HighsSolution& solution)
	{
		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;
		using Slab_Quadric = geometry::Slab_Quadric;

		compute_stability_ratio(nm);
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(nm, "sphere_radius");
		auto stability_ratio = get_attribute<double, NonManifoldEdge>(nm, "stability_ratio");
		auto stability_color = get_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		auto position = get_attribute<Vec3, NonManifoldVertex>(nm, "position");
		auto sphere_info = get_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info");
		auto fixed_vertex = add_attribute<bool, NonManifoldVertex>(nm, "fixed_vertex");
		int number_vertex_remain = 0;
		foreach_cell(nm, [&](NonManifoldVertex nv) -> bool {
			Vec3 pos = value<Vec3>(nm, position, nv);
			if (solution.col_value[index_of(nm, nv)] > 1e-5)
			{
				value<bool>(nm, fixed_vertex, nv) = true;
				number_vertex_remain++;
			}
			else
				value<bool>(nm, fixed_vertex, nv) = false;
			
			return true;
		});
		QMatHelper helper(0.00001, nm, position, sphere_info, stability_color, stability_ratio, sphere_radius, fixed_vertex);

		helper.initial_slab_mesh();
		helper.initial_boundary_mesh();
		helper.initial_collapse_queue();
		helper.simplify(number_vertex_remain, false);

		foreach_cell(nm, [&](NonManifoldVertex nv) {
			Vec3 pos = value<Vec3>(nm, position, nv);
			return true;
		});
		
		remove_attribute<NonManifoldVertex>(nm, fixed_vertex);
		nonmanifold_provider_->emit_connectivity_changed(nm);
		nonmanifold_provider_->emit_attribute_changed(nm, position.get());
		nonmanifold_provider_->emit_attribute_changed(nm, sphere_radius.get());
	}
	
 	 HighsSolution point_selection_by_coverage_axis(SURFACE& surface, NONMANIFOLD& mv, double dilation_factor)
	{
		typedef Eigen::SparseMatrix<double> SpMat; 
		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;
		auto inner_position = get_attribute<Vec3, NonManifoldVertex>(mv, "position");
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(mv, "sphere_radius");
		auto sample_position = get_attribute<Vec3, Vertex>(surface, "position");
		auto inner_point_nb = nb_cells<NonManifoldVertex>(mv);
		auto sample_point_nb = nb_cells<Vertex>(surface);
		foreach_cell(mv, [&](NonManifoldVertex nv) {
			foreach_cell(surface, [&](Vertex v) { 
				if (inside_sphere(value<Vec3>(surface, sample_position, v), value<Vec3>(mv, inner_position, nv),
								  value<double>(mv, sphere_radius, nv) + dilation_factor))
				{
					triplets.push_back(T(index_of(surface, v), index_of(mv, nv), 1.0));
				}
				return true; });
			return true;
		});
		SpMat A(sample_point_nb, inner_point_nb);
		A.setFromTriplets(triplets.begin(), triplets.end());
		A.makeCompressed();
		HighsModel model;
		model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
		model.lp_.num_col_ = A.cols();
		model.lp_.num_row_ = A.rows();
		model.lp_.sense_ = ObjSense::kMinimize;
		// Adding decision variable bounds
		HighsVarType type = HighsVarType::kInteger;
		model.lp_.col_cost_ = std::vector<double>(model.lp_.num_col_, 1.0);
		model.lp_.col_lower_ = std::vector<double>(model.lp_.num_col_, 0.0);
		model.lp_.col_upper_ = std::vector<double>(model.lp_.num_col_, 1.0);
		model.lp_.row_lower_ = std::vector<double>(model.lp_.num_row_, 1.0);
		model.lp_.row_upper_ = std::vector<double>(model.lp_.num_row_, 1e30);
		model.lp_.integrality_ = std::vector<HighsVarType>(model.lp_.num_col_, type);
		
		model.lp_.a_matrix_.num_col_ = model.lp_.num_col_;
		model.lp_.a_matrix_.num_row_ = model.lp_.num_row_;
		model.lp_.a_matrix_.start_.resize(A.cols() + 1);
		model.lp_.a_matrix_.index_.resize(A.nonZeros());
		model.lp_.a_matrix_.value_.resize(A.nonZeros());
		
		// Copy the data from A to the vectors
		std::copy(A.outerIndexPtr(), A.outerIndexPtr() + A.cols() + 1, model.lp_.a_matrix_.start_.begin());
		std::copy(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros(), model.lp_.a_matrix_.index_.begin());
		std::copy(A.valuePtr(), A.valuePtr() + A.nonZeros(), model.lp_.a_matrix_.value_.begin());
		Highs highs;
		HighsStatus status = highs.passModel(model);
		HighsSolution solution; 
		highs.setOptionValue("time_limit", 60);
		if (status == HighsStatus::kOk)
		{
			highs.run();
			
			assert(status == HighsStatus::kOk);
			solution = highs.getSolution();

 			/*Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
			
 			for (int i = 0; i < solution.col_value.size(); ++i)
 			{
 				x[i] = solution.col_value[i];
 			}

			std::cout << " solution.row_value.size(): " << solution.row_value.size() << std::endl;
			for (int i = 0; i < solution.row_value.size(); ++i)
			{
				std::cout << solution.row_value[i] << " ";
			}
			std::cout << std::endl;
			std::cout << "compare " << std::endl;
			std::cout << A * x << std::endl;
			 Get the primal solution values
			
 			auto ca_sphere_radius = get_attribute<double, NonManifoldVertex>(*ca, "sphere_radius");
 			foreach_cell(*ca, [&](NonManifoldVertex v) { 
 				value<double>(*ca, sphere_radius, v) += dilation_factor;
 				return true;
 			});

 			auto ca_sphere_center = get_attribute<Vec3, NonManifoldVertex>(*ca, "position");
 			uint32 outside_count = 0;
 			foreach_cell(surface, [&](Vertex v) { 
 				Vec3 pos = value<Vec3>(surface, sample_position, v);
 				bool inside = false;
 				foreach_cell(*ca, [&](NonManifoldVertex nv) {
 					Vec3 capos = value<Vec3>(*ca, ca_sphere_center, nv);
 					double radius = value<double>(*ca, ca_sphere_radius, nv);
 					inside |= inside_sphere(pos, capos, radius);
 					return true;
 				});
 				if (!inside)
 				{
 					outside_count++;
 					std::cout << "Fault! " << outside_count << std::endl;
 				}
 					return true;
 			});

			construct_complete_power_diagram(ca_complet, power_point, point_info);*/
		}
		return solution;
	}
	
protected:
	void init() override
	{
		point_provider_ = static_cast<ui::MeshProvider<POINT>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINT>::name} + ")"));
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
		nonmanifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_mesh_, "Surface", [&](SURFACE& m) {
			selected_surface_mesh_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
		imgui_mesh_selector(nonmanifold_provider_, selected_medial_axis_, "Medial_axis", [&](NONMANIFOLD& nm) {
			selected_medial_axis_ = &nm;
			nonmanifold_provider_->mesh_data(nm).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_surface_mesh_)
		{
			if (ImGui::Button("Power shape"))
			{
				compute_power_shape(*selected_surface_mesh_);
				
			}
			if (ImGui::Button("Original Medial Axis"))
				compute_original_power_diagram(*selected_surface_mesh_);
			if (selected_medial_axis_)
			{
				if (ImGui::Button("Compute stability ratio"))
					compute_stability_ratio(*selected_medial_axis_);
				static int32 number_vertex_remain = 1; 
				static float k = 1e-5;
				ImGui::DragInt("Vertices to delete", &number_vertex_remain, 1, 0,
								  nb_cells<NonManifoldVertex>(*selected_medial_axis_));
				ImGui::DragFloat("K", &k, 1e-5, 0.0f, 1.0f, "%.5f"); 
				 if (ImGui::Button("QMAT"))
				{
					collapse_non_manifold_using_QMat(*selected_medial_axis_, number_vertex_remain, k);
					
				}
				static float dilation_factor = 0.02f;
				ImGui::DragFloat("Dilation factor", &dilation_factor, 0.001f, 0.0f, 1.0f, "%.4f");
				if (ImGui::Button("Coverage Axis"))
				{
					solution = point_selection_by_coverage_axis(*selected_surface_mesh_, *selected_medial_axis_, dilation_factor);
					
				}
				if (solution.col_value.size() > 0)
				{

					if (ImGui::Button("Collpase"))
						coverage_axis_collapse(*selected_medial_axis_, solution);
					if (ImGui::Button("PD"))
						coverage_axis_PD(*selected_surface_mesh_, *selected_medial_axis_, solution, dilation_factor);
				}
			}
				
				
			
		}
	}

private:
	POINT* surface_sample_ = nullptr;
	SURFACE* selected_surface_mesh_ = nullptr;
	NONMANIFOLD* selected_medial_axis_ = nullptr;
	std::shared_ptr<NonManifoldAttribute<double>> stability_ratio_ = nullptr;
	MeshProvider<POINT>* point_provider_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
	Regular medial_axis;
	HighsSolution solution;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POWER_SHAPE_H_
