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
// import CGAL
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/optimal_bounding_box.h>

#include <Highs.h>

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
	using Vb = CGAL::Triangulation_vertex_base_with_info_3<uint32_t, K>;
	// Delaunay
	using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
	using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
	using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
	// Regular
	using Vb0 = CGAL::Regular_triangulation_vertex_base_3<K>;
	using RVb = CGAL::Triangulation_vertex_base_with_info_3<std::pair<uint32_t, bool>, K, Vb0>;
	using RCb = CGAL::Regular_triangulation_cell_base_3<K>;
	using RTds = CGAL::Triangulation_data_structure_3<RVb, RCb>;
	using Regular = CGAL::Regular_triangulation_3<K, RTds, CGAL::Fast_location>;

	using Point = K::Point_3;
	using Weight = K::FT;
	using Weight_Point = K::Weighted_point_3;

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
	PowerShape(const App& app)
		: Module(app, "PowerShape"), selected_surface_mesh_(nullptr)
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
	struct face_hash
	{
		std::size_t operator()(const std::tuple<uint32, uint32, uint32>& triangle) const
		{
			std::size_t seed = 0;
			seed ^= std::hash<uint32>{}(std::get<0>(triangle)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= std::hash<uint32>{}(std::get<1>(triangle)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= std::hash<uint32>{}(std::get<2>(triangle)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			return seed;
		}
	};

	struct face_equal
	{
		bool operator()(const std::tuple<uint32, uint32, uint32>& lhs,
						const std::tuple<uint32, uint32, uint32>& rhs) const
		{
			std::array<uint32, 3> lhsIndices = {std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs)};
			std::array<uint32, 3> rhsIndices = {std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs)};

			std::sort(lhsIndices.begin(), lhsIndices.end());
			std::sort(rhsIndices.begin(), rhsIndices.end());

			return lhsIndices == rhsIndices;
		}
	};

	bool pointInside(Tree& tree, Point& query)
	{
		// Initialize the point-in-polyhedron tester
		Point_inside inside_tester(tree);

		// Determine the side and return true if inside!
		return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
	}

public:
	Delaunay compute_delaunay_tredrahedron(SURFACE& surface, Cgal_Surface_mesh& csm, Tree& tree)
	{
		cgogn::io::PointImportData samples;
		std::vector<Point> Delaunay_tri_point;
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
		double offset = (md.bb_max_ - md.bb_min_).norm() *0.7;
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

		// Sampling the mesh surface
		std::vector<Point> mesh_samples;
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(
			csm, std::back_inserter(mesh_samples),
		    //CGAL::parameters::use_monte_carlo_sampling(true).number_of_points_per_area_unit(50));
			CGAL::parameters::use_grid_sampling(true).grid_spacing(1));

		// 	Add bounding box vertices in the sample points set
		for (auto& p : cube_corners)
		{
			Delaunay_tri_point.emplace_back(p[0], p[1], p[2]);
			samples.vertex_position_.emplace_back(p[0], p[1], p[2]);
		}

		// Add sampled vertices into the volume data to construct the delauney tredrahedron
		for (auto& s : mesh_samples)
		{
			Delaunay_tri_point.emplace_back(s[0], s[1], s[2]);
			samples.vertex_position_.emplace_back(s[0], s[1], s[2]);
		}

		uint32 nb_vertices = Delaunay_tri_point.size();

		//		auto start_timer = std::chrono::high_resolution_clock::now();

		// Indices info for constructing volume data in Cgogn
		std::vector<unsigned> indices;
		indices.reserve(nb_vertices);
		for (unsigned int i = 0; i < nb_vertices; ++i)
			indices.push_back(i);
		//Collect point data
		samples.reserve(nb_vertices);
		cgogn::io::import_point_data(*surface_sample, samples);
		auto position = get_attribute<Vec3, PointVertex>(*surface_sample, "position");
		if (position)
			point_provider_->set_mesh_bb_vertex_position(*surface_sample, position);
		// Construct delauney tredrahedron using CGAL
		Delaunay tri(boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.begin(), indices.begin())),
					 boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.end(), indices.end())));
		return tri;
	}
	
	void compute_inner_voronoi(Delaunay& tri, Tree& tree, std::vector<Weight_Point>& power_point,
							   std::vector<std::pair<uint32, bool>>& point_info,
							   std::unordered_map<uint32, uint32>& inside_indices)
	{
		uint32 count = 0, inside_vertices_count = 0;
		double dis;
		std::unordered_set<Point, point_hash> voronoi_vertices;
		std::vector<Delaunay::Cell_handle> cells;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			cells.clear();
			tri.finite_incident_cells(vit, std::back_inserter(cells));
			if (cells.size())
			{
				for (auto c = cells.begin(); c != cells.end(); ++c)
				{
					Point centroid = tri.dual(*c);
					dis = CGAL::squared_distance(centroid, vit->point());
					if (voronoi_vertices.find(centroid) == voronoi_vertices.end())
					{
						voronoi_vertices.insert(centroid);
						power_point.push_back(Weight_Point(centroid, dis));
						if (pointInside(tree, centroid))
						{
							point_info.push_back({count, true});
							inside_indices.insert({count++, inside_vertices_count++});
						}
						else if (!pointInside(tree, centroid))
						{
							point_info.push_back({count++, false});
						
						}
					}
				}
			}
		}
	}

	void compute_inner_poles(Delaunay& tri, Tree& tree, std::vector<Weight_Point>& power_point,
							 std::vector<std::pair<uint32, bool>>& point_info,
							 std::unordered_map<uint32, uint32>& inside_indices)
	{
		std::vector<std::pair<uint32, bool>> power_indices;
		// Find inside and outside poles
		uint32 count = 0, inside_poles_count = 0;
		double dis;
		Point farthest_inside_point, farthest_outside_point;
		double farthest_inside_distance, farthest_outside_distance;
		std::unordered_set<Point, point_hash> poles;
		 // Use for the construction of medial axis
		std::vector<Delaunay::Cell_handle> cells;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			farthest_inside_distance = 0;
			farthest_outside_distance = 0;
			cells.clear();
			tri.finite_incident_cells(vit, std::back_inserter(cells));
			if (cells.size())
			{
				for (auto c = cells.begin(); c != cells.end(); ++c)
				{
					Point centroid = tri.dual(*c);
					dis = CGAL::squared_distance(centroid, vit->point());
					if (pointInside(tree, centroid) && dis > farthest_inside_distance)
					{
						farthest_inside_point = centroid;
						farthest_inside_distance = dis;
					}
					else if (!pointInside(tree, centroid) && dis > farthest_outside_distance)
					{
						farthest_outside_point = centroid;
						farthest_outside_distance = dis;
					}
				}
			}
			if (farthest_inside_distance != 0)
			{
				if (poles.find(farthest_inside_point) == poles.end())
				{
					poles.insert(farthest_inside_point);
					power_point.push_back(Weight_Point(farthest_inside_point, farthest_inside_distance));
					point_info.push_back({count, true});
					inside_indices.insert({count++, inside_poles_count++});
				}
			}
			if (farthest_outside_distance != 0)
			{
				if (poles.find(farthest_outside_point) == poles.end())
				{
					poles.insert(farthest_outside_point);
					power_point.push_back(Weight_Point(farthest_outside_point, farthest_outside_distance));
					point_info.push_back({count++, false});
				}
			}
		}
	}
	void construct_complete_power_diagram(NONMANIFOLD* mv, std::vector<Weight_Point>& power_point,
										  std::vector<std::pair<uint32, bool>>& point_info)
	{
		cgogn::io::IncidenceGraphImportData Complete_Power_shape_data;

		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		uint32 edge_count = 0;

		medial_axis = Regular(boost::make_zip_iterator(boost::make_tuple(power_point.begin(), point_info.begin())),
							  boost::make_zip_iterator(boost::make_tuple(power_point.end(), point_info.end())));

		for (size_t idx = 0; idx < power_point.size(); ++idx)
		{
			Complete_Power_shape_data.vertex_position_.push_back(
				{power_point[idx].x(), power_point[idx].y(), power_point[idx].z()});
		}
		uint32 v, v1, v2;
		for (auto fit = medial_axis.finite_facets_begin(); fit != medial_axis.finite_facets_end(); ++fit)
		{
			v = fit->second;
			Complete_Power_shape_data.faces_nb_edges_.push_back(3);
			for (size_t i = 0; i < 4; i++)
			{
				if (i != v)
				{
					for (size_t j = i + 1; j < 4; j++)
					{
						if (j != v)
						{
							// Add edge
							v1 = fit->first->vertex(i)->info().first;
							v2 = fit->first->vertex(j)->info().first;
							if (edge_indices.find({v1, v2}) == edge_indices.end())
							{
								edge_indices.insert({{v1, v2}, edge_count});
								Complete_Power_shape_data.edges_vertex_indices_.push_back(v1);
								Complete_Power_shape_data.edges_vertex_indices_.push_back(v2);
								edge_count++;
							}
							// Add face
							Complete_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v1, v2}]);
						}
					}
				}
			}
			
		}
		uint32 complete_power_nb_vertices = Complete_Power_shape_data.vertex_position_.size();
		uint32 complete_power_nb_edges = Complete_Power_shape_data.edges_vertex_indices_.size() / 2;
		uint32 complete_power_nb_faces = Complete_Power_shape_data.faces_nb_edges_.size()/3;

		Complete_Power_shape_data.set_parameter(complete_power_nb_vertices, complete_power_nb_edges,
												complete_power_nb_faces);

		import_incidence_graph_data(*mv, Complete_Power_shape_data);
		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
		for (uint32 i = 0u; i < power_point.size(); ++i)
		{
			uint32 vertex_id = point_info[i].first;
			(*sphere_raidus)[vertex_id] = std::sqrt(power_point[i].weight());
		}

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}
	void constrcut_inner_power_diagram(NONMANIFOLD* mv, std::vector<Weight_Point>& power_point,
									   std::vector<std::pair<uint32, bool>>& point_info,
									   std::unordered_map<uint32, uint32>& inside_indices)
	{

		cgogn::io::IncidenceGraphImportData Inner_Power_shape_data;

		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		std::unordered_set<std::tuple<uint32, uint32, uint32>, face_hash, face_equal> face_sets;
		uint32 edge_count = 0;

		medial_axis = Regular(boost::make_zip_iterator(boost::make_tuple(power_point.begin(), point_info.begin())),
							  boost::make_zip_iterator(boost::make_tuple(power_point.end(), point_info.end())));

		for (size_t idx = 0; idx < power_point.size(); ++idx)
		{
			// if the point is inside
			if (point_info[idx].second)
			{
				Inner_Power_shape_data.vertex_position_.push_back(
					{power_point[idx].x(), power_point[idx].y(), power_point[idx].z()});
			}
		}
		bool inside;
		uint32 v, v_ind1, v_ind2, v1, v2;
		for (auto eit = medial_axis.finite_edges_begin(); eit != medial_axis.finite_edges_end(); ++eit)
		{
			v_ind1 = eit->second;
			v_ind2 = eit->third;
			inside = eit->first->vertex(v_ind1)->info().second && eit->first->vertex(v_ind2)->info().second;
			if (inside)
			{
				//Add edge
				v1 = inside_indices[eit->first->vertex(v_ind1)->info().first];
				v2 = inside_indices[eit->first->vertex(v_ind2)->info().first];
				if (edge_indices.find({v1, v2}) == edge_indices.end())
				{
					edge_indices.insert({{v1, v2}, edge_count});
					Inner_Power_shape_data.edges_vertex_indices_.push_back(v1);
					Inner_Power_shape_data.edges_vertex_indices_.push_back(v2);
					edge_count++;
				}
			}
		}
		std::tuple<uint32, uint32, uint32> face_tuple;
		std::vector<uint32> face_indices; 
		for (auto fit = medial_axis.finite_facets_begin(); fit != medial_axis.finite_facets_end(); ++fit)
		{
			face_indices.clear();
			inside = true;
			v = fit->second;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != v)
				{
					inside &= fit->first->vertex(idx)->info().second;
					face_indices.push_back(inside_indices[fit->first->vertex(idx)->info().first]);
				}
			}
			if (inside)
			{
				face_tuple = std::make_tuple(face_indices[0], face_indices[1], face_indices[2]);
				if (face_sets.find(face_tuple) == face_sets.end())
				{
					Inner_Power_shape_data.faces_nb_edges_.push_back(3);
					for (size_t i = 0; i < 4; i++)
					{
						if (i != v)
						{
							for (size_t j = i + 1; j < 4; j++)
							{
								if (j != v)
								{
									v1 = inside_indices[fit->first->vertex(i)->info().first];
									v2 = inside_indices[fit->first->vertex(j)->info().first];
									Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v1, v2}]);
								}
							}
						}
					}
					face_sets.insert(face_tuple);
				}
			}
		}
		uint32 inner_power_nb_vertices = Inner_Power_shape_data.vertex_position_.size();
		uint32 inner_power_nb_edges = Inner_Power_shape_data.edges_vertex_indices_.size() / 2;
		uint32 inner_power_nb_faces = Inner_Power_shape_data.faces_nb_edges_.size();

		Inner_Power_shape_data.set_parameter(inner_power_nb_vertices, inner_power_nb_edges, inner_power_nb_faces);

		import_incidence_graph_data(*mv, Inner_Power_shape_data);
		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
		for (uint32 i = 0u; i < power_point.size(); ++i)
		{
			uint32 vertex_id = inside_indices[point_info[i].first];
			// The radius is sqrt of the weight!
			(*sphere_raidus)[vertex_id] = std::sqrt(power_point[i].weight());
		}

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}
	//void constrcut_inner_power_diagram(NONMANIFOLD* mv, std::vector<Weight_Point>& power_point,
	//							 std::vector<std::pair<uint32, bool>>& point_info,
	//							 std::unordered_map<uint32, uint32>& inside_indices)
	//{
	//	
	//	cgogn::io::IncidenceGraphImportData Inner_Power_shape_data;
	//	
	//	std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
	//	uint32 edge_count = 0;

	//	medial_axis = Regular(boost::make_zip_iterator(boost::make_tuple(power_point.begin(), point_info.begin())),
	//				boost::make_zip_iterator(boost::make_tuple(power_point.end(), point_info.end())));

	//	for (size_t idx = 0; idx < power_point.size(); ++idx)
	//	{
	//		//if the point is inside
	//		if (point_info[idx].second)
	//		{
	//			Inner_Power_shape_data.vertex_position_.push_back(
	//				{power_point[idx].x(), power_point[idx].y(), power_point[idx].z()});
	//		}	
	//	}
	//	bool inside;
	//	uint32 v, v1, v2;
	//	for (auto fit = medial_axis.finite_facets_begin(); fit != medial_axis.finite_facets_end(); ++fit)
	//	{
	//		inside = true;
	//		v = fit->second;
	//		// If face is inside
	//		for (size_t idx = 0; idx < 4; ++idx)
	//		{
	//			if (idx != v)
	//			{
	//				inside &= fit->first->vertex(idx)->info().second;
	//			}
	//		}
	//		if (inside)
	//		{
	//			Inner_Power_shape_data.faces_nb_edges_.push_back(3);
	//			for (size_t i = 0; i < 4; i++)
	//			{
	//				if (i != v)
	//				{
	//					for (size_t j = i + 1; j < 4; j++)
	//					{
	//						if (j != v)
	//						{
	//							// Add edge
	//							v1 = inside_indices[fit->first->vertex(i)->info().first];
	//							v2 = inside_indices[fit->first->vertex(j)->info().first];
	//							if (edge_indices.find({v1, v2}) == edge_indices.end())
	//							{
	//								edge_indices.insert({{v1, v2}, edge_count});
	//								Inner_Power_shape_data.edges_vertex_indices_.push_back(v1);
	//								Inner_Power_shape_data.edges_vertex_indices_.push_back(v2);
	//								edge_count++;
	//							}
	//							// Add face
	//							Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v1, v2}]);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//	uint32 inner_power_nb_vertices = Inner_Power_shape_data.vertex_position_.size();
	//	uint32 inner_power_nb_edges = Inner_Power_shape_data.edges_vertex_indices_.size() / 2;
	//	uint32 inner_power_nb_faces = Inner_Power_shape_data.faces_nb_edges_.size();

	//	Inner_Power_shape_data.set_parameter(inner_power_nb_vertices, inner_power_nb_edges, inner_power_nb_faces);

	//	import_incidence_graph_data(*mv, Inner_Power_shape_data);
	//	auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
	//	for (uint32 i = 0u; i < power_point.size(); ++i)
	//	{
	//		uint32 vertex_id = inside_indices[point_info[i].first];
	//		//The radius is sqrt of the weight!
	//		(*sphere_raidus)[vertex_id] = std::sqrt(power_point[i].weight());
	//	}

	//	std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
	//		get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
	//	if (mv_vertex_position)
	//		nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

	//	nonmanifold_provider_->emit_connectivity_changed(*mv);
	//}

	void compute_power_shape(SURFACE& surface)
	{
		surface_sample = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "surface_samples");
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_inner_power_shape");
		NONMANIFOLD* mp = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_complete_power_shape");
		Cgal_Surface_mesh csm;

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
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Delaunay tri = compute_delaunay_tredrahedron(surface, csm, tree);
		std::vector<Weight_Point> Power_point;
		std::vector<std::pair<uint32, bool>> Point_info;
		std::unordered_map<uint32, uint32> Inside_indices;
		compute_inner_poles(tri, tree, Power_point, Point_info, Inside_indices);
		constrcut_inner_power_diagram(mv, Power_point, Point_info, Inside_indices);
		construct_complete_power_diagram(mp, Power_point, Point_info);
	}

	void compute_original_power_diagram(SURFACE& surface)
	{
		surface_sample = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "surface_samples");
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_inner_medial_axis");
		NONMANIFOLD* mp = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_complete_medial_axis");
		Cgal_Surface_mesh csm;

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
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Delaunay tri = compute_delaunay_tredrahedron(surface,csm, tree);
		std::vector<Weight_Point> Power_point;
		std::vector<std::pair<uint32, bool>> Point_info;
		std::unordered_map<uint32, uint32> Inside_indices;
		compute_inner_voronoi(tri, tree, Power_point, Point_info, Inside_indices);
 		constrcut_inner_power_diagram(mv, Power_point, Point_info, Inside_indices);
		construct_complete_power_diagram(mp, Power_point, Point_info);
	}

	bool inside_sphere(const Vec3& point, const Vec3& center, double radius)
	{
		return (point - center).norm() <= radius;
	}

 	void compute_stability_ratio(NONMANIFOLD& mv)
	{
		//TO do : better way to get attributes?
		IncidenceGraph& ig = static_cast<IncidenceGraph&>(mv);
		stability_ratio_ = add_attribute<double, NonManifoldEdge>(ig, "stability_ratio");
		auto stability_color = add_attribute<Vec3, NonManifoldEdge>(ig, "stability_color");
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(ig, "sphere_radius");
		auto position = get_attribute<Vec3, NonManifoldVertex>(ig, "position");
		
		parallel_foreach_cell(mv, [&](NonManifoldEdge e) -> bool {
			// TO do : better way to get attributes?
			auto [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];
			const Vec3& v1_p = value<Vec3>(ig, position, v1);
			const Vec3& v2_p = value<Vec3>(ig, position, v2);
			const double& r1 = value<double>(ig, sphere_radius, v1);
			const double& r2 = value<double>(ig, sphere_radius, v2);
			const double center_dist = (v1_p - v2_p).norm();
			double dis = std::max(0.0, (center_dist - std::abs(r1 - r2)));
			double stability = dis / center_dist;
			//std::cout << "Edge: " << e.index_ << ", Stability ratio: " << stability << std::flush;
			(*stability_ratio_)[e.index_] = stability;
			if (stability <= 0.5)
			{
				(*stability_color)[e.index_] = Vec3(0, stability, (0.5 - stability));
			}
			else
			{
				(*stability_color)[e.index_] = Vec3(stability-0.5, (1-stability),0);
			}
			return true;
		});
	}

	 void collapse_non_manifold_using_QMat(NONMANIFOLD& nm, uint32 number_vertices_erase)
	{
		using EdgeQueue = std::multimap<Scalar, NonManifoldEdge>;
		using EdgeQueueIt = typename EdgeQueue::const_iterator;
		using EdgeInfo = std::pair<bool, EdgeQueueIt>; // {valid, iterator}
		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;

		uint32 count = 0;
		EdgeQueue queue;
		auto edge_queue_it = add_attribute<EdgeInfo, NonManifoldEdge>(nm, "__non_manifold_edge_queue_it");
		auto position = get_attribute<Vec3, NonManifoldVertex>(nm, "position");
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(nm, "sphere_radius");
		auto sphere_info = add_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info");
		foreach_cell(nm,
					 [&](NonManifoldVertex v) {
						 Vec3 p = value<Vec3>(nm, position, v);
			value<Vec4>(nm, sphere_info, v) = Vec4(p[0], p[1], p[2], value<double>(nm,sphere_radius,v));
			return true;
		});
		QMatHelper helper(nm, sphere_info.get()); 

		 foreach_cell(nm, [&](NonManifoldEdge e) -> bool {
			Vec4 sphere_opt = helper.edge_optimal(e);
			value<EdgeInfo>(nm, edge_queue_it, e) = {
				true, 
				queue.emplace(helper.edge_cost(e,sphere_opt),e)};
			return true;
		});

		 while (!queue.empty() && count < number_vertices_erase)
		{
			auto it = queue.begin();
			NonManifoldEdge e = (*it).second;
			queue.erase(it);
			value<EdgeInfo>(nm, edge_queue_it, e).first = false;

			std::vector<NonManifoldFace> ifaces = incident_faces(nm, e);
			//TO do : add cone error to collapse edge
			if (ifaces.size() == 0)
				continue;

			auto [v, removed_edges] = collapse_edge(nm, e);
			Vec4 opt = helper.edge_optimal(e);
			value<Vec4>(nm, sphere_info, v) = opt;
			for (NonManifoldEdge re : removed_edges)
			{
				EdgeInfo einfo = value<EdgeInfo>(nm, edge_queue_it, re);
				if (einfo.first)
					queue.erase(einfo.second);
			}

			foreach_incident_edge(nm, v, [&](NonManifoldEdge ie) -> bool {
				EdgeInfo einfo = value<EdgeInfo>(nm, edge_queue_it, ie);
				if (einfo.first)
					queue.erase(einfo.second);
				value<EdgeInfo>(nm, edge_queue_it, ie) = {true, queue.emplace(helper.edge_cost(ie, opt), e)};
				return true;
			});
			++count;
		}

		foreach_cell(nm, [&](NonManifoldVertex v) -> bool {
			value<Vec3>(nm, position, v) = value<Vec4>(nm, sphere_info, v).head<3>();
			return true;
		});

		remove_attribute<NonManifoldEdge>(nm, edge_queue_it);
		remove_attribute<NonManifoldVertex>(nm, sphere_info);

		nonmanifold_provider_->emit_connectivity_changed(nm);
		nonmanifold_provider_->emit_attribute_changed(nm, position.get());
	}

	 void collapse_non_manifold_using_stability_ratio(NONMANIFOLD& nm, uint32 number_vertices_erase)
	{
		using EdgeQueue = std::multimap<Scalar, NonManifoldEdge>;
		using EdgeQueueIt = typename EdgeQueue::const_iterator;
		using EdgeInfo = std::pair<bool, EdgeQueueIt>; // {valid, iterator}

		auto position = get_attribute<Vec3, NonManifoldVertex>(nm, "position");
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(nm, "sphere_radius");
		auto stability_color = get_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		uint32 count = 0;
		EdgeQueue queue;
		auto edge_queue_it = add_attribute<EdgeInfo, NonManifoldEdge>(nm, "__non_manifold_edge_queue_it");

		foreach_cell(nm, [&](NonManifoldEdge e) -> bool{
			value<EdgeInfo>(nm, edge_queue_it, e) = {true, queue.emplace(value<double>(nm, stability_ratio_, e), e)};
			return true;
		});

		while (!queue.empty() && count < number_vertices_erase)
		{
			auto it = queue.begin();
			NonManifoldEdge e = (*it).second;
			queue.erase(it);
			value<EdgeInfo>(nm, edge_queue_it, e).first = false;
			if (value<double>(nm, stability_ratio_, e) > 0.90)
				break;
			auto [v, removed_edges] = collapse_edge(nm, e);
			//remove_edge_stability(nm, e);
			for (NonManifoldEdge re : removed_edges)
			{
				EdgeInfo einfo = value<EdgeInfo>(nm, edge_queue_it, re);
				if (einfo.first)
					queue.erase(einfo.second);
			}
// 			foreach_incident_edge(nm, v, [&](NonManifoldEdge ie) -> bool {
// 				EdgeInfo einfo = value<EdgeInfo>(nm, edge_queue_it, ie);
// 				if (einfo.first)
// 					queue.erase(einfo.second);
// 				std::vector<NonManifoldVertex> iv = incident_vertices(nm, ie);
// 				const Vec3& v1_p = value<Vec3>(nm, position, iv[0]);
// 				const Vec3& v2_p = value<Vec3>(nm, position, iv[1]);
// 				const double& r1 = value<double>(nm, sphere_radius, iv[0]);
// 				const double& r2 = value<double>(nm, sphere_radius, iv[1]);
// 				const double center_dist = (v1_p - v2_p).norm();
// 				double dis = std::max(0.0, (center_dist - std::abs(r1 - r2)));
// 				double stability = dis / center_dist;
// 				// std::cout << "Edge: " << e.index_ << ", Stability ratio: " << stability << std::flush;
// 				value<double>(nm, stability_ratio,ie) = stability;
// 				if (stability <= 0.5)
// 				{
// 					value<Vec3>(nm, stability_color, ie) = Vec3(0, stability, (0.5 - stability));
// 				}
// 				else
// 				{
// 					value<Vec3>(nm, stability_color, ie) = Vec3(stability - 0.5, (1 - stability), 0);
// 				}
// 				value<EdgeInfo>(nm, edge_queue_it, ie) = {
// 					true, queue.emplace(stability, ie)};
// 				return true;
// 			});
			++count;
		}

		remove_attribute<NonManifoldEdge>(nm, edge_queue_it);

		nonmanifold_provider_->emit_connectivity_changed(nm);
		nonmanifold_provider_->emit_attribute_changed(nm, position.get());
// 		nonmanifold_provider_->emit_attribute_changed(nm, stability_ratio.get());
// 		nonmanifold_provider_->emit_attribute_changed(nm, stability_color.get());
	}

 	 void point_selection_by_coverage_axis(SURFACE& surface, NONMANIFOLD& mv, double dilation_factor)
	{

		auto inner_position = get_attribute<Vec3, NonManifoldVertex>(mv, "position");
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(mv, "sphere_radius");
		auto sample_position = get_attribute<Vec3, Vertex>(surface, "position");
		auto inner_point_nb = nb_cells<NonManifoldVertex>(mv);
		auto sample_point_nb = nb_cells<Vertex>(surface);
		Eigen::MatrixXi A(sample_point_nb, inner_point_nb);
		foreach_cell(mv, [&](NonManifoldVertex nv) {
			foreach_cell(surface, [&](Vertex v) {
				A(index_of(mv,nv),
				index_of(surface,v)) =
					inside_sphere(value<Vec3>(mv, inner_position, nv), value<Vec3>(surface, sample_position, v),
								  value<double>(mv, sphere_radius, nv) * (1 + dilation_factor));
				return true;
			});
			return true;
		});
		HighsLp lp;
		lp.num_col_ = A.cols();
		lp.num_row_ = A.rows();

		// Adding decision variable bounds
		HighsVarType type = HighsVarType::kInteger;
		lp.col_lower_ = std::vector<double>(lp.num_col_, 0.0);
		lp.col_upper_ = std::vector<double>(lp.num_col_, 1.0);
		lp.row_lower_ = std::vector<double>(lp.num_row_, 1.0);
		lp.row_upper_ = std::vector<double>(lp.num_row_, 1e30);
		lp.integrality_ = std::vector<HighsVarType>(lp.num_col_, type);

		lp.a_matrix_.num_col_ = lp.num_col_;
		lp.a_matrix_.num_row_ = lp.num_row_;
		int currentStart = 0;
		for (int col = 0; col < lp.a_matrix_.num_col_; ++col)
		{
			lp.a_matrix_.start_.push_back(currentStart);

			for (int row = 0; row < lp.a_matrix_.num_row_; ++row)
			{
				double value = A(row, col);
				if (value != 0.0)
				{
					lp.a_matrix_.index_.push_back(row);
					lp.a_matrix_.value_.push_back(value);
					currentStart++;
				}
			}
		}

		lp.col_cost_ = std::vector<double>(lp.num_col_, 1.0);

		Highs highs;
		HighsStatus status = highs.passModel(lp);

		if (status == HighsStatus::kOk)
		{
			highs.run();
			auto solutions = highs.getSolution();
		}
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
		imgui_mesh_selector(nonmanifold_provider_, selected_medial_axis, "Medial_axis", [&](NONMANIFOLD& nm) {
			selected_medial_axis = &nm;
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
			if (selected_medial_axis)
			{
				if (ImGui::Button("Compute stablility ratio"))
				{
					compute_stability_ratio(*selected_medial_axis);
				}
				static int32 percent_vertices_to_remove = 10;
				ImGui::SliderInt("% vertices to keep", &percent_vertices_to_remove, 1, 99);
				 if (ImGui::Button("QMAT"))
				{
					collapse_non_manifold_using_stability_ratio(*selected_medial_axis,
																0.01*percent_vertices_to_remove *
																	nb_cells<NonManifoldVertex>(*selected_medial_axis));
					
				}
				if (ImGui::Button("Coverage Axis"))
				{
					static float dilation_factor = 0.1f;
					ImGui::SliderFloat("Dilation factor", &dilation_factor, 0.0f, 1.0f, "%.4f");
					point_selection_by_coverage_axis(*selected_surface_mesh_, * selected_medial_axis, dilation_factor );
				}
				
			}
		}
	}

private:
	POINT* surface_sample = nullptr;
	SURFACE* selected_surface_mesh_ = nullptr;
	NONMANIFOLD* selected_medial_axis = nullptr;
	std::shared_ptr<NonManifoldAttribute<double>> stability_ratio_ = nullptr;
	MeshProvider<POINT>* point_provider_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
	Regular medial_axis;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POWER_SHAPE_H_
