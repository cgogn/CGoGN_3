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

#include <cgogn/geometry/types/vector_traits.h>

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

namespace cgogn
{

namespace ui
{

template <typename SURFACE, typename NONMANIFOLD>
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

	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NonManifoldFace = typename mesh_traits<NONMANIFOLD>::Face;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

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
	struct Pole
	{
		Point _p;
		Pole(Point& p)
		{
			_p = p;
		}

		bool operator==(const Pole& other) const
		{
			return _p.x() == other._p.x() && _p.y() == other._p.y() && _p.z() == other._p.z();
		}
	};

	struct pole_hash
	{
		std::size_t operator()(const Pole& p) const
		{
			return ((std::hash<double>()(p._p.x()) ^ (std::hash<double>()(p._p.y()) << 1)) >> 1) ^
				   (std::hash<double>()(p._p.z()) << 1);
		}
	};

	struct edge_hash
	{
		std::size_t operator()(const std::pair<uint32, uint32>& edge) const
		{
			return std::hash<uint32>()(edge.first) + std::hash<uint32>()(edge.second);
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
	void compute_power_shape(SURFACE& surface)
	{
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_power_shape");
		
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

		std::vector<Point> Delaunay_tri_point;
		std::vector<Weight_Point> Power_point;

		// 1. Compute the Voronoi diagram of the sample points S.
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		// Compute the Bounding box of mesh
		std::array<Point, 8> obb_points;
		Point acc(0, 0, 0);
		CGAL::oriented_bounding_box(csm, obb_points, CGAL::parameters::use_convex_hull(true));
		for (size_t i = 0; i < obb_points.size(); i++)
		{
			acc += K::Vector_3(obb_points[i].x(), obb_points[i].y(), obb_points[i].z());
		}
		std::array<double, 3> center{acc.x() / 8, acc.y() / 8, acc.z() / 8}; //TODO
		// Create a large box surrounding object so that the Voronoi vertices are bounded
		double offset = 2.0;//TODO
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
			// CGAL::parameters::use_monte_carlo_sampling(true).number_of_points_per_area_unit(0.01));
			CGAL::parameters::use_grid_sampling(true).grid_spacing(0.1));

		// 	Add bounding box vertices in the sample points set
		for (auto& p : cube_corners)
		{
			Delaunay_tri_point.emplace_back(p[0], p[1], p[2]);
		}

		// Add sampled vertices into the volume data to construct the delauney tredrahedron
		for (auto& s : mesh_samples)
		{
			Delaunay_tri_point.emplace_back(s[0], s[1], s[2]);
		}

		uint32 nb_vertices = Delaunay_tri_point.size();

//		auto start_timer = std::chrono::high_resolution_clock::now();

		// Indices info for constructing volume data in Cgogn
		std::vector<unsigned> indices;
		indices.reserve(nb_vertices);
		for (unsigned int i = 0; i < nb_vertices; ++i)
			indices.push_back(i);

		// Construct delauney tredrahedron using CGAL
		Delaunay tri(boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.begin(), indices.begin())),
					 boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.end(), indices.end())));

// 		auto end_timer = std::chrono::high_resolution_clock::now();
// 		std::chrono::duration<double> elapsed_seconds = end_timer - start_timer;
// 		std::cout << "CGAL delaunay for " << nb_vertices << " points in " << elapsed_seconds.count() << std::endl;

//		start_timer = std::chrono::high_resolution_clock::now();

		// Indices of poles and a bool to indicate if it is inside of the object
		std::vector<std::pair<uint32, bool>> power_indices;
		// Find inside and outside poles
		uint32 count = 0, inside_poles_count = 0;
		Point farthest_inside_point, farthest_outside_point;
		double farthest_inside_distance, farthest_outside_distance;
		std::unordered_set<Pole, pole_hash> poles;
		std::unordered_map<uint32, uint32> inside_poles_indices; // Use for the construction of medial axis
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			Delaunay::Vertex_handle v = vit;
			std::vector<Delaunay::Facet> facets;
			farthest_inside_distance = 0;
			farthest_outside_distance = 0;
			tri.finite_incident_facets(vit, std::back_inserter(facets));
			if (facets.size() != 0)
			{
				for (auto f = facets.begin(); f != facets.end(); ++f)
				{
					K::Object_3 o = tri.dual(*f);
					if (const K::Segment_3* s = CGAL::object_cast<K::Segment_3>(&o))
					{
						auto start = s->start(), source = s->source();
						if (f == facets.begin())
						{
							auto dis = CGAL::squared_distance(start, v->point());
							if (pointInside(tree, start) && dis > farthest_inside_distance)
							{
								farthest_inside_point = start;
								farthest_inside_distance = dis;
							}
							else if (!pointInside(tree, start) && dis > farthest_outside_distance)
							{
								farthest_outside_point = start;
								farthest_outside_distance = dis;
							}
						}
						auto dis = CGAL::squared_distance(source, v->point());
						if (pointInside(tree, source) && dis > farthest_inside_distance)
						{
							farthest_inside_point = source;
							farthest_inside_distance = dis;
						}
						else if (!pointInside(tree, source) && dis > farthest_outside_distance)
						{
							farthest_outside_point = source;
							farthest_outside_distance = dis;
						}
					}
				}
			}
			if (farthest_inside_distance != 0)
			{
				Pole inside_pole = Pole(farthest_inside_point);
				if (poles.find(inside_pole) == poles.end())
				{
					poles.insert(inside_pole);
					Power_point.push_back(Weight_Point(farthest_inside_point, farthest_inside_distance));
					power_indices.push_back({count, true});
					inside_poles_indices.insert({count++, inside_poles_count++});
				}
			}
			if (farthest_outside_distance != 0)
			{
				Pole outside_pole = Pole(farthest_outside_point);
				if (poles.find(outside_pole) == poles.end())
				{
					poles.insert(outside_pole);
					Power_point.push_back(Weight_Point(farthest_outside_point, farthest_outside_distance));
					power_indices.push_back({count++, false});
				}
			}
		}
		// Build weighted delaunay tredrahedron

		cgogn::io::IncidenceGraphImportData Power_shape_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash> edge_indices;
		uint32 edge_count = 0;

		Regular reg(boost::make_zip_iterator(boost::make_tuple(Power_point.begin(), power_indices.begin())),
					boost::make_zip_iterator(boost::make_tuple(Power_point.end(), power_indices.end())));
		// Add vertices
		for (size_t idx = 0; idx < Power_point.size(); ++idx)
		{
			if (power_indices[idx].second)
				Power_shape_data.vertex_position_.push_back(
					{Power_point[idx].x(), Power_point[idx].y(), Power_point[idx].z()});
		}

		bool inside;
		uint32 v, v1, v2;
		int mask;
		for (auto fit = reg.finite_facets_begin(); fit != reg.finite_facets_end(); ++fit)
		{
			inside = true;
			v = fit->second;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != v)
				{
					inside &= fit->first->vertex(idx)->info().second;
				}
			}
			if (inside)
			{
				Power_shape_data.faces_nb_edges_.push_back(3);

				for (size_t i = 0; i < 4; i++)
				{
					if (i != v)
					{
						for (size_t j = i + 1; j < 4; j++)
						{
							if (j != v)
							{

								// Add edge
								v1 = inside_poles_indices[fit->first->vertex(i)->info().first];
								v2 = inside_poles_indices[fit->first->vertex(j)->info().first];
								if (edge_indices.find({v1, v2}) == edge_indices.end())
								{
									edge_indices.insert({{v1, v2}, edge_count});
									Power_shape_data.edges_vertex_indices_.push_back(v1);
									Power_shape_data.edges_vertex_indices_.push_back(v2);
									edge_count++;
								}
								// Add face
								Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v1, v2}]);
							}
						}
					}
				}
			}
		}
		uint32 power_nb_vertices = Power_shape_data.vertex_position_.size();
		uint32 power_nb_edges = Power_shape_data.edges_vertex_indices_.size() / 2;
		uint32 power_nb_faces = Power_shape_data.faces_nb_edges_.size();

		Power_shape_data.set_parameter(power_nb_vertices, power_nb_edges, power_nb_faces);
		

		import_incidence_graph_data(*mv, Power_shape_data);

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position = get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}

protected:
	void init() override
	{
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

		if (selected_surface_mesh_)
		{
			if (ImGui::Button("Power shape"))
				compute_power_shape(*selected_surface_mesh_);
		}
	}

private:
	SURFACE* selected_surface_mesh_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POWER_SHAPE_H_
