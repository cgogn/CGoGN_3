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

#include <cgogn/modeling/algos/graph_to_hex.h>

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/halfedge.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/io/surface/surface_import.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/modeling/algos/subdivision.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/fitting.h>
#include <cgogn/geometry/functions/intersection.h>
#include <cgogn/geometry/functions/projection.h>

#include <fstream>
#include <iostream>

#include <cmath>
#include <numeric>

namespace cgogn
{

namespace modeling
{

///////////
// CMaps //
///////////

std::tuple<GAttributes, M2Attributes, M3Attributes> graph_to_hex(Graph& g, CMap2& m2, CMap3& m3)
{
	bool okay;

	GraphData gData;
	GAttributes gAttribs;
	M2Attributes m2Attribs;
	M3Attributes m3Attribs;

	okay = get_graph_data(g, gData);
	std::cout << uint32(gData.intersections.size()) << " intersections" << std::endl;
	std::cout << uint32(gData.branches.size()) << " branches" << std::endl;
	// for (auto b : gData.branches)
	// {
	// 	uint32 nb_edges = 0;
	// 	Dart cur = b.first.dart;
	// 	Dart next = alpha0(g, cur);
	// 	do
	// 	{
	// 		++nb_edges;
	// 		cur = alpha1(g, next);
	// 		next = alpha0(g, cur);
	// 	} while (alpha_1(g, cur) != b.second.dart);
	// 	std::cout << index_of(g, Graph::Vertex(b.first.dart)) << " - " << index_of(g, Graph::Vertex(b.second.dart))
	// 			  << " / nb_edges: " << nb_edges << std::endl;
	// }
	if (!okay)
		std::cout << "error graph_to_hex: get_graph_data" << std::endl;
	else
		std::cout << "graph_to_hex (/): got graph data" << std::endl;

	if (okay)
		okay = add_graph_attributes(g, gAttribs);
	if (!okay)
		std::cout << "error graph_to_hex: add_graph_attributes" << std::endl;
	else
		std::cout << "graph_to_hex (/): added graph attributes" << std::endl;

	if (okay)
		okay = add_cmap2_attributes(m2, m2Attribs);
	if (!okay)
		std::cout << "error graph_to_hex: add_cmap2_attributes" << std::endl;
	else
		std::cout << "graph_to_hex (/): added cmap2 attributes" << std::endl;

	if (okay)
		okay = build_contact_surfaces(g, gAttribs, m2, m2Attribs);
	if (!okay)
		std::cout << "error graph_to_hex: build_contact_surfaces" << std::endl;
	else
		std::cout << "graph_to_hex (/): contact surfaces built" << std::endl;

	if (okay)
		okay = create_intersection_frames(g, gAttribs, m2, m2Attribs);
	if (!okay)
		std::cout << "error graph_to_hex: create_intersections_frames" << std::endl;
	else
		std::cout << "graph_to_hex (/): create_intersections_frames completed" << std::endl;

	if (okay)
		okay = propagate_frames(g, gAttribs, gData, m2);
	if (!okay)
		std::cout << "error graph_to_hex: propagate_frames" << std::endl;
	else
		std::cout << "graph_to_hex (/): propagate_frames completed" << std::endl;

	if (okay)
		okay = set_contact_surfaces_geometry(g, gAttribs, m2, m2Attribs);
	if (!okay)
		std::cout << "error graph_to_hex: set_contact_surfaces_geometry" << std::endl;
	else
		std::cout << "graph_to_hex (/): set_contact_surfaces_geometry completed" << std::endl;

	if (okay)
		okay = build_branch_sections(g, gAttribs, m2, m2Attribs, m3);
	if (!okay)
		std::cout << "error graph_to_hex: build_branch_sections" << std::endl;
	else
		std::cout << "graph_to_hex (/): build_branch_sections completed" << std::endl;

	if (okay)
		okay = sew_branch_sections(m2, m2Attribs, m3);
	if (!okay)
		std::cout << "error graph_to_hex: sew_sections" << std::endl;
	else
		std::cout << "graph_to_hex (/): sew_sections completed" << std::endl;

	if (okay)
	{
		add_cmap3_attributes(m3, m3Attribs);
		okay = set_volumes_geometry(m2, m2Attribs, m3, m3Attribs);
	}
	if (!okay)
		std::cout << "error graph_to_hex: set_volumes_geometry" << std::endl;
	else
		std::cout << "graph_to_hex (/): set_volumes_geometry completed" << std::endl;

	// bloat(m3, g, gAttribs);

	// uint32 nb_vertices_b = 0;
	// uint32 nb_edges_b = 0;
	// uint32 nb_faces_b = 0;
	// foreach_cell(m3, [&](CMap3::Edge e3) -> bool {
	// 	if (is_incident_to_boundary(m3, e3))
	// 		++nb_edges_b;
	// 	return true;
	// });
	// foreach_cell(m3, [&](CMap3::Face e3) -> bool {
	// 	if (is_incident_to_boundary(m3, e3))
	// 		++nb_faces_b;
	// 	return true;
	// });
	// foreach_cell(m3, [&](CMap3::Vertex e3) -> bool {
	// 	if (is_incident_to_boundary(m3, e3))
	// 		++nb_vertices_b;
	// 	std::cout << degree(m3, e3) << std::endl;
	// 	return true;
	// });

	// std::cout << "boundary cells: " << nb_vertices_b << " : " << nb_edges_b << " : " << nb_faces_b << std::endl;

	// parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
	// 	CMap2::Volume contact_surface(value<Dart>(g, gAttribs.vertex_contact_surface, v));

	// 	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	// 	// if (degree(g, v) == 1)
	// 	// 	value<Vec3>(m2, m2Attribs.volume_center, contact_surface) =
	// 	// 		center + (0.25 * (center - value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, v.dart)))));
	// 	// else
	// 	value<Vec3>(m2, m2Attribs.volume_center, contact_surface) = center;

	// 	Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v) * 1.1;

	// 	if (degree(g, v) < 3)
	// 	{
	// 		Graph::HalfEdge h(v.dart);
	// 		Dart csf = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h);
	// 		Mat3 frame = value<Mat3>(g, gAttribs.halfedge_frame, h);

	// 		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(csf)) = center - frame.col(1) * radius;
	// 		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, csf))) = center + frame.col(0) * radius;
	// 		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi<1, 1>(m2, csf))) =
	// 			center + frame.col(1) * radius;
	// 		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, csf))) = center - frame.col(0) * radius;
	// 	}
	// 	else
	// 	{
	// 		foreach_incident_vertex(m2, contact_surface, [&](CMap2::Vertex v2) -> bool {
	// 			Vec3 pos = value<Vec3>(m2, m2Attribs.vertex_position, v2);
	// 			geometry::project_on_sphere(pos, center, radius);
	// 			value<Vec3>(m2, m2Attribs.vertex_position, v2) = pos;
	// 			return true;
	// 		});
	// 	}

	// 	foreach_incident_edge(m2, contact_surface, [&](CMap2::Edge e) -> bool {
	// 		std::vector<CMap2::Vertex> vertices = incident_vertices(m2, e);
	// 		Vec3 mid = 0.5 * (value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) +
	// 						  value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]));
	// 		geometry::project_on_sphere(mid, center, radius);
	// 		value<Vec3>(m2, m2Attribs.edge_mid, e) = mid;
	// 		return true;
	// 	});
	// 	return true;
	// });

	return {gAttribs, m2Attribs, m3Attribs};
}

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

void index_volume_cells(CMap2& m, CMap2::Volume vol)
{
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			if (index_of(m, CMap2::HalfEdge(d)) == INVALID_INDEX)
				set_index(m, CMap2::HalfEdge(d), new_index<CMap2::HalfEdge>(m));
		}
		if (is_indexed<CMap2::Vertex>(m))
		{
			if (index_of(m, CMap2::Vertex(d)) == INVALID_INDEX)
				set_index(m, CMap2::Vertex(d), new_index<CMap2::Vertex>(m));
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			if (index_of(m, CMap2::Edge(d)) == INVALID_INDEX)
				set_index(m, CMap2::Edge(d), new_index<CMap2::Edge>(m));
		}
		if (is_indexed<CMap2::Face>(m))
		{
			if (index_of(m, CMap2::Face(d)) == INVALID_INDEX)
				set_index(m, CMap2::Face(d), new_index<CMap2::Face>(m));
		}
		return true;
	});
	if (is_indexed<CMap2::Volume>(m))
	{
		if (index_of(m, vol) == INVALID_INDEX)
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}
}

void sew_volumes(CMap3& m, Dart d0, Dart d1)
{
	cgogn_message_assert(codegree(m, CMap3::Face(d0)) == codegree(m, CMap3::Face(d1)),
						 "The faces to sew do not have the same codegree");
	Dart it0 = d0;
	Dart it1 = d1;
	do
	{
		cgogn_message_assert(phi3(m, it0) == it0 && phi3(m, it1) == it1, "The faces to sew are already sewn");
		phi3_sew(m, it0, it1);
		it0 = phi1(m, it0);
		it1 = phi_1(m, it1);
	} while (it0 != d0);
}

void unsew_volumes(CMap3& m, Dart d)
{
	Dart it = d;
	do
	{
		cgogn_message_assert(phi3(m, it) != it, "The faces to unsew are already not sewn");
		phi3_unsew(m, it);
		it = phi1(m, it);
	} while (it != d);
}

Dart add_branch_section(CMap3& m3)
{
	std::vector<Dart> D = {add_prism(static_cast<CMap2&>(m3), 4).dart, add_prism(static_cast<CMap2&>(m3), 4).dart,
						   add_prism(static_cast<CMap2&>(m3), 4).dart, add_prism(static_cast<CMap2&>(m3), 4).dart};

	sew_volumes(m3, phi2(m3, D[0]), phi2(m3, phi_1(m3, (D[1]))));
	sew_volumes(m3, phi2(m3, D[1]), phi2(m3, phi_1(m3, (D[2]))));
	sew_volumes(m3, phi2(m3, D[2]), phi2(m3, phi_1(m3, (D[3]))));
	sew_volumes(m3, phi2(m3, D[3]), phi2(m3, phi_1(m3, (D[0]))));

	return D[0];
}

void shift_frame(Mat3& frame, uint32 nb_shifts)
{
	for (uint32 i = 0; i < nb_shifts; ++i)
	{
		Vec3 R = frame.col(1);
		Vec3 S = -frame.col(0);
		frame.col(0) = R;
		frame.col(1) = S;
	}
}

void dualize_volume(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs, const Graph& g, GAttributes& gAttribs)
{
	// set the new phi1
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		Dart dd = phi2(m, phi_1(m, d));
		(*(m.phi1_))[d.index] = dd;
		return true;
	});

	// set the new phi_1
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		Dart next = phi1(m, d);
		(*(m.phi_1_))[next.index] = d;
		return true;
	});

	DartMarkerStore<CMap2> face_marker(m);
	const Graph::Vertex gv = value<Graph::Vertex>(m, m2Attribs.volume_gvertex, vol);
	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, gv);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (!face_marker.is_marked(d))
		{
			CMap2::Face f(d);
			foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
				face_marker.mark(d);
				return true;
			});
			if (is_indexed<CMap2::Face>(m))
				set_index(m, f, new_index<CMap2::Face>(m)); // give a new index to the face
			// darts of the new face orbit are darts of the old dual vertex
			// get the outgoing graph branch from the old dual vertex
			Dart branch = value<Dart>(m, m2Attribs.dual_vertex_graph_branch, CMap2::Vertex(d));
			// store the contact surface face on the outgoing branch
			value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(branch)) = d;
			return true;
		}
		return true;
	});

	DartMarkerStore<CMap2> vertex_marker(m);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (!vertex_marker.is_marked(d))
		{
			CMap2::Vertex v(d);

			std::vector<Vec3> points;
			foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
				vertex_marker.mark(d);
				points.push_back(value<Vec3>(m, m2Attribs.vertex_position, CMap2::Vertex(d)) - center);
				return true;
			});
			Vec3 b;
			if (points.size() == 2)
				b = slerp(points[0], points[1], 0.5, true) + center;
			else
				b = spherical_barycenter(points, 10) + center;

			set_index(m, v, new_index<CMap2::Vertex>(m));	  // give a new index to the vertex
			value<Vec3>(m, m2Attribs.vertex_position, v) = b; // set the position to the computed position
			return true;
		}
		return true;
	});
}

// void bloat(CMap3& m3, Scalar s)
void bloat(CMap3& m3, const Graph& g, const GAttributes& gAttribs)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position3 = get_attribute<Vec3, CMap3::Vertex>(m3, "position");

	foreach_cell(g, [&](Graph::Vertex gv) -> bool {
		if (degree(g, gv) == 1)
		{
			Vec3 p0 = value<Vec3>(g, gAttribs.vertex_position, gv);
			std::vector<Graph::Vertex> neigh = adjacent_vertices_through_edge(g, gv);
			Vec3 p1 = value<Vec3>(g, gAttribs.vertex_position, neigh[0]);
			Vec3 offset = (p0 - p1).normalized() * value<Scalar>(g, gAttribs.vertex_radius, gv);
			CMap3::Vertex v3 =
				CMap3::Vertex(value<Dart>(g, gAttribs.halfedge_volume_connection, Graph::HalfEdge(gv.dart)));
			value<Vec3>(m3, vertex_position3, v3) += offset;
		}
		return true;
	});

	// foreach_cell(m3, [&](CMap3::Vertex v0) -> bool {
	// if (!is_incident_to_boundary(m3, v0))
	//{
	//	Vec3 p0 = value<Vec3>(m3, vertex_position3, v0);
	//	foreach_adjacent_vertex_through_edge(m3, v0, [&](CMap3::Vertex v1) -> bool {
	//		if ()
	//		return true;
	//	});
	//}
	// return true;
	//});
}

void catmull_clark_approx(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, uint32 iterations)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_delta = add_attribute<Vec3, CMap2::Vertex>(m, "delta");
	std::shared_ptr<CMap2::Attribute<Vec3>> incident_faces_mid =
		add_attribute<Vec3, CMap2::Vertex>(m, "incident_faces_mid");

	for (uint32 i = 0; i < iterations; ++i)
	{
		CellCache<CMap2> cache_init_vertices(m);
		cache_init_vertices.template build<CMap2::Vertex>();
		CellCache<CMap2> cache_edge_vertices(m);
		CellCache<CMap2> cache_face_vertices(m);

		cgogn::modeling::quadrangulate_all_faces(
			m,
			[&](CMap2::Vertex v) -> void {
				cache_edge_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				value<Vec3>(m, incident_faces_mid, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += 0.5 * value<Vec3>(m, vertex_position, v2);
					value<Vec3>(m, vertex_delta, v) -= 0.25 * value<Vec3>(m, vertex_position, v2);
					return true;
				});
			},
			[&](CMap2::Vertex v) -> void {
				cache_face_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_position, v2);
					return true;
				});
				value<Vec3>(m, vertex_position, v) /= degree(m, v);

				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_delta, v2) += 0.25 * value<Vec3>(m, vertex_position, v);
					value<Vec3>(m, incident_faces_mid, v2) += 0.5 * value<Vec3>(m, vertex_position, v);
					return true;
				});
			});

		foreach_cell(cache_init_vertices, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
			Vec3 F = {0, 0, 0};
			Vec3 R = {0, 0, 0};
			Scalar n = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				F += value<Vec3>(m, incident_faces_mid, v2);
				R += value<Vec3>(m, vertex_position, v2);
				n += 1;
				return true;
			});

			value<Vec3>(m, vertex_delta, v) = (-3 * n * value<Vec3>(m, vertex_position, v) + F + 2 * R) / (n * n);
			return true;
		});

		foreach_cell(m, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_delta, v);
			return true;
		});
	}

	remove_attribute<CMap2::Vertex>(m, vertex_delta);
	remove_attribute<CMap2::Vertex>(m, incident_faces_mid);
}

void catmull_clark_inter(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, uint32 iterations)
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position2 = add_attribute<Vec3, CMap2::Vertex>(m, "position2");
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_delta = add_attribute<Vec3, CMap2::Vertex>(m, "delta");
	std::shared_ptr<CMap2::Attribute<Vec3>> incident_faces_mid =
		add_attribute<Vec3, CMap2::Vertex>(m, "incident_faces_mid");
	std::shared_ptr<CMap2::Attribute<Vec3>> incident_init_mid =
		add_attribute<Vec3, CMap2::Vertex>(m, "incident_init_mid");

	for (uint32 i = 0; i < iterations; ++i)
	{
		CellCache<CMap2> cache_init_vertices(m);
		cache_init_vertices.template build<CMap2::Vertex>();
		CellCache<CMap2> cache_edge_vertices(m);
		CellCache<CMap2> cache_face_vertices(m);

		foreach_cell(m, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
			return true;
		});

		cgogn::modeling::quadrangulate_all_faces(
			m,
			[&](CMap2::Vertex v) -> void {
				cache_edge_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				value<Vec3>(m, incident_faces_mid, v) = {0, 0, 0};
				value<Vec3>(m, incident_init_mid, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += 0.5 * value<Vec3>(m, vertex_position, v2);
					return true;
				});
			},
			[&](CMap2::Vertex v) -> void {
				cache_face_vertices.add(v);
				value<Vec3>(m, vertex_position, v) = {0, 0, 0};
				value<Vec3>(m, vertex_delta, v) = {0, 0, 0};
				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_position, v2);
					return true;
				});
				value<Vec3>(m, vertex_position, v) /= degree(m, v);

				foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
					value<Vec3>(m, vertex_delta, v2) += 0.25 * value<Vec3>(m, vertex_position, v);
					value<Vec3>(m, incident_faces_mid, v2) += 0.5 * value<Vec3>(m, vertex_position, v);
					return true;
				});
			});

		foreach_cell(cache_init_vertices, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, incident_faces_mid, v) = {0, 0, 0};
			Scalar n = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				value<Vec3>(m, incident_faces_mid, v) += value<Vec3>(m, incident_faces_mid, v2);
				n += 1;
				return true;
			});
			value<Vec3>(m, incident_faces_mid, v) /= n;

			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				value<Vec3>(m, vertex_delta, v2) -= 0.25 * value<Vec3>(m, incident_faces_mid, v);
				value<Vec3>(m, incident_init_mid, v2) += 0.5 * value<Vec3>(m, incident_faces_mid, v);
				return true;
			});
			return true;
		});

		foreach_cell(cache_face_vertices, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_delta, v);
			Scalar n = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap2::Vertex v2) -> bool {
				value<Vec3>(m, vertex_delta, v) -= 2 * value<Vec3>(m, incident_faces_mid, v2);
				value<Vec3>(m, vertex_delta, v) -= value<Vec3>(m, incident_init_mid, v2);
				n += 1;
				return true;
			});
			value<Vec3>(m, vertex_delta, v) += 3 * n * value<Vec3>(m, vertex_position, v);
			value<Vec3>(m, vertex_delta, v) /= n * n;
			return true;
		});

		foreach_cell(m, [&](CMap2::Vertex v) -> bool {
			value<Vec3>(m, vertex_position, v) += value<Vec3>(m, vertex_delta, v);
			return true;
		});
	}

	remove_attribute<CMap2::Vertex>(m, vertex_position2);
	remove_attribute<CMap2::Vertex>(m, vertex_delta);
	remove_attribute<CMap2::Vertex>(m, incident_faces_mid);
	remove_attribute<CMap2::Vertex>(m, incident_init_mid);
}

void padding(CMap3& m3, DartMarker<CMap3>* face_marker)
{
	Dart d_boundary;
	for (Dart d = m3.begin(), end = m3.end(); d != end && d_boundary.is_nil(); d = m3.next(d))
		if (is_boundary(m3, d))
			d_boundary = d;

	std::vector<CMap3::Volume> prisms;

	// visit all faces of the boundary and check if a prism must be inserted
	DartMarker visited_face(m3);
	foreach_dart_of_orbit(m3, CMap3::Volume(d_boundary), [&](Dart db) -> bool {
		if (!visited_face.is_marked(db))
		{
			uint32 n = 0;
			foreach_dart_of_orbit(m3, CMap3::Face2(db), [&](Dart df) -> bool {
				++n;
				visited_face.mark(df);
				return true;
			});

			Dart db3 = phi3(m3, db);
			unsew_volumes(m3, db);
			if (!face_marker || !face_marker->is_marked(db3))
			{
				CMap3::Volume p = add_prism(static_cast<CMap2&>(m3), n);
				sew_volumes(m3, db, p.dart);
				sew_volumes(m3, db3, phi<2, 1, 1, 2>(m3, p.dart));
				prisms.push_back(p);
			}
		}
		return true;
	});

	// sew neighboring prisms together
	DartMarker visited_edge(m3);
	foreach_dart_of_orbit(m3, CMap3::Volume(d_boundary), [&](Dart db) -> bool {
		if (!visited_edge.is_marked(db))
		{
			Dart db2 = phi2(m3, db);
			visited_edge.mark(db);
			visited_edge.mark(db2);
			if (phi3(m3, db) != db && phi3(m3, db2) != db2)
				sew_volumes(m3, phi<3, 2>(m3, db), phi<3, 2>(m3, db2));
			phi3_unsew(m3, db);
			phi3_unsew(m3, db2);
		}
		return true;
	});

	// rebuild boundary
	remove_volume(static_cast<CMap2&>(m3), CMap2::Volume(d_boundary));
	close(m3);

	// set indices
	// CMap2 cells (HalfEdge, Vertex2, Edge2, Face2, Volume) already indiced in add_prism
	for (CMap3::Volume vol : prisms)
	{
		foreach_dart_of_orbit(m3, vol, [&](Dart d) -> bool {
			if (is_indexed<CMap3::Vertex>(m3))
			{
				if (index_of(m3, CMap3::Vertex(d)) == INVALID_INDEX)
				{
					uint32 id = INVALID_INDEX;
					foreach_dart_of_orbit(m3, CMap3::Vertex(d), [&](Dart dv) -> bool {
						if (!is_boundary(m3, dv))
							id = index_of(m3, CMap3::Vertex(dv));
						return (id == INVALID_INDEX);
					});
					set_index(m3, CMap3::Vertex(d), (id == INVALID_INDEX ? new_index<CMap3::Vertex>(m3) : id));
				}
			}
			if (is_indexed<CMap3::Edge>(m3))
			{
				if (index_of(m3, CMap3::Edge(d)) == INVALID_INDEX)
				{
					uint32 id = INVALID_INDEX;
					foreach_dart_of_orbit(m3, CMap3::Edge(d), [&](Dart de) -> bool {
						if (!is_boundary(m3, de))
							id = index_of(m3, CMap3::Edge(de));
						return (id == INVALID_INDEX);
					});
					set_index(m3, CMap3::Edge(d), (id == INVALID_INDEX ? new_index<CMap3::Edge>(m3) : id));
				}
			}
			if (is_indexed<CMap3::Face>(m3))
			{
				if (index_of(m3, CMap3::Face(d)) == INVALID_INDEX)
				{
					uint32 id = INVALID_INDEX;
					if (!is_boundary(m3, phi3(m3, d)))
						id = index_of(m3, CMap3::Face(phi3(m3, d)));
					set_index(m3, CMap3::Face(d), (id == INVALID_INDEX ? new_index<CMap3::Face>(m3) : id));
				}
			}
			return true;
		});
	}

	// set positions
	auto vertex_position = get_attribute<Vec3, CMap3::Vertex>(m3, "position");

	for (CMap3::Volume vol : prisms)
	{
		Dart it = vol.dart;
		do
		{
			value<Vec3>(m3, vertex_position, CMap3::Vertex(it)) =
				value<Vec3>(m3, vertex_position, CMap3::Vertex(phi<2, 1, 1>(m3, it)));
			it = phi1(m3, it);
		} while (it != vol.dart);
	}
}

/*****************************************************************************/
/* Fiber wise subdvision                                                          */
/*****************************************************************************/

CMap3::Edge find_fiber_dir(CMap3& m3, CMap3::Face f)
{
	Dart dir0 = f.dart;
	Dart dir1 = phi1(m3, dir0);

	uint32 counter = 0;
	Dart d0 = dir0, d1 = dir1;
	do
	{
		++counter;
		d0 = phi<3, 2, 3, 1, 1>(m3, d0);
		d1 = phi<3, 2, 3, 1, 1>(m3, d1);
	} while (d0 != dir0 && d1 != dir1);

	return CMap3::Edge(d0 == dir0 ? dir0 : dir1);
}

uint32 get_ring_size(CMap3& m3, CMap3::Edge e)
{
	Dart d0 = e.dart;
	uint32 n = 0;
	Dart d = d0;
	do
	{
		++n;
		d = phi<3, 2, 3, 1, 1>(m3, d);
	} while (d != d0);
	return n;
}

bool unchecked_ring(CMap3& m3, CMap3::Edge e, uint32 ring_size, CellMarker<CMap3, CMap3::Edge>& visited_edge)
{
	if (visited_edge.is_marked(e) /* || check_ring(m, e, ring_size)*/)
		return false; // already explored || not a ring

	Dart d0 = e.dart;
	uint32 n = 0;
	Dart d = d0;
	do
	{
		visited_edge.mark(CMap3::Edge(d));
		++n;
		d = phi<3, 2, 3, 1, 1>(m3, d);
	} while (d != d0 && n < ring_size);
	return (d == d0 && n == ring_size); // true if looped around ring, false otherwise
}

void cut_slice(CMap3& m3, CMap3::Attribute<Vec3>* vertex_position, CellCache<CMap3>& slice)
{
	cut_all_edges(slice, [&](CMap3::Vertex v) {
		std::vector<CMap3::Vertex> av = adjacent_vertices_through_edge(m3, v);
		cgogn::value<Vec3>(m3, vertex_position, v) =
			0.5 * (cgogn::value<Vec3>(m3, vertex_position, av[0]) + cgogn::value<Vec3>(m3, vertex_position, av[1]));
	});

	foreach_cell(slice, [&](CMap3::Face f) -> bool {
		Dart d0 = phi1(m3, f.dart);
		Dart d1 = phi_1(m3, phi_1(m3, f.dart));
		cgogn::cut_face(m3, CMap3::Vertex(d0), CMap3::Vertex(d1));
		return true;
	});

	foreach_cell(slice, [&](CMap3::Volume w) -> bool {
		Dart d0 = phi1(m3, w.dart);
		std::vector<Dart> path;
		Dart d1 = d0;
		do
		{
			path.push_back(d1);
			d1 = phi<1, 2, 1>(m3, d1);
		} while (d1 != d0);

		cgogn::cut_volume(m3, path);
		return true;
	});
}

CellCache<CMap3> get_slice(CMap3& m3, CMap3::Edge e)
{
	CellCache<CMap3> slice(m3);

	cgogn::ui::CellsSet<CMap3, CMap3::Volume> slice_volumes(m3, "slice_w");
	cgogn::ui::CellsSet<CMap3, CMap3::Face> slice_faces(m3, "slice_f");
	cgogn::ui::CellsSet<CMap3, CMap3::Edge> slice_edges(m3, "slice_e");

	std::vector<Dart> pending;
	pending.push_back(e.dart);

	for (uint32 i = 0; i < pending.size(); ++i)
	{
		Dart d0 = pending[i];
		if (!slice_volumes.contains(CMap3::Volume(d0)))
		{
			slice_volumes.select(CMap3::Volume(d0));
			Dart d1 = d0;
			do
			{
				slice_edges.select(CMap3::Edge(d1));
				slice_faces.select(CMap3::Face(d1));
				foreach_incident_volume(m3, CMap3::Edge(d1), [&](CMap3::Volume w) -> bool {
					if (!is_boundary(m3, w.dart) && !slice_volumes.contains(w))
						pending.push_back(w.dart);
					return true;
				});
				d1 = phi<1, 1, 2>(m3, d1);
			} while (d1 != d0);
		}
	}

	slice_volumes.foreach_cell([&](CMap3::Volume w) { slice.add(w); });
	slice_faces.foreach_cell([&](CMap3::Face f) { slice.add(f); });
	slice_edges.foreach_cell([&](CMap3::Edge e) { slice.add(e); });

	return slice;
}

void volume_fiber_spread(CMap3& m3, CellCache<CMap3>& surface_fibers, CellMarker<CMap3, CMap3::Edge>& edge_fibers)
{
	foreach_cell(surface_fibers, [&](CMap3::Edge e0) -> bool {
		CellCache<CMap3> slice = get_slice(m3, e0);
		foreach_cell(slice, [&](CMap3::Edge e1) -> bool {
			edge_fibers.mark(e1);
			return true;
		});
		return true;
	});
}

CellCache<CMap3> surface_fiber_spread(CMap3& m, CMap3::Edge e0)
{
	CellMarker<CMap3, CMap3::Edge> visited_edge(m);
	CellCache<CMap3> fibers_cache(m);
	uint32 ring_size = get_ring_size(m, e0);
	unchecked_ring(m, e0, ring_size, visited_edge);

	std::vector<CMap3::Edge> fibers;
	fibers.push_back(e0);
	for (uint32 i = 0; i < fibers.size(); ++i)
	{
		Dart d0 = fibers[i].dart;
		fibers_cache.add(fibers[i]);

		Dart d = d0;
		do
		{
			CMap3::Edge e1 = CMap3::Edge(phi<1, 3, 2, 3, 1>(m, d));
			CMap3::Edge e_1 = CMap3::Edge(phi<3, 1, 2, 1, 3>(m, d));
			if (unchecked_ring(m, e1, ring_size, visited_edge))
				fibers.push_back(e1);
			if (unchecked_ring(m, e_1, ring_size, visited_edge))
				fibers.push_back(e_1);
			d = phi<3, 2, 3, 1, 1>(m, d);
		} while (d != d0);
	}
	return fibers_cache;
}

void mark_mesh_fibers(CMap3& m3, CMap3::Edge e, CellMarker<CMap3, CMap3::Edge>& edge_fibers)
{
	CellCache<CMap3> surface_fibers = surface_fiber_spread(m3, e);
	volume_fiber_spread(m3, surface_fibers, edge_fibers);
}

void fiber_aligned_subdivision(CMap3& m, CellMarker<CMap3, CMap3::Edge>& fibers)
{
	// find all non fibrous edges -> they will be cut
	// find all non fibrous faces -> they will be quadrangulated
	// sort volumes between fibrous and non fibrous
	//		fibrous will be cut in 4
	//		non fibrous will be cut in 8
	//
	// cut all non fibrous edges
	// quadrangulate all non fibrous faces
	// "cut" all volumes by adding leaflets

	CellMarker<CMap3, CMap3::Edge> new_edges(m);
	CellMarker<CMap3, CMap3::Edge> new_vertices(m);

	CellCache<CMap3> edges_to_cut(m);
	edges_to_cut.template build<CMap3::Edge>([&](CMap3::Edge e) { return !fibers.is_marked(e); });

	CellCache<CMap3> faces_to_quad(m);
	CellCache<CMap3> faces_to_bi(m);
	faces_to_quad.template build<CMap3::Face>([&](CMap3::Face f) {
		bool fiber = false;
		foreach_incident_edge(m, f, [&](CMap3::Edge e) -> bool {
			fiber = fibers.is_marked(e);
			if (fiber)
				faces_to_bi.add(CMap3::Face(e.dart));
			return !(fiber);
		});
		return !(fiber);
	});

	CellCache<CMap3> volumes_to_oct(m);
	CellCache<CMap3> volumes_to_quad(m);

	foreach_cell(m, [&](CMap3::Volume w) {
		bool fiber = false;
		foreach_incident_edge(m, w, [&](CMap3::Edge we) -> bool {
			fiber = fibers.is_marked(we);
			if (fiber)
				volumes_to_quad.add(CMap3::Volume(we.dart));
			// volumes are grabbed by the fiber for orientation
			return !fiber;
		});
		if (!fiber)
			volumes_to_oct.add(w);
		return true;
	});

	auto vertex_position = get_attribute<Vec3, CMap3::Vertex>(m, "position");

	quadrangulate_all_faces(
		faces_to_quad,
		[&](CMap3::Vertex v) {
			std::vector<CMap3::Vertex> av = adjacent_vertices_through_edge(m, v);
			cgogn::value<Vec3>(m, vertex_position, v) =
				0.5 * (cgogn::value<Vec3>(m, vertex_position, av[0]) + cgogn::value<Vec3>(m, vertex_position, av[1]));
		},
		[&](CMap3::Vertex v) {
			Vec3 center;
			center.setZero();
			uint32 count = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](CMap3::Vertex av) -> bool {
				center += cgogn::value<Vec3>(m, vertex_position, av);
				++count;
				return true;
			});
			center /= Scalar(count);
			cgogn::value<Vec3>(m, vertex_position, v) = center;
		});

	foreach_cell(faces_to_bi, [&](CMap3::Face f) -> bool {
		CMap3::Vertex v0(phi<1, 1>(m, f.dart));
		CMap3::Vertex v1(phi_1(m, f.dart));
		cut_face(m, v0, v1);
		return true;
	});

	foreach_cell(volumes_to_oct, [&](CMap3::Volume w) -> bool {
		CMap3::Vertex v = octosect_hex(m, w);
		Vec3 center;
		center.setZero();
		uint32 count = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&](CMap3::Vertex av) -> bool {
			center += cgogn::value<Vec3>(m, vertex_position, av);
			++count;
			return true;
		});
		center /= Scalar(count);
		cgogn::value<Vec3>(m, vertex_position, v) = center;

		return true;
	});

	foreach_cell(volumes_to_quad, [&](CMap3::Volume w) -> bool {
		quadrisect_hex(m, w);
		return true;
	});
}

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool subdivide_graph(Graph& g)
{
	auto vertex_position = get_attribute<Vec3, Graph::Vertex>(g, "position");
	if (!vertex_position)
	{
		std::cout << "The graph has no vertex position attribute" << std::endl;
		return false;
	}

	auto vertex_radius = get_attribute<Scalar, Graph::Vertex>(g, "radius");
	if (!vertex_radius)
	{
		std::cout << "The graph has no vertex radius attribute" << std::endl;
		return false;
	}

	CellCache<Graph> cache(g);
	cache.template build<Graph::Edge>();

	foreach_cell(cache, [&](Graph::Edge eg) -> bool {
		const Vec3& P0 = value<Vec3>(g, vertex_position, Graph::Vertex(eg.dart));
		const Vec3& Pn = value<Vec3>(g, vertex_position, Graph::Vertex(alpha0(g, eg.dart)));
		const Scalar R0 = value<Scalar>(g, vertex_radius, Graph::Vertex(eg.dart));
		const Scalar Rn = value<Scalar>(g, vertex_radius, Graph::Vertex(alpha0(g, eg.dart)));

		Scalar avg_radius = (R0 + Rn) / 2;
		Vec3 edge = Pn - P0;
		Scalar D = edge.norm();
		edge = edge.normalized();

		uint32 n = uint32(D / avg_radius);
		// if (D > 2 * avg_radius)
		//{
		//	n = 2;
		if (n > 1)
		{
			// Scalar y = (Rn - R0) / D;
			Scalar ratio = pow(Rn / R0, 1.0 / Scalar(n));

			Scalar sum_ratio = 0;
			Scalar alpha = 1;
			std::vector<Scalar> ratios_sums;
			for (uint32 i = 0; i < n; ++i)
			{
				sum_ratio += alpha;
				ratios_sums.push_back(sum_ratio);
				alpha *= ratio;
			}
			ratios_sums.pop_back();

			Vec3 D0 = D / sum_ratio * edge;
			std::vector<Vec3> Pi;
			for (Scalar r : ratios_sums)
			{
				Pi.push_back(P0 + r * D0);
			}

			Dart d = eg.dart;
			alpha = ratio * R0;
			for (Vec3 P : Pi)
			{
				d = cut_edge(g, Graph::Edge(d), true).dart;
				value<Vec3>(g, vertex_position, Graph::Vertex(d)) = P;
				value<Scalar>(g, vertex_radius, Graph::Vertex(d)) = alpha;
				alpha *= ratio;
				d = alpha1(g, d);
			}
		}
		//}
		return true;
	});
	return true;
}

bool add_graph_attributes(Graph& g, GAttributes& gAttribs)
{
	gAttribs.vertex_position = get_attribute<Vec3, Graph::Vertex>(g, "position");
	if (!gAttribs.vertex_position)
	{
		std::cout << "The graph has no vertex position attribute" << std::endl;
		return false;
	}
	gAttribs.vertex_radius = get_attribute<Scalar, Graph::Vertex>(g, "radius");
	if (!gAttribs.vertex_radius)
	{
		std::cout << "The graph has no vertex radius attribute" << std::endl;
		return false;
	}
	gAttribs.vertex_contact_surface = add_attribute<Dart, Graph::Vertex>(g, "contact_surface");
	gAttribs.halfedge_contact_surface_face = add_attribute<Dart, Graph::HalfEdge>(g, "contact_surface_face");
	gAttribs.halfedge_frame = add_attribute<Mat3, Graph::HalfEdge>(g, "frame");
	gAttribs.halfedge_volume_connection = add_attribute<Dart, Graph::HalfEdge>(g, "volume_connection");
	return true;
}

bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs)
{
	m2Attribs.vertex_position = add_attribute<Vec3, CMap2::Vertex>(m2, "position");
	m2Attribs.dual_vertex_graph_branch = add_attribute<Dart, CMap2::Vertex>(m2, "graph_branch");
	m2Attribs.volume_center = add_attribute<Vec3, CMap2::Volume>(m2, "center");
	m2Attribs.volume_gvertex = add_attribute<Graph::Vertex, CMap2::Volume>(m2, "gvertex");
	m2Attribs.edge_mid = add_attribute<Vec3, CMap2::Edge>(m2, "edge_mid");
	m2Attribs.halfedge_volume_connection = add_attribute<Dart, CMap2::HalfEdge>(m2, "volume_connection");
	m2Attribs.ortho_scaffold = add_attribute<CMap2*, CMap2::Volume>(m2, "ortho_scaffold");
	return true;
}

bool add_cmap3_attributes(CMap3& m3, M3Attributes& m3Attribs)
{
	m3Attribs.vertex_position = cgogn::add_attribute<Vec3, CMap3::Vertex>(m3, "position");
	m3Attribs.volume_graph_connection = cgogn::add_attribute<Graph::HalfEdge, CMap3::Volume>(m3, "graph_connection");

	return true;
}

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	bool res = true;
	gAttribs.vertex_contact_surface->fill(Dart());

	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		uint32 d = degree(g, v);
		if (d == 1)
		{
			build_contact_surface_1(g, gAttribs, m2, m2Attribs, v);
			return res;
		}
		if (d == 2)
		{
			build_contact_surface_2(g, gAttribs, m2, m2Attribs, v);
			return res;
		}
		// if (d >= 3 && d <= 6)
		// {
		// 	// check if branches directions are cube-friendly
		// 	// yes -> add a 8-hex intersection block
		// 	bool ortho = build_contact_surface_ortho(g, gAttribs, m2, m2Attribs, v);
		// 	if (ortho)
		// 		return res;
		// }
		// other cases
		// (calls build_contact_surface_orange on planar configurations)
		build_contact_surface_n(g, gAttribs, m2, m2Attribs, v);
		return res;
	});

	return res;
}

void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	Dart d = add_face(m2, 4, true).dart;

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d;
	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(d)) = v;

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(v.dart)) = d;
}

void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, d1);
	phi2_sew(m2, phi1(m2, d0), phi_1(m2, d1));
	phi2_sew(m2, phi<1, 1>(m2, d0), phi<1, 1>(m2, d1));
	phi2_sew(m2, phi_1(m2, d0), phi1(m2, d1));

	index_volume_cells(m2, CMap2::Volume(d0));

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d0;
	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(d0)) = v;

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(v.dart)) = d0;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(alpha1(g, v.dart))) = phi<1, 2>(m2, d0);
}

void build_contact_surface_orange(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								  Graph::Vertex v)
{
	uint32 nbf = degree(g, v);

	// create the topology of the surface
	std::vector<Dart> faces;
	faces.reserve(nbf);
	for (uint32 i = 0; i < nbf; ++i)
		faces.push_back(add_face(static_cast<CMap1&>(m2), 4, false).dart);
	for (uint32 i = 0; i < nbf; ++i)
	{
		Dart d = faces[i];
		Dart e = faces[(i + 1) % nbf];
		phi2_sew(m2, d, phi1(m2, e));
		phi2_sew(m2, phi_1(m2, d), phi<1, 1>(m2, e));
	}

	index_volume_cells(m2, CMap2::Volume(faces[0]));

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = faces[0];

	// get the points on the sphere for each incident branch
	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);

	std::vector<Vec3> Ppos;
	Ppos.reserve(nbf);
	std::vector<Dart> Pdart;
	Pdart.reserve(nbf);
	foreach_dart_of_orbit(g, v, [&](Dart d) {
		Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		geometry::project_on_sphere(p, center, 1);
		Ppos.push_back(p);
		Pdart.push_back(d);
		return true;
	});

	// get the best fitting plane normal and build a local frame based on this normal and first branch
	std::pair<Vec3, Scalar> plane = geometry::plane_fitting(Ppos);
	Vec3 L1 = (Ppos[0] - center).normalized();
	geometry::project_on_plane(L1, plane.first, 0);
	Vec3 L3 = plane.first;
	Vec3 L2 = L3.cross(L1);
	Mat3 L;
	L << L1[0], L1[1], L1[2], L2[0], L2[1], L2[2], L3[0], L3[1], L3[2];

	// sort the incident branches on the plane in CW order
	std::vector<uint32> permutation(nbf);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&](uint32 i, uint32 j) -> bool {
		Vec3 proji = (Ppos[i] - center).normalized();
		geometry::project_on_plane(proji, plane.first, 0);
		proji = L * proji;
		Scalar anglei = atan2(proji[1], proji[0]);

		Vec3 projj = (Ppos[j] - center).normalized();
		geometry::project_on_plane(projj, plane.first, 0);
		projj = L * projj;
		Scalar anglej = atan2(projj[1], projj[0]);

		if (anglei >= 0)
		{
			if (anglej >= 0)
				return anglei < anglej;
			else
				return true;
		}
		else
		{
			if (anglej >= 0)
				return false;
			else
				return anglei < anglej;
		}
	});

	// apply the permutation to branches point & dart
	std::vector<Vec3> sorted_Ppos(nbf);
	std::vector<Dart> sorted_Pdart(nbf);
	std::transform(permutation.begin(), permutation.end(), sorted_Ppos.begin(), [&](uint32 i) { return Ppos[i]; });
	std::transform(permutation.begin(), permutation.end(), sorted_Pdart.begin(), [&](uint32 i) { return Pdart[i]; });

	// put the geometry on the surface mesh vertices
	Vec3 Q1 = center + plane.first;
	Vec3 Q2 = center - plane.first;
	for (uint32 i = 0; i < nbf; ++i)
	{
		value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(sorted_Pdart[i])) = faces[i];
		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(faces[i])) =
			center + (sorted_Ppos[(i + 1) % nbf] - sorted_Ppos[i]).normalized().cross(plane.first);
	}
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, faces[0]))) = Q1;
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, faces[0]))) = Q2;

	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(faces[0])) = v;
}

void build_contact_surface_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	uint32 nbf = degree(g, v);

	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);

	std::vector<Vec3> Ppos;
	Ppos.reserve(nbf);

	foreach_dart_of_orbit(g, v, [&](Dart d) -> bool {
		Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		geometry::project_on_sphere(p, center, 1);
		Ppos.push_back(p);
		return true;
	});
	std::pair<Vec3, Scalar> plane = geometry::plane_fitting(Ppos);
	bool planar = true;
	for (const Vec3& p : Ppos)
	{
		Scalar dist = geometry::distance_plane_point(plane.first, plane.second, p);
		if (dist > 0.15)
		{
			planar = false;
			break;
		}
	}
	if (planar)
	{
		build_contact_surface_orange(g, gAttribs, m2, m2Attribs, v);
		return;
	}

	// compute the n points on the sphere
	// generate Delaunay mesh from the n points
	// store the graph branch on their respective delaunay vertex (m2Attribs.dual_vertex_graph_branch)
	// modify connectivity until all vertices are valence 4
	// call dualize_volume

	Dart vol_dart = convex_hull_around_vertex(g, v, m2, m2Attribs, Ppos);

	index_volume_cells(m2, CMap2::Volume(vol_dart));
	value<Graph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(vol_dart)) = v;

	vol_dart = remesh(m2, CMap2::Volume(vol_dart), m2Attribs);
	dualize_volume(m2, CMap2::Volume(vol_dart), m2Attribs, g, gAttribs);

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = vol_dart;

	return;
}

bool build_contact_surface_ortho(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 Graph::Vertex gv)
{
	Mat3 frame = Mat3();
	bool found_frame = find_inter_frame(g, gv, gAttribs, frame);
	if (!found_frame)
		return false;

	std::cout << "ortho intersection " << degree(g, gv) << std::endl;

	Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, gv);

	Vec3& p = value<Vec3>(g, gAttribs.vertex_position, gv);

	std::vector<Vec3> corners = {(frame.col(0) - frame.col(1) - frame.col(2)).normalized() * radius * 1.15 + p,
								 (frame.col(0) + frame.col(1) - frame.col(2)).normalized() * radius * 1.15 + p,
								 (-frame.col(0) + frame.col(1) - frame.col(2)).normalized() * radius * 1.15 + p,
								 (-frame.col(0) - frame.col(1) - frame.col(2)).normalized() * radius * 1.15 + p,
								 (frame.col(0) - frame.col(1) + frame.col(2)).normalized() * radius * 1.15 + p,
								 (frame.col(0) + frame.col(1) + frame.col(2)).normalized() * radius * 1.15 + p,
								 (-frame.col(0) + frame.col(1) + frame.col(2)).normalized() * radius * 1.15 + p,
								 (-frame.col(0) - frame.col(1) + frame.col(2)).normalized() * radius * 1.15 + p};

	/// create support
	CMap2* scaffold = new CMap2();
	CMap2::Volume w = add_prism(*scaffold, 4);
	auto scaffold_position = add_attribute<Vec3, CMap2::Vertex>(*scaffold, "position");

	Dart d0 = w.dart;
	std::vector<Dart> vertices = {
		d0,									  // 0
		phi1(*scaffold, d0),				  // 1
		phi<1, 2, 1, 1>(*scaffold, d0),		  // 2
		phi<2, 1, 2>(*scaffold, d0),		  // 3
		phi_1(*scaffold, d0),				  // 4
		phi<1, 1>(*scaffold, d0),			  // 5
		phi<1, 1, 2, 1, 1>(*scaffold, d0),	  // 6
		phi<1, 1, 2, 1, 1, 1>(*scaffold, d0), // 7
	};

	for (uint32 i = 0; i < 8; ++i)
		value<Vec3>(*scaffold, scaffold_position, CMap2::Vertex(vertices[i])) = corners[i];

	/// setup geometry for m3 use
	auto scaffold_position_edge = add_attribute<Vec3, CMap2::Edge>(*scaffold, "position");
	foreach_cell(*scaffold, [&](CMap2::Edge e) -> bool {
		Vec3 center = {0, 0, 0};
		foreach_incident_vertex(*scaffold, e, [&](CMap2::Vertex v) -> bool {
			center += value<Vec3>(*scaffold, scaffold_position, v);
			return true;
		});
		center /= 2;
		value<Vec3>(*scaffold, scaffold_position_edge, e) = (center - p).normalized() * radius * 1.05 + p;
		return true;
	});
	auto scaffold_position_face = add_attribute<Vec3, CMap2::Face>(*scaffold, "position");
	foreach_cell(*scaffold, [&](CMap2::Face f) -> bool {
		Vec3 center = {0, 0, 0};
		foreach_incident_vertex(*scaffold, f, [&](CMap2::Vertex v) -> bool {
			center += value<Vec3>(*scaffold, scaffold_position, v);
			return true;
		});
		center /= codegree(*scaffold, f);
		value<Vec3>(*scaffold, scaffold_position_face, f) = (center - p).normalized() * radius + p;
		return true;
	});
	auto scaffold_position_volume = add_attribute<Vec3, CMap2::Volume>(*scaffold, "position");

	value<Vec3>(*scaffold, scaffold_position_volume, w) = p;

	// associate branches with cube faces
	auto scaffold_face_branch = add_attribute<CMap2::HalfEdge, CMap2::Face>(*scaffold, "face_branch");
	std::vector<CMap2::Face> faces = {
		CMap2::Face(phi2(*scaffold, vertices[0])),		   CMap2::Face(vertices[0]),
		CMap2::Face(phi<2, 1, 2>(*scaffold, vertices[6])), CMap2::Face(phi2(*scaffold, vertices[6])),
		CMap2::Face(phi<2, 1, 2>(*scaffold, vertices[0])), CMap2::Face(vertices[6])};

	CellCache<CMap2> cache_active_faces(*scaffold);
	Vec3 A = value<Vec3>(g, gAttribs.vertex_position, gv);
	std::vector<uint32> face_axis = {1, 3, 2, 4, 5, 0};
	std::vector<uint32> face_axis_used = {false, false, false, false, false, false};
	bool same_face_twice = false;
	// std::cout << "active faces: " << std::endl;
	foreach_dart_of_orbit(g, gv, [&](Dart d) -> bool {
		Vec3 B = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		Vec3 AB = (B - A).normalized();
		uint32 faceid = 0;
		Scalar maxdot = 0;
		for (uint32 i = 0; i < 3; ++i)
		{
			Scalar dot = AB.dot(frame.col(i));
			if (abs(dot) > maxdot)
			{
				if (dot > 0)
					faceid = 2 * i;
				else
					faceid = 2 * i + 1;
				maxdot = abs(dot);
			}
			// Scalar dot = round(AB.dot(frame.col(i)));
		}
		// std::cout << "maxdot = " << maxdot << std::endl;
		// std::cout << "faceid = " << faceid << std::endl;
		if (maxdot > 0.4)
		{
			// std::cout << faceid << ",";
			value<Graph::HalfEdge>(*scaffold, scaffold_face_branch, faces[face_axis[faceid]]) = Graph::HalfEdge(d);
			cache_active_faces.add(faces[face_axis[faceid]]);
			if (face_axis_used[faceid])
				same_face_twice = true;
			face_axis_used[faceid] = true;
		}

		// if (dot == 1)
		// {
		// 	value<Graph::HalfEdge>(*scaffold, scaffold_face_branch, faces[face_axis[2 * i]]) = Graph::HalfEdge(d);
		// 	cache_active_faces.add(faces[face_axis[2 * i]]);
		// 	if (face_axis_used[2 * i])
		// 		same_face_twice = true;
		// 	face_axis_used[2 * i] = true;
		// 	std::cout << 2 * i << ",";
		// }
		// if (dot == -1)
		// {
		// 	value<Graph::HalfEdge>(*scaffold, scaffold_face_branch, faces[face_axis[2 * i + 1]]) = Graph::HalfEdge(d);
		// 	cache_active_faces.add(faces[face_axis[2 * i + 1]]);
		// 	if (face_axis_used[2 * i + 1])
		// 		same_face_twice = true;
		// 	face_axis_used[2 * i + 1] = true;
		// 	std::cout << 2 * i + 1 << ",";
		// }
		return true;
	});
	std::cout << std::endl;
	if (same_face_twice)
	{
		delete scaffold;
		std::cout << "abort.." << std::endl;
		return false;
	}

	/// start building connection surfaces
	auto scaffold_cs_connection = add_attribute<Dart, CMap2::HalfEdge>(*scaffold, "scaffold_connection");

	CellMarker<CMap2, CMap2::Vertex> active_vertices(*scaffold);
	foreach_cell(cache_active_faces, [&](CMap2::Face f) -> bool {
		Dart dcs = add_face(static_cast<CMap1&>(m2), 4, false).dart;
		Dart df = f.dart;
		do
		{
			value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(df)) = dcs;
			active_vertices.mark(CMap2::Vertex(df));
			dcs = phi1(m2, dcs);
			df = phi1(*scaffold, df);
		} while (df != f.dart);

		return true;
	});

	Dart d_ind; // capture of a dart of the volume being built
	foreach_cell(*scaffold, [&](CMap2::Vertex v) -> bool {
		if (active_vertices.is_marked(v))
		{
			Dart d0;
			foreach_dart_of_orbit(*scaffold, v, [&](Dart d2) -> bool {
				if (!value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(d2)).is_nil())
					d0 = d2;
				return d0.is_nil();
			});

			// if (d_ind.is_nil())
			// 	d_ind = d0;

			std::vector<Dart> path;
			Dart d = d0;
			do
			{
				active_vertices.unmark(CMap2::Vertex(d));
				path.push_back(value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(d)));
				if (value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(phi2(*scaffold, d))).is_nil())
					d = phi1(*scaffold, d);
				else
					d = phi<2, 1>(*scaffold, d);
			} while (d != d0);

			Dart d_new = add_face(static_cast<CMap1&>(m2), uint32(path.size()), false).dart;
			for (uint32 i = 0; i < path.size(); ++i)
			{
				phi2_sew(m2, path[i], d_new);
				d_new = phi_1(m2, d_new);
			}

			d_ind = d_new;
		}
		return true;
	});

	DartMarkerStore<CMap2> vertex_marker(m2);
	foreach_dart_of_orbit(m2, CMap2::Volume(d_ind), [&](Dart d) -> bool {
		if (!vertex_marker.is_marked(d))
		{
			CMap2::Vertex v(d);
			foreach_dart_of_orbit(m2, v, [&](Dart d) -> bool {
				vertex_marker.mark(d);
				return true;
			});
		}
		return true;
	});

	index_volume_cells(m2, CMap2::Volume(d_ind));

	value<CMap2*>(m2, m2Attribs.ortho_scaffold, CMap2::Volume(d_ind)) = scaffold;
	value<Dart>(g, gAttribs.vertex_contact_surface, gv) = d_ind;

	/// add connection graph halfedge -> contact surface face
	foreach_cell(*scaffold, [&](CMap2::Face fs) -> bool {
		Graph::HalfEdge g_he = value<Graph::HalfEdge>(*scaffold, scaffold_face_branch, fs);
		if (g_he.is_valid())
		{
			value<Dart>(g, gAttribs.halfedge_contact_surface_face, g_he) =
				value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(fs.dart));
		}
		return true;
	});

	// // debug
	// auto cs_position = add_attribute<Vec3, CMap2::Vertex>(m2, "position");
	foreach_cell(*scaffold, [&](CMap2::Edge e) -> bool {
		Dart d0 = value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(e.dart));
		Dart d1 = value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(phi2(*scaffold, e.dart)));
		if (!d0.is_nil())
		{
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, d0))) =
				value<Vec3>(*scaffold, scaffold_position_edge, e);
		}
		else if (!d1.is_nil())
		{
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, d1))) =
				value<Vec3>(*scaffold, scaffold_position_edge, e);
		}

		return true;
	});

	return true;
}

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	bool res = true;
	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		if (degree(g, v) > 2)
			res = create_intersection_frame_n(g, gAttribs, m2, m2Attribs, v);
		return res;
	});
	return res;
}

bool create_intersection_frame_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 Graph::Vertex v)
{
	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	foreach_dart_of_orbit(g, v, [&](Dart d) -> bool {
		Dart d0 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(d));
		Dart d1 = phi<1, 1>(m2, d0);

		Vec3 R, S, T, diag, temp;
		T = (value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d))) - center).normalized();
		diag = (value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d1)) -
				value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d0)))
				   .normalized();
		R = diag.cross(T).normalized();
		S = T.cross(R).normalized();

		Mat3& f = value<Mat3>(g, gAttribs.halfedge_frame, Graph::HalfEdge(d));
		f.col(0) = R;
		f.col(1) = S;
		f.col(2) = T;

		return true;
	});
	return true;
}

bool create_extremity_frame(const Graph& g, GAttributes& gAttribs, Graph::Vertex v)
{
	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	Vec3 R, S, T, temp;
	T = (value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, v.dart))) - center).normalized();
	temp = Vec3(T[1], -T[0], T[2]);

	R = temp.cross(T).normalized();
	S = T.cross(R).normalized();

	Mat3& f = value<Mat3>(g, gAttribs.halfedge_frame, Graph::HalfEdge(v.dart));
	f.col(0) = R;
	f.col(1) = S;
	f.col(2) = T;

	return true;
}

bool propagate_frames(const Graph& g, GAttributes& gAttribs, const GraphData& gData, CMap2& m2)
{
	for (auto& branch : gData.branches)
	{
		if (degree(g, Graph::Vertex(branch.first.dart)) > 1)
		{
			if (degree(g, Graph::Vertex(branch.second.dart)) == 1)
				propagate_frame_n_1(g, gAttribs, branch.first);
			else
				propagate_frame_n_n(g, gAttribs, m2, branch.first);
		}
		else
		{
			if (degree(g, Graph::Vertex(branch.second.dart)) == 1)
				create_extremity_frame(g, gAttribs, Graph::Vertex(branch.second.dart));

			propagate_frame_n_1(g, gAttribs, branch.second);
		}
	}
	return true;
}

void propagate_frame_n_1(const Graph& g, GAttributes& gAttribs, Graph::HalfEdge h_from_start)
{
	Dart d_from = h_from_start.dart;
	Dart d_to = alpha0(g, d_from);
	Graph::HalfEdge h_from(d_from);
	Graph::HalfEdge h_to(d_to);
	Graph::Vertex v_from(d_from);
	Graph::Vertex v_to(d_to);

	uint32 valence = degree(g, v_to);

	Mat3 U = value<Mat3>(g, gAttribs.halfedge_frame, h_from);
	Mat3 U_;

	while (valence == 2)
	{
		Dart next_d_from = alpha1(g, d_to);
		Dart next_d_to = alpha0(g, next_d_from);
		Graph::HalfEdge next_h_from(next_d_from);
		Graph::HalfEdge next_h_to(next_d_to);
		Graph::Vertex next_v_from(next_d_from);
		Graph::Vertex next_v_to(next_d_to);

		Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
		Scalar c1, c2;

		Ri = U.col(0);
		Ti = U.col(2);
		v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
		c1 = v1.dot(v1);
		RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
		TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;

		v3 = value<Vec3>(g, gAttribs.vertex_position, next_v_to) - value<Vec3>(g, gAttribs.vertex_position, v_to);
		Ti1 = (v1.normalized() + v3.normalized()).normalized();
		v2 = Ti1 - TiL;
		c2 = v2.dot(v2);
		Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
		Si1 = Ti1.cross(Ri1);

		U.col(0) = Ri1;
		U.col(1) = Si1;
		U.col(2) = Ti1;

		U_.col(0) = U.col(0);
		U_.col(1) = -U.col(1);
		U_.col(2) = -U.col(2);

		value<Mat3>(g, gAttribs.halfedge_frame, h_to) = U_;
		value<Mat3>(g, gAttribs.halfedge_frame, next_h_from) = U;

		d_from = next_d_from;
		d_to = next_d_to;
		v_from = next_v_from;
		v_to = next_v_to;
		h_from = next_h_from;
		h_to = next_h_to;

		valence = degree(g, v_to);
	}

	if (valence == 1)
	{
		Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
		Scalar c1, c2;

		Ri = U.col(0);
		Ti = U.col(2);
		v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
		c1 = v1.dot(v1);
		RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
		TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;
		Ti1 = v1.normalized();
		v2 = Ti1 - TiL;
		c2 = v2.dot(v2);
		Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
		Si1 = Ti1.cross(Ri1);

		U.col(0) = Ri1;
		U.col(1) = Si1;
		U.col(2) = Ti1;

		U_.col(0) = U.col(0);
		U_.col(1) = -U.col(1);
		U_.col(2) = -U.col(2);

		value<Mat3>(g, gAttribs.halfedge_frame, h_to) = U_;
	}
}

bool propagate_frame_n_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::HalfEdge h_from_start)
{
	Dart d_from = h_from_start.dart;
	Dart d_to = alpha0(g, d_from);
	Graph::HalfEdge h_from(d_from);
	Graph::HalfEdge h_to(d_to);
	Graph::Vertex v_from(d_from);
	Graph::Vertex v_to(d_to);

	uint32 valence = degree(g, v_to);

	Mat3 U = value<Mat3>(g, gAttribs.halfedge_frame, h_from);
	Mat3 U_;

	uint32 nb_e = 0;
	Vec3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
	Scalar c1, c2;

	while (valence == 2)
	{
		Dart next_d_from = alpha1(g, d_to);
		Dart next_d_to = alpha0(g, next_d_from);
		Graph::HalfEdge next_h_from(next_d_from);
		Graph::HalfEdge next_h_to(next_d_to);
		Graph::Vertex next_v_from(next_d_from);
		Graph::Vertex next_v_to(next_d_to);

		nb_e++;

		Ri = U.col(0);
		Ti = U.col(2);
		v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
		c1 = v1.dot(v1);
		RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
		TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;

		v3 = value<Vec3>(g, gAttribs.vertex_position, next_v_to) - value<Vec3>(g, gAttribs.vertex_position, v_to);
		Ti1 = (v1.normalized() + v3.normalized()).normalized();
		v2 = Ti1 - TiL;
		c2 = v2.dot(v2);
		Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
		Si1 = Ti1.cross(Ri1);

		U.col(0) = Ri1;
		U.col(1) = Si1;
		U.col(2) = Ti1;

		value<Mat3>(g, gAttribs.halfedge_frame, next_h_from) = U;

		d_from = next_d_from;
		d_to = next_d_to;
		v_from = next_v_from;
		v_to = next_v_to;
		h_from = next_h_from;
		h_to = next_h_to;

		valence = degree(g, v_to);
	}

	Ri = U.col(0);
	Ti = U.col(2);
	v1 = value<Vec3>(g, gAttribs.vertex_position, v_to) - value<Vec3>(g, gAttribs.vertex_position, v_from);
	c1 = v1.dot(v1);
	RiL = Ri - (2 / c1) * (v1.dot(Ri)) * v1;
	TiL = Ti - (2 / c1) * (v1.dot(Ti)) * v1;
	Ti1 = v1.normalized();
	v2 = Ti1 - TiL;
	c2 = v2.dot(v2);
	Ri1 = RiL - (2 / c2) * (v2.dot(RiL)) * v2;
	Si1 = Ti1.cross(Ri1);

	U.col(0) = Ri1;
	U.col(1) = Si1;
	U.col(2) = Ti1;

	U_.col(0) = U.col(0);
	U_.col(1) = -U.col(1);
	U_.col(2) = -U.col(2);

	Vec3 X = (U_.col(0) + U_.col(1)).normalized();
	Mat3 UE = value<Mat3>(g, gAttribs.halfedge_frame, h_to);
	Vec3 RE = UE.col(0), SE = UE.col(1);
	bool A = (RE.dot(X) >= 0);
	bool B = (SE.dot(X) >= 0);

	uint32 nb_shifts = 0;
	if (!A && B)
		nb_shifts = 1;
	else if (!A && !B)
		nb_shifts = 2;
	else if (A && !B)
		nb_shifts = 3;
	if (nb_shifts > 0)
	{
		shift_frame(UE, nb_shifts);
		value<Mat3>(g, gAttribs.halfedge_frame, h_to) = UE;
		Dart csface = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h_to);
		for (uint32 i = 0; i < nb_shifts; ++i)
			csface = phi1(m2, csface);
		value<Dart>(g, gAttribs.halfedge_contact_surface_face, h_to) = csface;
	}

	if (nb_e > 0)
	{
		Scalar cos0 = UE.col(0).dot(U_.col(0));
		Scalar cos1 = UE.col(1).dot(U_.col(0));
		Scalar angle = std::acos(cos0) * (cos1 > 0 ? 1 : -1);
		Scalar angle_step = -angle / Scalar(nb_e);

		d_from = h_from_start.dart;
		d_to = alpha0(g, d_from);
		h_from.dart = d_from;
		h_to.dart = d_to;
		v_from.dart = d_from;
		v_to.dart = d_to;

		valence = degree(g, v_to);
		uint32 step = 0;
		while (valence == 2)
		{
			step++;

			Dart next_d_from = alpha1(g, d_to);
			Dart next_d_to = alpha0(g, next_d_from);
			Graph::HalfEdge next_h_from(next_d_from);
			Graph::HalfEdge next_h_to(next_d_to);
			Graph::Vertex next_v_from(next_d_from);
			Graph::Vertex next_v_to(next_d_to);

			U = value<Mat3>(g, gAttribs.halfedge_frame, next_h_from);
			Eigen::AngleAxisd rot(-angle_step * step, U.col(2));
			U.col(0) = rot * U.col(0);
			U.col(1) = U.col(2).cross(U.col(0));

			U_.col(0) = U.col(0);
			U_.col(1) = -U.col(1);
			U_.col(2) = -U.col(2);

			value<Mat3>(g, gAttribs.halfedge_frame, h_to) = U_;
			value<Mat3>(g, gAttribs.halfedge_frame, next_h_from) = U;

			d_from = next_d_from;
			d_to = next_d_to;
			v_from = next_v_from;
			v_to = next_v_to;
			h_from = next_h_from;
			h_to = next_h_to;

			valence = degree(g, v_to);
		}
	}

	return true;
}

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const Graph& g, const GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		CMap2::Volume contact_surface(value<Dart>(g, gAttribs.vertex_contact_surface, v));

		const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
		// if (degree(g, v) == 1)
		// 	value<Vec3>(m2, m2Attribs.volume_center, contact_surface) =
		// 		center + (0.25 * (center - value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, v.dart)))));
		// else
		value<Vec3>(m2, m2Attribs.volume_center, contact_surface) = center;

		Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);

		if (degree(g, v) < 3)
		{
			Graph::HalfEdge h(v.dart);
			Dart csf = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h);
			Mat3 frame = value<Mat3>(g, gAttribs.halfedge_frame, h);

			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(csf)) = center - frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, csf))) = center + frame.col(0) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi<1, 1>(m2, csf))) =
				center + frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, csf))) = center - frame.col(0) * radius;
		}
		else
		{
			foreach_incident_vertex(m2, contact_surface, [&](CMap2::Vertex v2) -> bool {
				Vec3 pos = value<Vec3>(m2, m2Attribs.vertex_position, v2);
				geometry::project_on_sphere(pos, center, radius);
				value<Vec3>(m2, m2Attribs.vertex_position, v2) = pos;
				return true;
			});
		}

		foreach_incident_edge(m2, contact_surface, [&](CMap2::Edge e) -> bool {
			std::vector<CMap2::Vertex> vertices = incident_vertices(m2, e);
			Vec3 mid = 0.5 * (value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) +
							  value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]));
			geometry::project_on_sphere(mid, center, radius);
			value<Vec3>(m2, m2Attribs.edge_mid, e) = mid;
			return true;
		});
		return true;
	});

	return true;
}

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/

void insert_ortho_chunks(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		CMap2::Volume contact_surface(value<Dart>(g, gAttribs.vertex_contact_surface, v));
		CMap2* scaffold = value<CMap2*>(m2, m2Attribs.ortho_scaffold, contact_surface);
		if (scaffold)
		{
			auto scaffold_cs_connection = get_attribute<Dart, CMap2::HalfEdge>(*scaffold, "scaffold_connection");
			auto scaffold_hex_connection = add_attribute<Dart, CMap2::HalfEdge>(*scaffold, "hex_connection");
			foreach_cell(*scaffold, [&](CMap2::Vertex v) -> bool {
				Dart dv = v.dart;
				Dart d3 = add_prism(static_cast<CMap2&>(m3), 4).dart;
				do
				{
					value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(dv)) = d3;

					Dart cs_dart = value<Dart>(*scaffold, scaffold_cs_connection, CMap2::HalfEdge(dv));
					if (!cs_dart.is_nil())
						value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(phi2(m2, cs_dart))) =
							phi_1(m3, d3);

					dv = phi<2, 1>(*scaffold, dv);
					d3 = phi<2, 1>(m3, d3);
				} while (dv != v.dart);

				return true;
			});

			foreach_cell(*scaffold, [&](CMap2::Edge e) -> bool {
				Dart dc0 = value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(e.dart));
				Dart dc1 = value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(phi2(*scaffold, e.dart)));
				sew_volumes(m3, phi<1, 2>(m3, dc0), phi<1, 2, 1>(m3, dc1));
				return true;
			});
		}
		return true;
	});
}

bool build_branch_sections(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	insert_ortho_chunks(g, gAttribs, m2, m2Attribs, m3);

	foreach_cell(g, [&](Graph::Edge e) -> bool {
		std::vector<Graph::HalfEdge> halfedges = incident_halfedges(g, e);

		Dart m2f0 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, halfedges[0]);
		Dart m2f1 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, halfedges[1]);

		std::vector<Dart> F0 = {m2f0, phi1(m2, m2f0), phi<1, 1>(m2, m2f0), phi_1(m2, m2f0)};
		std::vector<Dart> F1 = {m2f1, phi1(m2, m2f1), phi<1, 1>(m2, m2f1), phi_1(m2, m2f1)};

		Dart m3d = add_branch_section(m3);
		std::vector<Dart> D0 = {m3d, phi<2, 3, 2, 1>(m3, m3d), phi<2, 3, 2, 1, 2, 3, 2, 1>(m3, m3d),
								phi<1, 1, 1, 2, 3, 2>(m3, m3d)};
		std::vector<Dart> D1 = {phi<2, 1, 1, 2>(m3, D0[0]), phi<2, 1, 1, 2>(m3, D0[1]), phi<2, 1, 1, 2>(m3, D0[2]),
								phi<2, 1, 1, 2>(m3, D0[3])};

		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);

		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[0])) = phi<1, 1>(m3, D1[1]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[1])) = phi<1, 1>(m3, D1[0]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[2])) = phi<1, 1>(m3, D1[3]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[3])) = phi<1, 1>(m3, D1[2]);

		value<Dart>(g, gAttribs.halfedge_volume_connection, halfedges[0]) = D0[0];
		value<Dart>(g, gAttribs.halfedge_volume_connection, halfedges[1]) = phi1(m3, D1[0]);
		return true;
	});

	return true;
}

bool sew_branch_sections(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		if (is_incident_to_boundary(m2, e))
			return true;

		std::vector<CMap2::HalfEdge> halfedges = incident_halfedges(m2, e);
		sew_volumes(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, halfedges[0]),
					phi1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, halfedges[1])));
		return true;
	});

	close(m3, false);

	return true;
}

bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, M3Attributes& m3Attribs)
{
	parallel_foreach_cell(m2, [&](CMap2::Volume w) -> bool {
		Dart m3d = phi_1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(w.dart)));
		value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.volume_center, w);
		return true;
	});

	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		Dart m3d = phi1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(e.dart)));
		value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.edge_mid, e);
		return true;
	});

	CellMarker<CMap3, CMap3::Vertex> cm(m3);
	for (Dart m2d = m2.begin(), end = m2.end(); m2d != end; m2d = m2.next(m2d))
	{
		if (!is_boundary(m2, m2d))
		{
			Dart m3d = value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(m2d));
			CMap3::Vertex m3v(m3d);
			if (!cm.is_marked(m3v))
			{
				cm.mark(m3v);
				value<Vec3>(m3, m3Attribs.vertex_position, m3v) =
					value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, m2d)));
			}
		}
	}

	parallel_foreach_cell(m2, [&](CMap2::Volume w) -> bool {
		CMap2* scaffold = value<CMap2*>(m2, m2Attribs.ortho_scaffold, w);
		if (scaffold)
		{
			auto scaffold_cs_connection = get_attribute<Dart, CMap2::HalfEdge>(*scaffold, "scaffold_connection");
			auto scaffold_position = get_attribute<Vec3, CMap2::Vertex>(*scaffold, "position");
			auto scaffold_position_edge = get_attribute<Vec3, CMap2::Edge>(*scaffold, "position");
			auto scaffold_position_face = get_attribute<Vec3, CMap2::Face>(*scaffold, "position");
			auto scaffold_position_volume = get_attribute<Vec3, CMap2::Volume>(*scaffold, "position");
			auto scaffold_hex_connection = get_attribute<Dart, CMap2::HalfEdge>(*scaffold, "hex_connection");
			// uint32 i = 0;
			foreach_cell(*scaffold, [&](CMap2::Vertex v2) -> bool {
				Dart d3 = value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(v2.dart));
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(d3)) =
					value<Vec3>(*scaffold, scaffold_position, v2);
				return true;
			});
			foreach_cell(*scaffold, [&](CMap2::Edge e2) -> bool {
				Dart d3 = value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(e2.dart));
				d3 = phi1(m3, d3);
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(d3)) =
					value<Vec3>(*scaffold, scaffold_position_edge, e2);

				return true;
			});
			foreach_cell(*scaffold, [&](CMap2::Face f2) -> bool {
				Dart d3 = value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(f2.dart));
				d3 = phi<1, 1>(m3, d3);
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(d3)) =
					value<Vec3>(*scaffold, scaffold_position_face, f2);
				return true;
			});
			foreach_cell(*scaffold, [&](CMap2::Volume w2) -> bool {
				Dart d3 = value<Dart>(*scaffold, scaffold_hex_connection, CMap2::HalfEdge(w2.dart));
				d3 = phi<1, 1, 2, 1, 1>(m3, d3);
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(d3)) =
					value<Vec3>(*scaffold, scaffold_position_volume, w2);
				return true;
			});
		}
		return true;
	});

	return true;
}

/*****************************************************************************/
/* utils				                                                     */
/*****************************************************************************/

bool dijkstra_topo(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
				   std::shared_ptr<CMap2::Attribute<uint32>> dist)
{
	foreach_incident_vertex(m2, CMap2::Volume(v0.dart), [&](CMap2::Vertex v) -> bool {
		value<Dart>(m2, previous, v) = Dart();
		value<uint32>(m2, dist, v) = UINT_MAX;
		return true;
	});

	DartMarker visited(m2);

	std::vector<CMap2::Vertex> vertices = {v0};
	value<uint32>(m2, dist, v0) = 0;
	CMap2::Vertex v_act;
	uint32 dist_act;
	while (vertices.size())
	{
		v_act = vertices[vertices.size() - 1];
		vertices.pop_back();
		dist_act = value<uint32>(m2, dist, v_act) + 1;

		std::vector<CMap2::Vertex> neighbors;
		foreach_dart_of_orbit(m2, v_act, [&](Dart d) -> bool {
			if (!visited.is_marked(d))
			{
				Dart d2 = phi2(m2, d);
				visited.mark(d);
				visited.mark(d2);

				uint32 dist_2 = value<uint32>(m2, dist, CMap2::Vertex(d2));
				if (dist_act < dist_2)
				{
					value<uint32>(m2, dist, CMap2::Vertex(d2)) = dist_act;
					value<Dart>(m2, previous, CMap2::Vertex(d2)) = d;
					neighbors.push_back(CMap2::Vertex(d2));
				}
			}

			return true;
		});

		vertices.insert(vertices.begin(), neighbors.begin(), neighbors.end());
	}

	return false;
}

Dart convex_hull_around_vertex(const Graph& g, Graph::Vertex v, CMap2& m2, M2Attributes& m2Attribs,
							   std::vector<Vec3>& Ppos)
{
	std::vector<uint32> Pid;
	Pid.reserve(Ppos.size());

	uint32 id = 0;
	foreach_dart_of_orbit(g, v, [&](Dart d) -> bool {
		uint32 vertex_id = new_index<CMap2::Vertex>(m2);
		Pid.push_back(vertex_id);
		(*m2Attribs.vertex_position)[vertex_id] = Ppos[id++]; // positions are in the same order in Ppos
		(*m2Attribs.dual_vertex_graph_branch)[vertex_id] = d;
		return true;
	});

	std::vector<uint32> faces_nb_vertices;
	std::vector<uint32> faces_vertex_indices;

	std::vector<uint32> indices;
	for (uint32 i = 0; i < Ppos.size() - 2; ++i)
	{
		for (uint32 j = i + 1; j < Ppos.size() - 1; ++j)
		{
			for (uint32 k = j + 1; k < Ppos.size(); ++k)
			{
				Vec3 t0 = Ppos[j] - Ppos[i];
				Vec3 t1 = Ppos[k] - Ppos[i];
				Vec3 n = t0.cross(t1);
				int sign = 0;

				for (uint32 m = 0; m < Ppos.size(); ++m)
				{
					if (m == i || m == j || m == k)
						continue;

					Vec3 vec = Ppos[m] - Ppos[i];
					Scalar d = vec.dot(n);

					if (!sign)
						sign = (d < 0 ? -1 : 1);
					else
					{
						if (sign != (d < 0 ? -1 : 1))
						{
							sign = 0;
							break;
						}
					}
				}

				if (sign != 0)
				{
					if (sign == 1)
					{
						indices.push_back(Pid[j]);
						indices.push_back(Pid[i]);
						indices.push_back(Pid[k]);
					}
					else
					{
						indices.push_back(Pid[i]);
						indices.push_back(Pid[j]);
						indices.push_back(Pid[k]);
					}

					faces_nb_vertices.push_back(3);
					faces_vertex_indices.insert(faces_vertex_indices.end(), indices.begin(), indices.end());
					indices.clear();
				}
			}
		}
	}

	auto darts_per_vertex = add_attribute<std::vector<Dart>, CMap2::Vertex>(m2, "__darts_per_vertex");

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(16u);

	Dart volume_dart;
	std::vector<Dart> all_darts;
	for (std::size_t i = 0u, end = faces_nb_vertices.size(); i < end; ++i)
	{
		uint32 nbv = faces_nb_vertices[i];

		vertices_buffer.clear();

		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = faces_vertex_indices[faces_vertex_index++];
			vertices_buffer.push_back(idx);
		}

		if (nbv > 2u)
		{
			CMap1::Face f = add_face(static_cast<CMap1&>(m2), nbv, false);
			Dart d = f.dart;
			for (uint32 j = 0u; j < nbv; ++j)
			{
				all_darts.push_back(d);
				const uint32 vertex_index = vertices_buffer[j];
				set_index<CMap2::Vertex>(m2, d, vertex_index);
				(*darts_per_vertex)[vertex_index].push_back(d);
				d = phi1(m2, d);
			}
		}
	}

	for (Dart d : all_darts)
	{
		if (phi2(m2, d) == d)
		{
			uint32 vertex_index = index_of(m2, CMap2::Vertex(d));

			const std::vector<Dart>& next_vertex_darts =
				value<std::vector<Dart>>(m2, darts_per_vertex, CMap2::Vertex(phi1(m2, d)));
			bool phi2_found = false;

			for (auto it = next_vertex_darts.begin(); it != next_vertex_darts.end() && !phi2_found; ++it)
			{
				if (index_of(m2, CMap2::Vertex(phi1(m2, *it))) == vertex_index)
				{
					if (phi2(m2, *it) == *it)
					{
						phi2_sew(m2, d, *it);
						phi2_found = true;
					}
				}
			}
		}
	}

	remove_attribute<CMap2::Vertex>(m2, darts_per_vertex);
	return all_darts[0];
}

Vec3 slerp(Vec3 A, Vec3 B, Scalar alpha, bool in)
{
	Scalar phi = geometry::angle(A, B) - (in ? 0 : 2 * M_PI);
	Scalar s0 = std::sin(phi * (1 - alpha));
	Scalar s1 = std::sin(phi * alpha);
	Scalar s2 = std::sin(phi);
	Vec3 sl = A * (s0 / s2);
	sl += B * (s1 / s2);
	return sl;
}

Scalar angle_on_sphere(Vec3 A, Vec3 B, Vec3 C)
{
	Vec3 sB = slerp(A, B, 0.01, true);
	Vec3 sC = slerp(A, C, 0.01, true);
	Vec3 AB = sB - A;
	Vec3 AC = sC - A;
	return geometry::angle(AB, AC);
}

Scalar edge_max_angle(CMap2& m2, CMap2::Edge e, M2Attributes& m2Attribs)
{
	Dart ed0 = e.dart;
	Dart ed1 = phi2(m2, ed0);

	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed0)));

	Scalar a0 = angle_on_sphere(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed1)));
	a0 += angle_on_sphere(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed1)));

	Scalar a1 = angle_on_sphere(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed0)));
	a1 += angle_on_sphere(A, B, C);

	return std::max(a0, a1);
}

Scalar min_cut_angle(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, M2Attributes& m2Attribs)
{
	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v0.dart)));
	Scalar a0 = angle_on_sphere(A, B, C);

	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v1.dart)));
	Scalar a1 = angle_on_sphere(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	B = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v1.dart)));

	Scalar a2 = angle_on_sphere(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v0.dart)));
	Scalar a3 = angle_on_sphere(A, B, C);

	return std::min({a0, a1, a2, a3});
}

Vec3 spherical_barycenter(std::vector<Vec3>& points, uint32 iterations)
{
	uint32 nb_pts = uint32(points.size());
	std::vector<Vec3> pts0(nb_pts);
	std::vector<Vec3> pts1(nb_pts);

	for (uint32 i = 0; i < nb_pts; ++i)
		pts0[i] = points[i];

	std::vector<Vec3> readvec = pts0;
	std::vector<Vec3> writevec = pts1;
	std::vector<Vec3> tempvec;

	for (uint32 it = 0; it < iterations; ++it)
	{
		for (uint32 i = 0; i < nb_pts; ++i)
		{
			writevec[i] = slerp(readvec[i], readvec[(i + 1) % nb_pts], Scalar(0.5), true);
		}

		tempvec = writevec;
		writevec = readvec;
		readvec = tempvec;
	}

	Vec3 bary = {0, 0, 0};
	for (uint32 i = 0; i < nb_pts; ++i)
	{
		bary += readvec[i];
	}
	bary /= nb_pts;

	Vec3 normal = {0, 0, 0};
	for (uint32 i = 0; i < nb_pts; ++i)
	{
		normal += points[i].cross(points[(i + 1) % nb_pts]);
	}
	normal.normalize();

	if (normal.dot(bary) < 0)
		bary *= -1;

	return bary;
}

Dart remesh(CMap2& m2, CMap2::Volume vol, M2Attributes& m2Attribs)
{
	Dart vol_dart = vol.dart;

	std::vector<CMap2::Vertex> valence_sup4;
	std::vector<CMap2::Vertex> valence_3;

	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
		uint32 valence = degree(m2, v);
		if (valence == 3)
			valence_3.push_back(v);
		if (valence > 4)
			valence_sup4.push_back(v);

		return true;
	});

	if (valence_sup4.size())
	{
		auto edge_angle_max = add_attribute<Scalar, CMap2::Edge>(m2, "angle_max");

		std::vector<CMap2::Edge> edges_n_n;
		std::vector<CMap2::Edge> edges_n_4;
		std::vector<CMap2::Edge> candidate_edges;

		do
		{
			edges_n_n.clear();
			edges_n_4.clear();
			foreach_incident_edge(m2, vol, [&](CMap2::Edge e) -> bool {
				auto vertices = incident_vertices(m2, e);

				uint32 deg_0 = degree(m2, vertices[0]);
				uint32 deg_1 = degree(m2, vertices[1]);
				uint32 deg_min = deg_0 < deg_1 ? deg_0 : deg_1;

				if (deg_min > 4)
				{
					edges_n_n.push_back(e);
				}
				else if (deg_min == 4 && deg_0 + deg_1 > 8)
				{
					edges_n_4.push_back(e);
				}

				return true;
			});

			candidate_edges = edges_n_n.size() ? edges_n_n : edges_n_4;
			if (!candidate_edges.size())
				break;
			for (CMap2::Edge e : candidate_edges)
			{
				value<Scalar>(m2, edge_angle_max, e) = edge_max_angle(m2, e, m2Attribs);
			}

			std::sort(candidate_edges.begin(), candidate_edges.end(), [&](CMap2::Edge e0, CMap2::Edge e1) -> bool {
				return value<Scalar>(m2, edge_angle_max, e0) < value<Scalar>(m2, edge_angle_max, e1);
			});

			CMap2::Edge prime_edge = candidate_edges[0];
			auto neigh_vertices = incident_vertices(m2, prime_edge);

			vol_dart = phi_1(m2, prime_edge.dart);
			merge_incident_faces(m2, prime_edge, true);

		} while (candidate_edges.size() > 1);
		remove_attribute<CMap2::Edge>(m2, edge_angle_max);
	}

	std::vector<std::pair<std::pair<CMap2::Vertex, CMap2::Vertex>, Scalar>> candidate_vertices_pairs;

	valence_3.clear();
	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
		if (degree(m2, v) == 3)
			valence_3.push_back(v);
		return true;
	});
	while (valence_3.size())
	{
		candidate_vertices_pairs.clear();
		std::vector<CMap2::Vertex> verts_3;
		verts_3.reserve(10);
		foreach_incident_face(m2, vol, [&](CMap2::Face f) -> bool {
			verts_3.clear();
			foreach_incident_vertex(m2, f, [&](CMap2::Vertex v) -> bool {
				if (degree(m2, v) == 3)
					verts_3.push_back(v);

				return true;
			});

			if (verts_3.size() >= 2)
				for (uint32 i = 0; i < verts_3.size() - 1; ++i)
				{
					for (uint32 j = i + 1; j < verts_3.size(); ++j)
					{
						candidate_vertices_pairs.push_back(
							{{verts_3[i], verts_3[j]}, min_cut_angle(m2, verts_3[i], verts_3[j], m2Attribs)});
					}
				}
			return true;
		});

		if (candidate_vertices_pairs.size())
		{
			std::sort(candidate_vertices_pairs.begin(), candidate_vertices_pairs.end(),
					  [&](auto pair0, auto pair1) -> bool { return pair0.second < pair1.second; });

			auto pair = candidate_vertices_pairs[0].first;

			cut_face(m2, pair.first, pair.second);
		}
		else
		{
			uint32 max_path_length = 0;
			std::shared_ptr<CMap2::Attribute<Dart>> max_previous;
			std::shared_ptr<CMap2::Attribute<uint32>> max_dist;
			std::pair<CMap2::Vertex, CMap2::Vertex> max_shortest_path;
			for (uint32 i = 0; i < valence_3.size(); ++i)
			{
				CMap2::Vertex v = valence_3[i];
				std::shared_ptr<CMap2::Attribute<Dart>> previous =
					add_attribute<Dart, CMap2::Vertex>(m2, "previous" + std::to_string(v.dart.index));
				std::shared_ptr<CMap2::Attribute<uint32>> dist =
					add_attribute<uint32, CMap2::Vertex>(m2, "dist" + std::to_string(v.dart.index));

				dijkstra_topo(m2, v, previous, dist);

				uint32 curr_min = UINT32_MAX;
				CMap2::Vertex curr_min_vert;
				for (CMap2::Vertex v2 : valence_3)
				{
					if (index_of(m2, v) == index_of(m2, v2))
						continue;

					uint32 new_min = value<uint32>(m2, dist, CMap2::Vertex(v2));
					if (new_min < curr_min)
					{

						curr_min = new_min; // value<uint32>(m2, dist, CMap2::Vertex(v2)); ????
						curr_min_vert = v2;
					}
				}

				if (curr_min > max_path_length)
				{
					if (max_previous)
						remove_attribute<CMap2::Vertex>(m2, max_previous);
					if (max_dist)
						remove_attribute<CMap2::Vertex>(m2, max_dist);

					max_previous = previous;
					max_dist = dist;
					max_shortest_path = {v, curr_min_vert};
					max_path_length = curr_min;
				}
			}

			// building path
			std::vector<Dart> path;
			path.reserve(max_path_length);
			Dart dp = value<Dart>(m2, max_previous, max_shortest_path.second);
			do
			{
				path.push_back(dp);
				dp = value<Dart>(m2, max_previous, CMap2::Vertex(dp));
			} while (dp.index != INVALID_INDEX);

			if (!(path.size() % 2))
			{
				std::vector<CMap2::Edge> first_edges;
				first_edges.push_back(CMap2::Edge(phi_1(m2, path.front())));
				first_edges.push_back(CMap2::Edge(phi2(m2, phi1(m2, path.front()))));
				first_edges.push_back(CMap2::Edge(phi1(m2, path.back())));
				first_edges.push_back(CMap2::Edge(phi2(m2, phi_1(m2, path.back()))));

				std::vector<std::pair<Scalar, uint32>> angles;
				angles.push_back({edge_max_angle(m2, first_edges[0], m2Attribs), 0});
				angles.push_back({edge_max_angle(m2, first_edges[1], m2Attribs), 1});
				angles.push_back({edge_max_angle(m2, first_edges[2], m2Attribs), 2});
				angles.push_back({edge_max_angle(m2, first_edges[3], m2Attribs), 3});

				std::sort(angles.begin(), angles.end(),
						  [&](std::pair<Scalar, uint32> p0, std::pair<Scalar, uint32> p1) -> bool {
							  return p0.first < p1.first;
						  });

				CMap2::Vertex v0, v1;
				CMap2::Edge e0;
				switch (angles[0].second)
				{
				case 0:
					v0 = CMap2::Vertex(phi1(m2, path.front()));
					v1 = CMap2::Vertex(phi_1(m2, path.front()));
					e0 = CMap2::Edge(v1.dart);
					path.erase(path.begin());
					break;
				case 1:
					v0 = CMap2::Vertex(phi2(m2, path.front()));
					v1 = CMap2::Vertex(phi<2, 1, 1>(m2, path.front()));
					e0 = CMap2::Edge(phi<2, 1>(m2, path.front()));
					path.erase(path.begin());
					break;

				case 2:
					v0 = CMap2::Vertex(path.back());
					v1 = CMap2::Vertex(phi<1, 1>(m2, path.back()));
					e0 = CMap2::Edge(phi1(m2, path.back()));
					path.pop_back();
					break;

				case 3:
					// v0 = CMap2::Vertex(phi_1(m2, phi2(m2, path.back())));
					v0 = CMap2::Vertex(phi<2, -1>(m2, path.back()));
					v1 = CMap2::Vertex(phi<2, 1>(m2, path.back()));
					e0 = CMap2::Edge(v0.dart);
					path.pop_back();
					break;
				default:
					break;
				}

				vol_dart = cut_face(m2, v0, v1, true).dart;
				merge_incident_faces(m2, e0, true);
			}

			for (uint32 i = 0; i < path.size(); ++i)
			{
				if (!(i % 2))
					cut_face(m2, CMap2::Vertex(path[i]), CMap2::Vertex(phi1(m2, path[i])), true);
				else
				{
					vol_dart = phi1(m2, path[i]);
					merge_incident_faces(m2, CMap2::Edge(path[i]), true);
				}
			}
		}

		valence_3.clear();
		foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
			if (degree(m2, v) == 3)
				valence_3.push_back(v);
			return true;
		});
	}

	return vol_dart;
}
/*
void mark_tranversal_faces(CMap3& m3, CMap2& m2, M2Attributes& m2Attribs, CellMarker<CMap3, CMap3::Face>& cm)
{
	foreach_cell(m2, [&](CMap2::Edge e2) -> bool {
		Dart m3d = value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(e2.dart));
		cm.mark(CMap3::Face(m3d));
		return true;
	});
}

void subdivide_length_wise(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
						   GAttributes& gAttribs)
{
	using MESH = CMap3;
	using Vertex = MESH::Vertex;
	using Edge = MESH::Edge;
	using Face = MESH::Face;
	using Volume = MESH::Volume;

	CellMarker<MESH, Vertex> new_vertices(m3);
	CellMarker<MESH, Edge> edge2cut(m3);
	CellMarker<MESH, Face> face2cut(m3);
	CellMarker<MESH, Volume> volume2cut(m3);
	CellCache<MESH> cache_edge2cut(m3);
	CellCache<MESH> cache_face2cut(m3);
	CellCache<MESH> cache_sideface2cut(m3);
	CellCache<MESH> cache_volume2cut(m3);

	foreach_cell(m3, [&](Face f) -> bool {
		if (trans_faces.is_marked(f))
		{
			cache_face2cut.add(f);
			foreach_incident_edge(m3, f, [&](Edge e) -> bool {
				if (!edge2cut.is_marked(e))
				{
					edge2cut.mark(e);
					cache_edge2cut.add(e);
				}
				return true;
			});
		}
		else
		{
			cache_sideface2cut.add(f);
		}
		return true;
	});

	cut_all_edges(cache_edge2cut, [&](Vertex v) {
		std::vector<Vertex> vertices = adjacent_vertices_through_edge(m3, v);
		Vec3 mid = (value<Vec3>(m3, m3Attribs.vertex_position, vertices[0]) +
					value<Vec3>(m3, m3Attribs.vertex_position, vertices[1])) *
				   Scalar(0.5);
		value<Vec3>(m3, m3Attribs.vertex_position, v) = mid;
		new_vertices.mark(v);
	});

	foreach_cell(cache_face2cut, [&](Face f) -> bool {
		Vertex center_vertex = quadrangulate_face(m3, f);
		value<Vec3>(m3, m3Attribs.vertex_position, center_vertex) = {0, 0, 0};
		foreach_adjacent_vertex_through_edge(m3, center_vertex, [&](Vertex v) -> bool {
			value<Vec3>(m3, m3Attribs.vertex_position, center_vertex) += value<Vec3>(m3, m3Attribs.vertex_position, v);
			return true;
		});
		value<Vec3>(m3, m3Attribs.vertex_position, center_vertex) /= 4;

		foreach_incident_face(m3, center_vertex, [&](Face new_face) -> bool {
			trans_faces.mark(new_face);
			return true;
		});

		foreach_incident_volume(m3, center_vertex, [&](Volume w) -> bool {
			if (!volume2cut.is_marked(w))
			{
				volume2cut.mark(w);
				cache_volume2cut.add(w);
			}
			return true;
		});

		return true;
	});

	foreach_cell(cache_sideface2cut, [&](Face f) -> bool {
		std::vector<Vertex> vertices;
		vertices.reserve(2);
		if (!trans_faces.is_marked(f))
		{
			foreach_incident_vertex(m3, f, [&](Vertex v) -> bool {
				if (new_vertices.is_marked(v))
					vertices.push_back(v);
				return true;
			});
			cut_face(m3, vertices[0], vertices[1]);
			vertices.clear();
		}
		return true;
	});

	foreach_cell(cache_volume2cut, [&](Volume w) -> bool {
		Dart d0 = w.dart;
		// Dart d1 = phi<2, 1>(m3, d0);
		// Dart d2 = phi2(m3, phi_1(m3, d0));
		std::vector<Dart> path0;
		std::vector<Dart> path1;
		std::vector<Dart> path2;
		get_loop_path(m3, d0, path0);
		cut_volume(m3, path0);
		Edge e = cut_face(m3, Vertex(phi2(m3, path0[2])), Vertex(phi2(m3, path0[5])));
		Dart d1 = e.dart;
		Dart d2 = phi3(m3, e.dart);
		get_loop_path(m3, d1, path1);
		get_loop_path(m3, d2, path2);
		cut_volume(m3, path1);
		cut_volume(m3, path2);

		return true;
	});
}

void subdivide_width_wise(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
						  GAttributes& gAttribs)
{
	using MESH = CMap3;
	using Vertex = MESH::Vertex;
	using Edge = MESH::Edge;
	using Face = MESH::Face;
	using Volume = MESH::Volume;

	using GRAPH = Graph;
	using GVertex = GRAPH::Vertex;
	using GEdge = GRAPH::Edge;

	CellCache<GRAPH> graph_edges(g);
	graph_edges.template build<GEdge>();

	foreach_cell(graph_edges, [&](GEdge eg) -> bool {
		cut_chunk(m3, m3Attribs, trans_faces, g, gAttribs, eg, 0.5);
		return true;
	});
}

void cut_chunk(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
			   GAttributes& gAttribs, Graph::Edge eg, Scalar slice)
{
	using MESH = CMap3;
	using Vertex = MESH::Vertex;
	using Edge = MESH::Edge;
	using Face = MESH::Face;
	using Face2 = MESH::Face2;
	using Volume = MESH::Volume;

	using GRAPH = Graph;
	using GVertex = GRAPH::Vertex;
	using GEdge = GRAPH::Edge;
	using GHEdge = GRAPH::HalfEdge;

	CellMarker<MESH, Vertex> new_vertices(m3);
	CellMarker<MESH, Face2> visited_face2(m3);
	CellMarker<MESH, Face> visited_face(m3);
	CellMarker<MESH, Edge> visited_edge(m3);
	CellCache<MESH> cache_edge2cut(m3);
	CellCache<MESH> cache_face2cut(m3);
	CellCache<MESH> cache_vol2cut(m3);

	Dart dg0 = eg.dart;
	Dart d0 = value<Dart>(g, gAttribs.halfedge_volume_connection, GHEdge(dg0));

	std::vector<Dart> face2_stack;
	face2_stack.push_back(d0);
	visited_face2.mark(Face2(d0));

	uint32 nb_f = 0;
	uint32 nb_e = 0;
	for (uint32 i = 0; i < uint32(face2_stack.size()); ++i)
	{
		const Dart f2d = face2_stack[i];
		Dart it = f2d;
		cache_vol2cut.add(Volume(phi<2, 1>(m3, f2d)));
		do
		{
			Edge e = Edge(phi<2, 1>(m3, it));
			Face f = Face(phi2(m3, it));
			if (!visited_edge.is_marked(e))
			{
				visited_edge.mark(e);
				cache_edge2cut.add(e);
				++nb_e;
			}
			if (!visited_face.is_marked(f))
			{
				visited_face.mark(f);
				cache_face2cut.add(f);
				++nb_f;
			}

			const Dart adj = phi<2, 3, 2>(m3, it);
			if (!is_boundary(m3, adj) && !visited_face2.is_marked(Face2(adj)))
			{
				face2_stack.push_back(adj);
				visited_face2.mark(Face2(adj));
			}
			it = phi1(m3, it);
		} while (it != f2d);
	}

	foreach_cell(cache_edge2cut, [&](Edge e) -> bool {
		std::vector<Vertex> vertices = incident_vertices(m3, e);
		Vertex v_mid = cut_edge(m3, e);
		value<Vec3>(m3, m3Attribs.vertex_position, v_mid) =
			(1 - slice) * value<Vec3>(m3, m3Attribs.vertex_position, vertices[0]) +
			slice * value<Vec3>(m3, m3Attribs.vertex_position, vertices[1]);

		return true;
	});

	foreach_cell(cache_face2cut, [&](Face f) -> bool {
		cut_face(m3, Vertex(phi_1(m3, f.dart)), Vertex(phi<1, 1>(m3, f.dart)));
		return true;
	});

	foreach_cell(cache_vol2cut, [&](Volume vol) -> bool {
		std::vector<Dart> path;
		get_loop_path(m3, phi1(m3, vol.dart), path);
		Face new_face = cut_volume(m3, path);
		trans_faces.mark(new_face);
		return true;
	});

	Dart dg1 = alpha0(g, dg0);
	Dart d1 = value<Dart>(g, gAttribs.halfedge_volume_connection, GHEdge(dg1));

	Dart d0_0 = phi<2, 1, 1, 3, 2>(m3, d0);
	Dart d1_0 = phi<2, 1, 1, 3, 2>(m3, d1);

	GVertex gv = cut_edge(g, eg);
	value<Vec3>(g, gAttribs.vertex_position, gv) =
		(1 - slice) * value<Vec3>(g, gAttribs.vertex_position, GVertex(dg0)) +
		slice * value<Vec3>(g, gAttribs.vertex_position, GVertex(dg1));
	value<Dart>(g, gAttribs.halfedge_volume_connection, GHEdge(alpha0(g, dg0))) = d0_0;
	value<Dart>(g, gAttribs.halfedge_volume_connection, GHEdge(alpha0(g, dg1))) = d1_0;
}

void trisect_length_wise(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
						 GAttributes& gAttribs)
{
	using MESH = CMap3;
	using Vertex = MESH::Vertex;
	using Edge = MESH::Edge;
	using Face = MESH::Face;
	using Volume = MESH::Volume;

	CellMarker<MESH, Vertex> new_vertices(m3);
	CellMarker<MESH, Volume> volume_cut(m3);
	CellMarker<MESH, Edge> edge2cut(m3);
	CellMarker<MESH, Face> face2cut(m3);
	// CellMarker<MESH, Face> sideface2cut(m3);
	CellMarker<MESH, Volume> vol2cut(m3);
	CellCache<MESH> cache_edge2cut(m3);
	CellCache<MESH> cache_sideface2cut(m3);
	CellCache<MESH> cache_face2cut(m3);
	CellCache<MESH> cache_vol2cut(m3);
	auto centroids = cgogn::add_attribute<Vec3, Face>(m3, "centroids");

	foreach_cell(g, [&](Graph::Vertex gv) -> bool {
		Dart d0 = value<Dart>(g, gAttribs.halfedge_volume_connection, Graph::HalfEdge(gv.dart));
		Vertex v0 = Vertex(d0);
		foreach_incident_face(m3, v0, [&](Face f) -> bool {
			if (trans_faces.is_marked(f))
			{
				value<Vec3>(m3, centroids, f) = {0, 0, 0};
				uint32 n = 0;
				foreach_incident_vertex(m3, f, [&](Vertex v) -> bool {
					value<Vec3>(m3, centroids, f) += value<Vec3>(m3, m3Attribs.vertex_position, v);
					++n;
					return true;
				});
				value<Vec3>(m3, centroids, f) /= n;
				cache_face2cut.add(f);
			}
			return true;
		});

		return true;
	});

	foreach_cell(m3, [&](Edge e) -> bool {
		uint32 nb_ft = 0;
		foreach_incident_face(m3, e, [&](Face f) -> bool {
			if (trans_faces.is_marked(f))
			{
				++nb_ft;
			}
			return true;
		});
		if (nb_ft >= 2)
			cache_edge2cut.add(e);

		return true;
	});

	foreach_cell(cache_edge2cut, [&](Edge e) -> bool {
		std::vector<Vertex> vertices = incident_vertices(m3, e);
		Vertex v_mid = cut_edge(m3, e);
		value<Vec3>(m3, m3Attribs.vertex_position, v_mid) =
			0.5 * (value<Vec3>(m3, m3Attribs.vertex_position, vertices[0]) +
				   value<Vec3>(m3, m3Attribs.vertex_position, vertices[1]));
		new_vertices.mark(v_mid);
		return true;
	});

	foreach_cell(cache_face2cut, [&](Face f) -> bool {
		Dart d0 = f.dart;
		bool isboundary = is_boundary(m3, d0);
		if (!is_boundary(m3, d0) && !vol2cut.is_marked(Volume(d0)))
		{
			vol2cut.mark(Volume(d0));
			cache_vol2cut.add(Volume(d0));
		}
		d0 = phi<3, 1>(m3, d0);
		if (!is_boundary(m3, d0) && !vol2cut.is_marked(Volume(d0)))
		{
			vol2cut.mark(Volume(d0));
			cache_vol2cut.add(Volume(d0));
		}

		Vec3 center = value<Vec3>(m3, centroids, f);
		Dart d1 = phi1(m3, f.dart);
		Edge e = cut_face(m3, Vertex(phi1(m3, f.dart)), Vertex(phi_1(m3, f.dart)));
		std::vector<Vertex> vertices = incident_vertices(m3, e);
		Vertex center_vertex = cut_edge(m3, e);
		value<Vec3>(m3, m3Attribs.vertex_position, center_vertex) = center;

		cut_face(m3, Vertex(phi<1, 1>(m3, d1)), Vertex(phi_1(m3, d1)));
		foreach_incident_face(m3, center_vertex, [&](Face new_face) -> bool {
			trans_faces.mark(new_face);
			return true;
		});
		return true;
	});

	foreach_cell(m3, [&](Face f) -> bool {
		if (!trans_faces.is_marked(f) && codegree(m3, f) == 6)
		{
			cache_sideface2cut.add(f);
		}
		return true;
	});

	std::vector<Vertex> vertices;
	vertices.reserve(2);
	foreach_cell(cache_sideface2cut, [&](Face f) -> bool {
		foreach_incident_vertex(m3, f, [&](Vertex v) -> bool {
			if (new_vertices.is_marked(v))
				vertices.push_back(v);
			return true;
		});
		cut_face(m3, vertices[0], vertices[1]);
		vertices.clear();
		return true;
	});

	std::vector<Dart> path0, path1;
	path0.reserve(6);
	path1.reserve(4);
	foreach_cell(cache_vol2cut, [&](Volume w) -> bool {
		Dart d = phi1(m3, w.dart);
		path0.push_back(d);
		d = phi1(m3, d);
		path0.push_back(d);
		d = phi<1, 2, 1>(m3, d);
		path0.push_back(d);
		d = phi<1, 2, 1>(m3, d);
		path0.push_back(d);
		d = phi1(m3, d);
		path0.push_back(d);
		d = phi<1, 2, 1>(m3, d);
		path0.push_back(d);

		d = phi<2, 1>(m3, path0[1]);

		cut_volume(m3, path0);
		cut_face(m3, Vertex(phi2(m3, path0[0])), Vertex(phi2(m3, path0[3])));

		path1.push_back(d);
		d = phi<1, 2, 1>(m3, d);
		path1.push_back(d);
		d = phi<1, 2, 1>(m3, d);
		path1.push_back(d);
		d = phi<1, 2, 1>(m3, d);
		path1.push_back(d);

		cut_volume(m3, path1);

		path0.clear();
		path1.clear();
		return true;
	});

	remove_attribute<CMap3::Face>(m3, centroids);
}

void get_loop_path(CMap3& m3, Dart d0, std::vector<Dart>& path)
{
	Dart d = d0;
	do
	{
		path.push_back(d);
		d = phi<1, 2, 1>(m3, d);
	} while (d != d0);
}
*/
bool find_inter_frame(const Graph& g, Graph::Vertex gv, const GAttributes& gAttribs, Mat3& frame)
{
	Scalar eps = M_PI / 6;
	uint32 nb_points = degree(g, gv);

	std::vector<Vec3> directions;
	Vec3 A = value<Vec3>(g, gAttribs.vertex_position, gv);
	foreach_adjacent_vertex_through_edge(g, gv, [&](Graph::Vertex gv1) -> bool {
		directions.push_back((value<Vec3>(g, gAttribs.vertex_position, gv1) - A).normalized());
		return true;
	});

	std::vector<Scalar> angles(nb_points * nb_points);
	bool end = false;
	for (uint32 i = 0; i < nb_points && !end; ++i)
	{
		for (uint32 j = i + 1; j < nb_points && !end; ++j)
		{
			Scalar ang = cgogn::geometry::angle(directions[i], directions[j]);
			bool in_range = (ang < eps && ang > -eps) || (ang < (M_PI + eps) && ang > (M_PI - eps)) ||
							(ang < (M_PI / 2 + eps) && ang > (M_PI / 2 - eps));

			if (in_range)
				angles[i * nb_points + j] = round(directions[i].dot(directions[j]));
			else
				end = true;
		}
	}

	if (end)
		return false;

	uint32 current_id = 0;
	std::vector<uint32> axis_id(nb_points);
	for (uint32 j = 0; j < nb_points; ++j)
	{
		axis_id[j] = 0xffffffff;
	}

	for (uint32 i = 0; i < nb_points; ++i)
	{
		if (axis_id[i] != 0xffffffff)
			continue;

		bool new_axis = true;
		for (uint32 j = i; j < nb_points && new_axis; ++j)
		{
			if (angles[i * nb_points + j] == -1)
			{
				new_axis = false;
				axis_id[i] = current_id;
				axis_id[j] = current_id;
				++current_id;
			}
		}
		if (new_axis)
			axis_id[i] = current_id++;
	}

	std::vector<uint32> axis_set = {0, 0, 0};
	frame = Mat3::Zero();
	for (uint32 j = 0; j < nb_points; ++j)
	{
		frame.col(axis_id[j]) += axis_set[axis_id[j]]++ == 0 ? directions[j] : -directions[j];
		frame.col(axis_id[j]).normalize();
	}
	if (current_id == 2)
		frame.col(2) = frame.col(0).cross(frame.col(1));

	Vec3 x_y = frame.col(0).cross(frame.col(1)).normalized();
	bool up = frame.col(2).dot(x_y) > 0;
	if (!up)
		frame.col(2) *= -1;

	return true;
}

} // namespace modeling

} // namespace cgogn
