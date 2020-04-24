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

#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/core/types/cmap/cmap_ops.h>

namespace cgogn
{

namespace modeling
{

///////////
// CMaps //
///////////

bool graph_to_hex(Graph& g, CMap2& m2, CMap3& m3)
{
	bool okay;

	GData gData;
	GAttributes gAttribs;
	M2Attributes m2Attribs;

	okay = get_graph_data(g, gData);
	std::cout << uint32(gData.intersections.size()) << " inters" << std::endl;
	std::cout << uint32(gData.branches.size()) << " branches" << std::endl;
	for (auto b : gData.branches)
		std::cout << index_of(g, Graph::Vertex(b.first.dart)) << " - " << index_of(g, Graph::Vertex(b.second.dart))
				  << std::endl;
	if (!okay)
	{
		std::cout << "error graph_to_hex: get_graph_data" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): got graph data" << std::endl;

	okay = add_graph_attributes(g, gAttribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: add_graph_attributes" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): added graph attributes" << std::endl;

	okay = add_cmap2_attributes(m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: add_cmap2_attributes" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): added cmap2 attributes" << std::endl;

	okay = build_contact_surfaces(g, gAttribs, m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: build_contact_surfaces" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): intersection cores built" << std::endl;

	okay = create_intersection_frames(g, gAttribs, m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: create_intersections_frames" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): create_intersections_frames completed" << std::endl;

	okay = propagate_frames(g, gAttribs, gData, m2);
	if (!okay)
	{
		std::cout << "error graph_to_hex: propagate_frames" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): propagate_frames completed" << std::endl;

	okay = set_contact_surfaces_geometry(g, gAttribs, m2, m2Attribs);
	if (!okay)
	{
		std::cout << "error graph_to_hex: set_contact_surfaces_geometry" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): set_contact_surfaces_geometry completed" << std::endl;

	okay = build_branch_sections(g, gAttribs, m2, m2Attribs, m3);
	if (!okay)
	{
		std::cout << "error graph_to_hex: build_branch_sections" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): build_branch_sections completed" << std::endl;

	okay = sew_branch_sections(m2, m2Attribs, m3);
	if (!okay)
	{
		std::cout << "error graph_to_hex: sew_sections" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): sew_sections completed" << std::endl;

	okay = set_volumes_geometry(m2, m2Attribs, m3);
	if (!okay)
	{
		std::cout << "error graph_to_hex: set_volumes_geometry" << std::endl;
		return false;
	}
	else
		std::cout << "graph_to_hex (/): set_volumes_geometry completed" << std::endl;

	return okay;
}

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

void index_volume_cells(CMap2& m, CMap2::Volume vol)
{
	if (is_indexed<CMap2::Vertex>(m))
	{
		foreach_incident_vertex(m, vol, [&](CMap2::Vertex v) -> bool {
			set_index(m, v, new_index<CMap2::Vertex>(m));
			return true;
		});
	}
	if (is_indexed<CMap2::HalfEdge>(m))
	{
		foreach_incident_edge(m, vol, [&](CMap2::Edge e) -> bool {
			set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
			set_index(m, CMap2::HalfEdge(phi2(m, e.dart)), new_index<CMap2::HalfEdge>(m));
			return true;
		});
	}
	if (is_indexed<CMap2::Edge>(m))
	{
		foreach_incident_edge(m, vol, [&](CMap2::Edge e) -> bool {
			set_index(m, e, new_index<CMap2::Edge>(m));
			return true;
		});
	}
	if (is_indexed<CMap2::Face>(m))
	{
		foreach_incident_face(m, vol, [&](CMap2::Face f) -> bool {
			set_index(m, f, new_index<CMap2::Face>(m));
			return true;
		});
	}
	if (is_indexed<CMap2::Volume>(m))
		set_index(m, vol, new_index<CMap2::Volume>(m));
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

Graph::HalfEdge branch_extremity(const Graph& g, Graph::HalfEdge v1, CellMarker<Graph, Graph::Edge>& cm)
{
	Dart d = v1.dart;
	while (degree(g, Graph::Vertex(d)) == 2)
	{
		d = alpha0(g, alpha1(g, d));
		cm.mark(Graph::Edge(d));
	}
	return Graph::HalfEdge(d);
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

void project_on_sphere(Vec3& P, const Vec3& C, Scalar R)
{
	P = C + (P - C).normalized() * R;
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

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool get_graph_data(const Graph& g, GData& gData)
{
	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		if (degree(g, v) > 2)
			gData.intersections.push_back(v);
		return true;
	});

	CellMarker<Graph, Graph::Edge> cm(g);
	foreach_cell(g, [&](Graph::Edge e) -> bool {
		if (cm.is_marked(e))
			return true;
		cm.mark(e);
		std::vector<Graph::HalfEdge> halfedges = incident_halfedges(g, e);
		gData.branches.push_back({branch_extremity(g, halfedges[0], cm), branch_extremity(g, halfedges[1], cm)});
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
	if (!gAttribs.vertex_contact_surface)
	{
		std::cout << "Failed to add vertex_contact_surface attribute to graph" << std::endl;
		return false;
	}

	gAttribs.halfedge_contact_surface_face = add_attribute<Dart, Graph::HalfEdge>(g, "contact_surface_face");
	if (!gAttribs.halfedge_contact_surface_face)
	{
		std::cout << "Failed to add halfedge_contact_surface_face attribute to graph" << std::endl;
		return false;
	}

	gAttribs.halfedge_frame = add_attribute<Mat3, Graph::HalfEdge>(g, "frame");
	if (!gAttribs.halfedge_frame)
	{
		std::cout << "Failed to add halfedge_frame attribute to graph" << std::endl;
		return false;
	}

	return true;
}

bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs)
{
	m2Attribs.vertex_position = add_attribute<Vec3, CMap2::Vertex>(m2, "position");
	if (!m2Attribs.vertex_position)
	{
		std::cout << "Failed to add position attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.volume_center = add_attribute<Vec3, CMap2::Volume>(m2, "center");
	if (!m2Attribs.volume_center)
	{
		std::cout << "Failed to add volume_center attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.edge_mid = add_attribute<Vec3, CMap2::Edge>(m2, "edge_mid");
	if (!m2Attribs.edge_mid)
	{
		std::cout << "Failed to add edge_mid attribute to map2" << std::endl;
		return false;
	}

	m2Attribs.halfedge_volume_connection = add_attribute<Dart, CMap2::HalfEdge>(m2, "volume_connection");
	if (!m2Attribs.halfedge_volume_connection)
	{
		std::cout << "Failed to add halfedge_volume_connection attribute to map2" << std::endl;
		return false;
	}

	return true;
}

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs)
{
	bool res = true;
	parallel_foreach_cell(g, [&](Graph::Vertex v) -> bool {
		switch (degree(g, v))
		{
		case 1:
			build_contact_surface_1(g, gAttribs, m2, v);
			break;
		case 2:
			build_contact_surface_2(g, gAttribs, m2, v);
			break;
		case 3:
			build_contact_surface_3(g, gAttribs, m2, m2Attribs, v);
			break;
		default:
			std::cout << "build_contact_surfaces: intersections of degree > 3 not managed yet" << std::endl;
			res = false;
			break;
		}
		return res;
	});
	return res;
}

void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::Vertex v)
{
	Dart d = add_face(m2, 4, true).dart;

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d;

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(v.dart)) = d;
}

void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::Vertex v)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, d1);
	phi2_sew(m2, phi1(m2, d0), phi_1(m2, d1));
	phi2_sew(m2, phi<11>(m2, d0), phi<11>(m2, d1));
	phi2_sew(m2, phi_1(m2, d0), phi1(m2, d1));

	index_volume_cells(m2, CMap2::Volume(d0));

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d0;

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(v.dart)) = d0;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(alpha1(g, v.dart))) = phi<12>(m2, d0);
}

void build_contact_surface_3(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, Graph::Vertex v)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d2 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, phi1(m2, d1));
	phi2_sew(m2, phi1(m2, d0), d2);
	phi2_sew(m2, phi<11>(m2, d0), phi_1(m2, d2));
	phi2_sew(m2, phi_1(m2, d0), phi<11>(m2, d1));
	phi2_sew(m2, d1, phi1(m2, d2));
	phi2_sew(m2, phi_1(m2, d1), phi<11>(m2, d2));

	index_volume_cells(m2, CMap2::Volume(d0));

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = d0;

	const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, v);
	Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);

	std::vector<Vec3> Ppos;
	Ppos.reserve(3);
	std::vector<Dart> Pdart;
	Pdart.reserve(3);

	foreach_dart_of_orbit(g, v, [&](Dart d) {
		Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		project_on_sphere(p, center, radius);
		Ppos.push_back(p);
		Pdart.push_back(d);
		return true;
	});

	Vec3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
	std::vector<Vec3> Q{center + V * radius, center - V * radius};
	std::vector<Vec3> M{center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
						center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
						center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};

	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(Pdart[0])) = d0;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(Pdart[1])) = d1;
	value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(Pdart[2])) = d2;

	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d0)) = M[0];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d1)) = M[1];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d2)) = M[2];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, d0))) = Q[0];
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, d0))) = Q[1];
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
		Dart d1 = phi<11>(m2, d0);

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

bool propagate_frames(const Graph& g, GAttributes& gAttribs, const GData& gData, CMap2& m2)
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
			propagate_frame_n_1(g, gAttribs, branch.second);
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
		Scalar cos = UE.col(0).dot(U_.col(0));
		Scalar angle = cos > 1 ? std::acos(1) : std::acos(cos);
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
			Eigen::AngleAxisd rot(angle_step * step, U.col(2));
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
		value<Vec3>(m2, m2Attribs.volume_center, contact_surface) = center;

		Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, v);

		if (degree(g, v) < 3)
		{
			Graph::HalfEdge h(v.dart);
			Dart csf = value<Dart>(g, gAttribs.halfedge_contact_surface_face, h);
			Mat3 frame = value<Mat3>(g, gAttribs.halfedge_frame, h);

			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(csf)) = center - frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, csf))) = center + frame.col(0) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi<11>(m2, csf))) =
				center + frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, csf))) = center - frame.col(0) * radius;
		}

		foreach_incident_edge(m2, contact_surface, [&](CMap2::Edge e) -> bool {
			std::vector<CMap2::Vertex> vertices = incident_vertices(m2, e);
			Vec3 mid = 0.5 * (value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) +
							  value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]));
			project_on_sphere(mid, center, radius);
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

bool build_branch_sections(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	parallel_foreach_cell(g, [&](Graph::Edge e) -> bool {
		std::vector<Graph::HalfEdge> halfedges = incident_halfedges(g, e);

		Dart m2f0 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, halfedges[0]);
		Dart m2f1 = value<Dart>(g, gAttribs.halfedge_contact_surface_face, halfedges[1]);

		std::vector<Dart> F0 = {m2f0, phi1(m2, m2f0), phi<11>(m2, m2f0), phi_1(m2, m2f0)};
		std::vector<Dart> F1 = {m2f1, phi1(m2, m2f1), phi<11>(m2, m2f1), phi_1(m2, m2f1)};

		Dart m3d = add_branch_section(m3);
		std::vector<Dart> D0 = {m3d, phi<2321>(m3, m3d), phi<23212321>(m3, m3d), phi<111232>(m3, m3d)};
		std::vector<Dart> D1 = {phi<2112>(m3, D0[0]), phi<2112>(m3, D0[1]), phi<2112>(m3, D0[2]), phi<2112>(m3, D0[3])};

		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);

		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[0])) = phi<11>(m3, D1[1]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[1])) = phi<11>(m3, D1[0]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[2])) = phi<11>(m3, D1[3]);
		value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[3])) = phi<11>(m3, D1[2]);

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

bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3)
{
	auto m3pos = cgogn::add_attribute<Vec3, CMap3::Vertex>(m3, "position");

	parallel_foreach_cell(m2, [&](CMap2::Volume v) -> bool {
		Dart m3d = phi_1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(v.dart)));
		value<Vec3>(m3, m3pos, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.volume_center, v);
		return true;
	});

	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		Dart m3d = phi1(m3, value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(e.dart)));
		value<Vec3>(m3, m3pos, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.edge_mid, e);
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
				value<Vec3>(m3, m3pos, m3v) = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, m2d)));
			}
		}
	}
	return true;
}

} // namespace modeling

} // namespace cgogn
