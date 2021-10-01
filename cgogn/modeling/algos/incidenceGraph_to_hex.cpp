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

#include <cgogn/modeling/algos/incidenceGraph_to_hex.h>

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
#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>

namespace cgogn
{

namespace modeling
{

std::tuple<IGAttributes, M2Attributes, M3Attributes> incidenceGraph_to_hex(IncidenceGraph& ig, CMap2& m2/*, CMap3& m3*/)
{
	IncidenceGraphData igData;
	IGAttributes igAttribs;
	M2Attributes m2Attribs;
	M3Attributes m3Attribs;
	std::cout << " start" << std::endl;

	bool okay = get_incidenceGraph_data(ig, igData);
	std::cout << uint32(igData.intersections.size()) << " intersections" << std::endl;
	std::cout << uint32(igData.branches.size()) << " branches" << std::endl;
	std::cout << uint32(igData.efjunctures.size()) << " efjunctures" << std::endl;
	std::cout << uint32(igData.ffjunctures.size()) << " ffjunctures" << std::endl;
	std::cout << uint32(igData.leaflets.size()) << " leaflets" << std::endl;

	if (!okay)
		std::cout << "error incidenceGraph_to_hex: get_incidenceGraph_data" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): got incidenceGraph data" << std::endl;

	if(okay)
		okay = add_incidenceGraph_attributes(ig, igAttribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: add_incidenceGraph_attributes" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): added incidenceGraph attributes" << std::endl;

	if(okay)
		okay = compute_faces_geometry(ig, igData, igAttribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: compute_faces_geometry" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): computed incidence graph face geometries" << std::endl;


if (okay)
		okay = add_cmap2_attributes(m2, m2Attribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: add_cmap2_attributes" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): added cmap2 attributes" << std::endl;

	if (okay)
		okay = build_contact_surfaces(ig, igAttribs, igData, m2, m2Attribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: build_contact_surfaces" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): contact surfaces built" << std::endl;


	return {igAttribs, m2Attribs, m3Attribs};
}

bool add_incidenceGraph_attributes(IncidenceGraph& ig, IGAttributes& igAttribs)
{
	igAttribs.vertex_position = get_attribute<Vec3, IncidenceGraph::Vertex>(ig, "position");
	if (!igAttribs.vertex_position)
	{
		std::cout << "The incidence graph has no vertex position attribute" << std::endl;
		return false;
	}
	igAttribs.face_center = add_attribute<Vec3, IncidenceGraph::Face>(ig, "face_center");
	igAttribs.face_normal = add_attribute<Vec3, IncidenceGraph::Face>(ig, "face_normal");
	igAttribs.face_vertex_tangent = add_attribute<std::vector<Vec3>, IncidenceGraph::Face>(ig, "face_vertex_tangents");
	igAttribs.vertex_contact_surface = add_attribute<Dart, IncidenceGraph::Vertex>(ig, "vertex_contact_surface");
	igAttribs.halfedge_contact_surface_face = add_attribute<std::pair<Dart, Dart>, IncidenceGraph::Edge>(ig, "halfedge_contact_surface_face");
	
	return true;
}

bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs)
{
	m2Attribs.vertex_position = add_attribute<Vec3, CMap2::Vertex>(m2, "position");
	m2Attribs.dual_vertex_graph_branch = add_attribute<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>, CMap2::Vertex>(m2, "graph_branch");
	// m2Attribs.volume_center = add_attribute<Vec3, CMap2::Volume>(m2, "center");
	m2Attribs.volume_igvertex = add_attribute<IncidenceGraph::Vertex, CMap2::Volume>(m2, "gvertex");
	// m2Attribs.edge_mid = add_attribute<Vec3, CMap2::Edge>(m2, "edge_mid");
	// m2Attribs.halfedge_volume_connection = add_attribute<Dart, CMap2::HalfEdge>(m2, "volume_connection");
	// m2Attribs.ortho_scaffold = add_attribute<CMap2*, CMap2::Volume>(m2, "ortho_scaffold");
	return true;
}

std::pair<IncidenceGraph::Edge, IncidenceGraph::Edge> find_branch_extremities(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Edge e0, CellMarker<IncidenceGraph, IncidenceGraph::Edge>& cm)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	Edge firstEdge = e0;
	Edge currentEdge = e0;

	std::vector<Vertex> incidentVertices = incident_vertices(ig, e0);
	std::vector<Edge> incidentEdges;

	Vertex nextVertex = v0;
	uint32 vDegree;

	do {
		cm.mark(currentEdge);
		incidentVertices = incident_vertices(ig, currentEdge);
		nextVertex = incidentVertices[0].index_ != nextVertex.index_ ? incidentVertices[0] : incidentVertices[1];
		vDegree = degree(ig, nextVertex);
		if(vDegree == 2)
		{
			incidentEdges = incident_edges(ig, nextVertex);
			currentEdge = incidentEdges[0].index_ != currentEdge.index_ ? incidentEdges[0] : incidentEdges[1];
		}
	} while (vDegree == 2);

	return {firstEdge, currentEdge};
}

bool get_incidenceGraph_data(const IncidenceGraph& ig, IncidenceGraphData& ig_data)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	bool success = true;
	uint32 debugC = 0;
	foreach_cell(ig, [&](Vertex v) -> bool {
		// uint32 degree = pseudoDegree(ig, v);
		std::pair<uint32, uint32> info = pseudoDegree(ig, v);
	// std::cout << "incidenceGraph_data" << std::endl;
		
		if(info.first > 2) 
		{
			ig_data.intersections.push_back(v);
			return true;
		}

		if(info.first == 2)
		{
			if(info.second == 1)
				ig_data.efjunctures.push_back(v);
			else if(info.second == 0)
				ig_data.ffjunctures.push_back(v);
		}

		return true;
	});

	CellMarker<IncidenceGraph, Face> face_marker(ig);
	foreach_cell(ig, [&](Face f0) -> bool
	{
		if(face_marker.is_marked(f0))
			return true;
		face_marker.mark(f0);

		ig_data.leaflets.push_back(f0);
		std::vector<Face> toVisit = {f0};
		for(uint32 i = 0; i < toVisit.size(); ++i)
		{
			foreach_adjacent_face_through_edge(ig, toVisit[i], [&](Face f) -> bool {
				if(!face_marker.is_marked(f))
				{
					face_marker.mark(f);
					toVisit.push_back(f);
				}
			});
		}

		return true;
	});

	CellMarker<IncidenceGraph, Edge> edge_marker(ig);
	// foreach_cell(ig, [&](Vertex v) -> bool {

	// 	return true;
	// });
	for(Vertex v : ig_data.efjunctures)
	{
		Edge e0;
		foreach_incident_edge(ig, v, [&](Edge e) -> bool {
			if(degree(ig, e) == 0)
			{
				e0 = e;
				return false;
			}
			return true;
		});

		if(!edge_marker.is_marked(e0))
			ig_data.branches.push_back(find_branch_extremities(ig, v, e0, edge_marker));
	}

	for(Vertex v : ig_data.intersections)
	{
		foreach_incident_edge(ig, v, [&](Edge e) -> bool {
			if(degree(ig, e) == 0 && !edge_marker.is_marked(e))
				ig_data.branches.push_back(find_branch_extremities(ig, v, e, edge_marker));
			return true;
		});

	}

	return success;
}

/// Get a face-cell from each leaflet incident to the vertex
std::vector<IncidenceGraph::Face> incident_leaflets(const IncidenceGraph& ig, IncidenceGraph::Vertex v0)
{
	std::vector<IncidenceGraph::Face> leaflets;
	CellMarker<IncidenceGraph, IncidenceGraph::Face> face_marker(ig);

	foreach_incident_face(ig, v0, [&](IncidenceGraph::Face f0) -> bool {
		if(face_marker.is_marked(f0))
				return true;

		std::vector<IncidenceGraph::Face> leaflet = incident_leaflet(ig, v0, f0);
		for(IncidenceGraph::Face f : leaflet)
		{
			face_marker.mark(f);
		}
		if(leaflet.size())
			leaflets.push_back(leaflet[0]);

		return true;
	});

	return leaflets;
}

// Get all face-cells from a input leaflet f0
std::vector<IncidenceGraph::Face> incident_leaflet(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Face f0)
{
	std::vector<IncidenceGraph::Face> leaflet = {f0};
	CellMarker<IncidenceGraph, IncidenceGraph::Face> face_marker(ig);
	face_marker.mark(f0);

	for(uint32 i = 0; i < leaflet.size(); ++i)
	{
		foreach_adjacent_face_through_edge(ig, leaflet[i], [&](IncidenceGraph::Face f) -> bool {
			if(face_marker.is_marked(f))
				return true;
			face_marker.mark(f);

			if(contains_vertex(ig, v0, f))
				leaflet.push_back(f);

			return true;
		});
	}

	return leaflet;
}

// Get all edges belonging to the same leaflet as e0 that are also incident to v0
std::vector<IncidenceGraph::Edge> incident_leaflet_edges(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Edge e0)
{
	std::vector<IncidenceGraph::Edge> edges = {e0};
	CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);

	for(uint32 i = 0; i < edges.size(); ++i)
	{
		foreach_adjacent_edge_through_face(ig, edges[i], [&](IncidenceGraph::Edge e) -> bool {
			if(edge_marker.is_marked(e))
				return true;
			edge_marker.mark(e);

			std::vector<IncidenceGraph::Vertex> iv = incident_vertices(ig, e);

			if(v0.index_ == iv[0].index_ || v0.index_ == iv[1].index_)
				edges.push_back(e);
				
			return true;
		});
	}

	return edges;
}

/// Check if a face contains given vertex
bool contains_vertex(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Face f0)
{
	bool found = false;
	foreach_incident_vertex(ig, f0, [&](IncidenceGraph::Vertex v) -> bool {
		if(v.index_ == v0.index_)
			found = true;
		return !found;
	});

	return found;
}

void index_volume_cells_igh(CMap2& m, CMap2::Volume vol)
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

bool compute_faces_geometry(const IncidenceGraph& ig, const IncidenceGraphData& incidenceGraph_data, IGAttributes& igAttribs)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	foreach_cell(ig, [&](Face f0) -> bool {
		Vec3 center = {0, 0, 0};
		std::vector<Vec3> vertices;

		foreach_incident_vertex(ig, f0, [&](Vertex v) -> bool {
			Vec3 pos = value<Vec3>(ig, igAttribs.vertex_position, v);
			vertices.push_back(pos);
			center += pos;
		});
		center /= vertices.size();

		value<Vec3>(ig, igAttribs.face_center, f0) = center;
		
		value<Vec3>(ig, igAttribs.face_normal, f0) = {0, 0, 0};
		for(uint32 i = 9; i < vertices.size(); ++i)
		{
			Vec3 v0 = vertices[i] - center;
			Vec3 v1 = vertices[(i+1)%vertices.size()] - center;
			Vec3 n = v0.cross(v1);
			n.normalize();
			value<Vec3>(ig, igAttribs.face_normal, f0) += n;
		}
		value<Vec3>(ig, igAttribs.face_normal, f0).normalize();
	
		return true;
	});

	/// todo : face_vertex_tangents

	// for(IncidenceGraph::Vertex v : incidenceGraph_data.efjunctures)
	// {
	// 	Vec3 tangent = {0, 0, 0};
	// 	foreach_incident_face(ig, v, [&](IncidenceGraph::Face f) -> bool {

	// 		return true;
	// 	});
	// }

	return true;
}

bool build_contact_surfaces(const IncidenceGraph& ig, IGAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data, CMap2& m2, M2Attributes& m2Attribs)
{
	bool res = true;

	foreach_cell(ig, [&](IncidenceGraph::Vertex v) -> bool {
		std::pair<uint32, uint32> pd = pseudoDegree(ig, v);

		if (pd.second == 1)
		{
			build_contact_surface_1(ig, igAttribs, m2, m2Attribs, v);
			return res;
		}
		if (pd.first == 2)
		{
			build_contact_surface_2(ig, igAttribs, m2, m2Attribs, v);
			return res;
		}
		if (pd.first >= 3 && pd.first <= 6)
		{
			bool ortho = false;

			if(ortho)
				return res;
		}

		if (pd.first >= 3 )
			build_contact_surface_n(ig, igAttribs, m2, m2Attribs, v);

		return res;
	});


	return res;
}

void build_contact_surface_1(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, IncidenceGraph::Vertex v)
{
	Dart d = add_face(m2, 4, true).dart;

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = d;
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(d)) = v;

	IncidenceGraph::Edge e = incident_edges(ig, v)[0];
	std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, e);
	if(v.index_ == ivs[0].index_)
		value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first = d;
	else
		value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second = d;
}

void build_contact_surface_2(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, d1);
	phi2_sew(m2, phi1(m2, d0), phi_1(m2, d1));
	phi2_sew(m2, phi<11>(m2, d0), phi<11>(m2, d1));
	phi2_sew(m2, phi_1(m2, d0), phi1(m2, d1));

	index_volume_cells_igh(m2, CMap2::Volume(d0));

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = d0;
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(d0)) = v;

	CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);
	uint32 firstFace = 0;
	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e0) -> bool {
		std::vector<IncidenceGraph::Edge> leaflet = incident_leaflet_edges(ig, v, e0);
		for(IncidenceGraph::Edge e : leaflet)
		{
			std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, e);
			if(v.index_ == ivs[0].index_)
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first = firstFace == 0 ? d0 : d1;
			else
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second = firstFace == 0 ? d0 : d1;
		}
		++firstFace;
		return firstFace != 2;
	});
}

void build_contact_surface_orange(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs,
								  IncidenceGraph::Vertex v)
{

	uint32 nbf = pseudoDegree(ig, v).first;
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
		phi2_sew(m2, phi_1(m2, d), phi<11>(m2, e));
	}

	index_volume_cells_igh(m2, CMap2::Volume(faces[0]));

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = faces[0];

	// get the points on the sphere for each incident branch
	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);
	std::vector<Vec3> Ppos;
	Ppos.reserve(nbf);

	std::vector<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>> Pev;
	Pev.reserve(nbf);
	
	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
		std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, e);
		IncidenceGraph::Vertex v2 = ivs[0].index_ == v.index_ ? ivs[1] : ivs[0];
		Vec3 p = value<Vec3>(ig, igAttribs.vertex_position, v2);
		geometry::project_on_sphere(p, center, 1);
		Ppos.push_back(p);
		Pev.push_back({v, e});
		return true;

		// Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
		// geometry::project_on_sphere(p, center, 1);
		// Ppos.push_back(p);
		// Pdart.push_back(d);
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
	std::vector<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>> sorted_Pev(nbf);
	std::transform(permutation.begin(), permutation.end(), sorted_Ppos.begin(), [&](uint32 i) { return Ppos[i]; });
	std::transform(permutation.begin(), permutation.end(), sorted_Pev.begin(), [&](uint32 i) { return Pev[i]; });

	// put the geometry on the surface mesh vertices
	Vec3 Q1 = center + plane.first;
	Vec3 Q2 = center - plane.first;
	for (uint32 i = 0; i < nbf; ++i)
	{
		IncidenceGraph::Edge e = sorted_Pev[i].second;
		IncidenceGraph::Vertex v1 = sorted_Pev[i].first;
		std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, e);
		if(v1.index_ == ivs[0].index_)
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first = faces[i];
		else
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second = faces[i];
				
		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(faces[i])) =
			center + (sorted_Ppos[(i + 1) % nbf] - sorted_Ppos[i]).normalized().cross(plane.first);
	}
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, faces[0]))) = Q1;
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, faces[0]))) = Q2;

	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(faces[0])) = v;
}


void build_contact_surface_n(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, IncidenceGraph::Vertex v)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	uint32 nbf = pseudoDegree(ig, v).first;


	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

	std::vector<Vec3> Ppos;
	Ppos.reserve(nbf);

	foreach_incident_edge(ig, v, [&](Edge e) -> bool {
		std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, e);
		Vertex v2 = ivs[0].index_ == v.index_ ? ivs[1] : ivs[0];
		Vec3 p = value<Vec3>(ig, igAttribs.vertex_position, v2);
		geometry::project_on_sphere(p, center, 1);
		Ppos.push_back(p);
		return true;
	});


	std::pair<Vec3, Scalar> plane = geometry::plane_fitting(Ppos);
	bool planar = true;
	for (const Vec3& p : Ppos)
	{
		Scalar dist = geometry::distance_plane_point(plane.first, plane.second, p);
		if (dist > 0.05)
		{
			planar = false;
			break;
		}
	}
	if (planar)
	{
		build_contact_surface_orange(ig, igAttribs, m2, m2Attribs, v);
		return;
	}


	Dart vol_dart = convex_hull_around_vertex(ig, v, m2, m2Attribs, Ppos);


	index_volume_cells_igh(m2, CMap2::Volume(vol_dart));
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(vol_dart)) = v;

	vol_dart = remesh_igh(m2, CMap2::Volume(vol_dart), m2Attribs);
	dualize_volume(m2, CMap2::Volume(vol_dart), m2Attribs, ig, igAttribs);

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = vol_dart;

	return;
}

Dart convex_hull_around_vertex(const IncidenceGraph& ig, IncidenceGraph::Vertex v, CMap2& m2, M2Attributes& m2Attribs,
							   std::vector<Vec3>& Ppos)
{
	std::vector<uint32> Pid;
	Pid.reserve(Ppos.size());

	uint32 i = 0;
	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
		uint32 vertex_id = new_index<CMap2::Vertex>(m2);
		Pid.push_back(vertex_id);
		(*m2Attribs.vertex_position)[vertex_id] = Ppos[i]; // positions are in the same order in Ppos
		(*m2Attribs.dual_vertex_graph_branch)[vertex_id] = {v, e};
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

Vec3 slerp_igh(Vec3 A, Vec3 B, Scalar alpha, bool in)
{
	Scalar phi = geometry::angle(A, B) - (in ? 0 : 2 * M_PI);
	Scalar s0 = std::sin(phi * (1 - alpha));
	Scalar s1 = std::sin(phi * alpha);
	Scalar s2 = std::sin(phi);
	Vec3 sl = A * (s0 / s2);
	sl += B * (s1 / s2);
	return sl;
}

Scalar angle_on_sphere_igh(Vec3 A, Vec3 B, Vec3 C)
{
	Vec3 sB = slerp_igh(A, B, 0.01, true);
	Vec3 sC = slerp_igh(A, C, 0.01, true);
	Vec3 AB = sB - A;
	Vec3 AC = sC - A;
	return geometry::angle(AB, AC);
}

Scalar edge_max_angle_igh(CMap2& m2, CMap2::Edge e, M2Attributes& m2Attribs)
{
	Dart ed0 = e.dart;
	Dart ed1 = phi2(m2, ed0);

	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed0)));

	Scalar a0 = angle_on_sphere_igh(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed1)));
	a0 += angle_on_sphere_igh(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed1)));

	Scalar a1 = angle_on_sphere_igh(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed0)));
	a1 += angle_on_sphere_igh(A, B, C);

	return std::max(a0, a1);
}

Scalar min_cut_angle_igh(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, M2Attributes& m2Attribs)
{
	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v0.dart)));
	Scalar a0 = angle_on_sphere_igh(A, B, C);

	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v1.dart)));
	Scalar a1 = angle_on_sphere_igh(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	B = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v1.dart)));

	Scalar a2 = angle_on_sphere_igh(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v0.dart)));
	Scalar a3 = angle_on_sphere_igh(A, B, C);

	return std::min({a0, a1, a2, a3});
}

Vec3 spherical_barycenter_igh(std::vector<Vec3>& points, uint32 iterations)
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
			writevec[i] = slerp_igh(readvec[i], readvec[(i + 1) % nb_pts], Scalar(0.5), true);
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


Dart remesh_igh(CMap2& m2, CMap2::Volume vol, M2Attributes& m2Attribs)
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
				value<Scalar>(m2, edge_angle_max, e) = edge_max_angle_igh(m2, e, m2Attribs);
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
							{{verts_3[i], verts_3[j]}, min_cut_angle_igh(m2, verts_3[i], verts_3[j], m2Attribs)});
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

				dijkstra_topo_igh(m2, v, previous, dist);

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
				angles.push_back({edge_max_angle_igh(m2, first_edges[0], m2Attribs), 0});
				angles.push_back({edge_max_angle_igh(m2, first_edges[1], m2Attribs), 1});
				angles.push_back({edge_max_angle_igh(m2, first_edges[2], m2Attribs), 2});
				angles.push_back({edge_max_angle_igh(m2, first_edges[3], m2Attribs), 3});

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
					v1 = CMap2::Vertex(phi<211>(m2, path.front()));
					e0 = CMap2::Edge(phi<21>(m2, path.front()));
					path.erase(path.begin());
					break;

				case 2:
					v0 = CMap2::Vertex(path.back());
					v1 = CMap2::Vertex(phi<11>(m2, path.back()));
					e0 = CMap2::Edge(phi1(m2, path.back()));
					path.pop_back();
					break;

				case 3:
					v0 = CMap2::Vertex(phi_1(m2, phi2(m2, path.back())));
					v1 = CMap2::Vertex(phi<21>(m2, path.back()));
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


void dualize_volume(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs, const IncidenceGraph& ig, IGAttributes& igAttribs)
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
	const IncidenceGraph::Vertex gv = value<IncidenceGraph::Vertex>(m, m2Attribs.volume_igvertex, vol);
	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, gv);
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
			std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge> branch = 
				value<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>>(m, m2Attribs.dual_vertex_graph_branch, CMap2::Vertex(d));
			// store the contact surface face on the outgoing branch
			std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, branch.second);
			if(branch.first.index_ == ivs[0].index_)
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, branch.second).first = d;
			else
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, branch.second).second = d;

			return true;
		}
		return true;
	});

	// DartMarkerStore<CMap2> vertex_marker(m);
	// foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
	// 	if (!vertex_marker.is_marked(d))
	// 	{
	// 		CMap2::Vertex v(d);

	// 		std::vector<Vec3> points;
	// 		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
	// 			vertex_marker.mark(d);
	// 			points.push_back(value<Vec3>(m, m2Attribs.vertex_position, CMap2::Vertex(d)) - center);
	// 			return true;
	// 		});
	// 		Vec3 b;
	// 		if (points.size() == 2)
	// 			b = slerp(points[0], points[1], 0.5, true) + center;
	// 		else
	// 			b = spherical_barycenter(points, 10) + center;

	// 		set_index(m, v, new_index<CMap2::Vertex>(m));	  // give a new index to the vertex
	// 		value<Vec3>(m, m2Attribs.vertex_position, v) = b; // set the position to the computed position
	// 		return true;
	// 	}
	// 	return true;
	// });
}

bool dijkstra_topo_igh(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
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


} // namespace modeling

} // namespace cgogn
