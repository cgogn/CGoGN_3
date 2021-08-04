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
	// m2Attribs.dual_vertex_graph_branch = add_attribute<Dart, CMap2::Vertex>(m2, "graph_branch");
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
	std::cout << "incidenceGraph_data " << 0 << std::endl;
	foreach_cell(ig, [&](Vertex v) -> bool {
		// uint32 degree = pseudoDegree(ig, v);
		std::pair<uint32, uint32> info = pseudoDegree(ig, v);
		std::cout << "incidenceGraph_data " << 1 << std::endl;
	// std::cout << "incidenceGraph_data" << std::endl;
		
		if(info.first > 2) 
		{
			ig_data.intersections.push_back(v);
			return true;
		}
		std::cout << "incidenceGraph_data " << 2 << std::endl;

		if(info.first == 2)
		{
			if(info.second == 1)
				ig_data.efjunctures.push_back(v);
			else if(info.second == 0)
				ig_data.ffjunctures.push_back(v);
		}

		return true;
	});
	std::cout << "incidenceGraph_data " << 3 << std::endl;

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
	std::cout << "incidenceGraph_data " << 4 << std::endl;

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

	std::cout << "incidenceGraph_data " << 5 << std::endl;

	for(Vertex v : ig_data.intersections)
	{
		foreach_incident_edge(ig, v, [&](Edge e) -> bool {
			if(degree(ig, e) == 0 && !edge_marker.is_marked(e))
				ig_data.branches.push_back(find_branch_extremities(ig, v, e, edge_marker));
			return true;
		});

	}
	std::cout << "incidenceGraph_data " << 6 << std::endl;

	return success;
}

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
		if (pd.first == 1)
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

	index_volume_cells(m2, CMap2::Volume(d0));

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

void build_contact_surface_n(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, IncidenceGraph::Vertex v)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	uint32 nbf = degree(ig, v);

	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

	std::vector<Vec3> Ppos;
	Ppos.reserve(nbf);

	// foreach_dart_of_orbit(g, v, [&](Dart d) -> bool {
	// 	Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(alpha0(g, d)));
	// 	geometry::project_on_sphere(p, center, 1);
	// 	Ppos.push_back(p);
	// 	return true;
	// });
	foreach_incident_edge(ig, v, [&](Edge e) -> bool {
		std::vector<IncidenceGraph::Vertex> ivs = incident_vertices(ig, e);
		Vertex v2 = ivs[0].index_ == v.index_ ? ivs[1] : ivs[0];
		Vec3 p = value<Vec3>(ig, igAttribs.vertex_position, v2);
		geometry::project_on_sphere(p, center, 1);
		Ppos.push_back(p);
		return true;
	});

	// std::pair<Vec3, Scalar> plane = geometry::plane_fitting(Ppos);
	// bool planar = true;
	// for (const Vec3& p : Ppos)
	// {
	// 	Scalar dist = geometry::distance_plane_point(plane.first, plane.second, p);
	// 	if (dist > 0.15)
	// 	{
	// 		planar = false;
	// 		break;
	// 	}
	// }
	// if (planar)
	// {
	// 	build_contact_surface_orange(g, gAttribs, m2, m2Attribs, v);
	// 	return;
	// }


	Dart vol_dart = convex_hull_around_vertex(g, v, m2, m2Attribs, Ppos);

	index_volume_cells(m2, CMap2::Volume(vol_dart));
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_gvertex, CMap2::Volume(vol_dart)) = v;

	vol_dart = remesh(m2, CMap2::Volume(vol_dart), m2Attribs);
	dualize_volume(m2, CMap2::Volume(vol_dart), m2Attribs, ig, igAttribs);

	value<Dart>(g, gAttribs.vertex_contact_surface, v) = vol_dart;

	return;
}

Dart convex_hull_around_vertex(const IncidenceGraph& g, IncidenceGraph::Vertex v, CMap2& m2, M2Attributes& m2Attribs,
							   std::vector<Vec3>& Ppos)
{
	// std::vector<uint32> Pid;
	// Pid.reserve(Ppos.size());

	// uint32 i = 0;
	// foreach_dart_of_orbit(g, v, [&](Dart d) -> bool {
	// 	uint32 vertex_id = new_index<CMap2::Vertex>(m2);
	// 	Pid.push_back(vertex_id);
	// 	(*m2Attribs.vertex_position)[vertex_id] = Ppos[i]; // positions are in the same order in Ppos
	// 	(*m2Attribs.dual_vertex_graph_branch)[vertex_id] = d;
	// 	return true;
	// });

	// std::vector<uint32> faces_nb_vertices;
	// std::vector<uint32> faces_vertex_indices;

	// std::vector<uint32> indices;
	// for (uint32 i = 0; i < Ppos.size() - 2; ++i)
	// {
	// 	for (uint32 j = i + 1; j < Ppos.size() - 1; ++j)
	// 	{
	// 		for (uint32 k = j + 1; k < Ppos.size(); ++k)
	// 		{
	// 			Vec3 t0 = Ppos[j] - Ppos[i];
	// 			Vec3 t1 = Ppos[k] - Ppos[i];
	// 			Vec3 n = t0.cross(t1);
	// 			int sign = 0;

	// 			for (uint32 m = 0; m < Ppos.size(); ++m)
	// 			{
	// 				if (m == i || m == j || m == k)
	// 					continue;

	// 				Vec3 vec = Ppos[m] - Ppos[i];
	// 				Scalar d = vec.dot(n);

	// 				if (!sign)
	// 					sign = (d < 0 ? -1 : 1);
	// 				else
	// 				{
	// 					if (sign != (d < 0 ? -1 : 1))
	// 					{
	// 						sign = 0;
	// 						break;
	// 					}
	// 				}
	// 			}

	// 			if (sign != 0)
	// 			{
	// 				if (sign == 1)
	// 				{
	// 					indices.push_back(Pid[j]);
	// 					indices.push_back(Pid[i]);
	// 					indices.push_back(Pid[k]);
	// 				}
	// 				else
	// 				{
	// 					indices.push_back(Pid[i]);
	// 					indices.push_back(Pid[j]);
	// 					indices.push_back(Pid[k]);
	// 				}

	// 				faces_nb_vertices.push_back(3);
	// 				faces_vertex_indices.insert(faces_vertex_indices.end(), indices.begin(), indices.end());
	// 				indices.clear();
	// 			}
	// 		}
	// 	}
	// }

	// auto darts_per_vertex = add_attribute<std::vector<Dart>, CMap2::Vertex>(m2, "__darts_per_vertex");

	// uint32 faces_vertex_index = 0u;
	// std::vector<uint32> vertices_buffer;
	// vertices_buffer.reserve(16u);

	// Dart volume_dart;
	// std::vector<Dart> all_darts;
	// for (std::size_t i = 0u, end = faces_nb_vertices.size(); i < end; ++i)
	// {
	// 	uint32 nbv = faces_nb_vertices[i];

	// 	vertices_buffer.clear();

	// 	for (uint32 j = 0u; j < nbv; ++j)
	// 	{
	// 		uint32 idx = faces_vertex_indices[faces_vertex_index++];
	// 		vertices_buffer.push_back(idx);
	// 	}

	// 	if (nbv > 2u)
	// 	{
	// 		CMap1::Face f = add_face(static_cast<CMap1&>(m2), nbv, false);
	// 		Dart d = f.dart;
	// 		for (uint32 j = 0u; j < nbv; ++j)
	// 		{
	// 			all_darts.push_back(d);
	// 			const uint32 vertex_index = vertices_buffer[j];
	// 			set_index<CMap2::Vertex>(m2, d, vertex_index);
	// 			(*darts_per_vertex)[vertex_index].push_back(d);
	// 			d = phi1(m2, d);
	// 		}
	// 	}
	// }

	// for (Dart d : all_darts)
	// {
	// 	if (phi2(m2, d) == d)
	// 	{
	// 		uint32 vertex_index = index_of(m2, CMap2::Vertex(d));

	// 		const std::vector<Dart>& next_vertex_darts =
	// 			value<std::vector<Dart>>(m2, darts_per_vertex, CMap2::Vertex(phi1(m2, d)));
	// 		bool phi2_found = false;

	// 		for (auto it = next_vertex_darts.begin(); it != next_vertex_darts.end() && !phi2_found; ++it)
	// 		{
	// 			if (index_of(m2, CMap2::Vertex(phi1(m2, *it))) == vertex_index)
	// 			{
	// 				if (phi2(m2, *it) == *it)
	// 				{
	// 					phi2_sew(m2, d, *it);
	// 					phi2_found = true;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	// remove_attribute<CMap2::Vertex>(m2, darts_per_vertex);
	// return all_darts[0];
}

void dualize_volume(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs, const IncidenceGraph& ig, IGAttributes& igAttribs)
{
	// 	// set the new phi1
	// foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
	// 	Dart dd = phi2(m, phi_1(m, d));
	// 	(*(m.phi1_))[d.index] = dd;
	// 	return true;
	// });

	// // set the new phi_1
	// foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
	// 	Dart next = phi1(m, d);
	// 	(*(m.phi_1_))[next.index] = d;
	// 	return true;
	// });

	// DartMarkerStore<CMap2> face_marker(m);
	// const Graph::Vertex gv = value<Graph::Vertex>(m, m2Attribs.volume_gvertex, vol);
	// const Vec3& center = value<Vec3>(g, gAttribs.vertex_position, gv);
	// foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
	// 	if (!face_marker.is_marked(d))
	// 	{
	// 		CMap2::Face f(d);
	// 		foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
	// 			face_marker.mark(d);
	// 			return true;
	// 		});
	// 		if (is_indexed<CMap2::Face>(m))
	// 			set_index(m, f, new_index<CMap2::Face>(m)); // give a new index to the face
	// 		// darts of the new face orbit are darts of the old dual vertex
	// 		// get the outgoing graph branch from the old dual vertex
	// 		Dart branch = value<Dart>(m, m2Attribs.dual_vertex_graph_branch, CMap2::Vertex(d));
	// 		// store the contact surface face on the outgoing branch
	// 		value<Dart>(g, gAttribs.halfedge_contact_surface_face, Graph::HalfEdge(branch)) = d;
	// 		return true;
	// 	}
	// 	return true;
	// });

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

} // namespace modeling

} // namespace cgogn
