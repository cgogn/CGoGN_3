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

#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// add_edge(MESH& m, typename mesh_traits<MESH>::Vertex v0, typename mesh_traits<MESH>::Vertex v1);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

IncidenceGraph::Edge add_edge(IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1)
{
	using Edge = IncidenceGraph::Edge;

	Edge e = add_cell<Edge>(ig);
	(*ig.edge_incident_vertices_)[e.index_] = {v0, v1};
	(*ig.edge_incident_faces_)[e.index_].clear();
	(*ig.vertex_incident_edges_)[v0.index_].push_back(e);
	(*ig.vertex_incident_edges_)[v1.index_].push_back(e);

	return e;
}

/*****************************************************************************/

// template <typename MESH>
// void
// remove_edge(MESH& m, typename mesh_traits<MESH>::Edge e);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

void remove_edge(IncidenceGraph& ig, IncidenceGraph::Edge e)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	while ((*ig.edge_incident_faces_)[e.index_].size() > 0)
		remove_face(ig, (*ig.edge_incident_faces_)[e.index_].back());

	auto [v0, v1] = (*ig.edge_incident_vertices_)[e.index_];
	remove_edge_in_vertex(ig, v0, e);
	remove_edge_in_vertex(ig, v1, e);

	remove_cell<Edge>(ig, e);
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// cut_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

IncidenceGraph::Vertex cut_edge(IncidenceGraph& ig, IncidenceGraph::Edge e0, bool /*set_indices*/)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	auto [v0, v1] = (*ig.edge_incident_vertices_)[e0.index_];
	Vertex v = add_cell<Vertex>(ig);
	(*ig.edge_incident_vertices_)[e0.index_] = {v0, v};
	Edge e1 = add_edge(ig, v, v1);
	for (Face f : (*ig.edge_incident_faces_)[e0.index_])
	{
		(*ig.face_incident_edges_)[f.index_].push_back(e1);
		sort_face_edges(ig, f); // TODO: could do more efficient (insert)
	}
	return v;
}

///////////
// Graph //
///////////

Graph::Vertex cut_edge(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart e0 = e.dart;
	Dart e1 = alpha0(g, e0);

	Dart v0 = add_dart(g);
	Dart v1 = add_dart(g);

	alpha1_sew(g, v0, v1);
	alpha0_unsew(g, e0);
	alpha0_sew(g, e0, v0);
	alpha0_sew(g, e1, v1);

	if (set_indices)
	{
		if (is_indexed<Graph::Vertex>(g))
			set_index(g, Graph::Vertex(v0), new_index<Graph::Vertex>(g));
		if (is_indexed<Graph::HalfEdge>(g))
		{
			set_index(g, Graph::HalfEdge(v0), new_index<Graph::HalfEdge>(g));
			set_index(g, Graph::HalfEdge(v1), new_index<Graph::HalfEdge>(g));
		}
		if (is_indexed<Graph::Edge>(g))
		{
			copy_index<Graph::Edge>(g, v0, e0);
			set_index(g, Graph::Edge(e1), new_index<Graph::Edge>(g));
		}
	}

	return Graph::Vertex(v0);
}

///////////
// CMap1 //
///////////

CMap1::Vertex cut_edge(CMap1& m, CMap1::Edge e, bool set_indices)
{
	Dart d = add_dart(m);
	phi1_sew(m, e.dart, d);
	CMap1::Vertex v(d);

	if (set_indices)
	{
		if (is_indexed<CMap1::Vertex>(m))
			set_index(m, v, new_index<CMap1::Vertex>(m));
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (is_indexed<CMap1::Face>(m))
			copy_index<CMap1::Face>(m, d, e.dart);
	}

	return v;
}

///////////
// CMap2 //
///////////

CMap2::Vertex cut_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart d1 = e.dart;
	Dart d2 = phi2(m, d1);
	phi2_unsew(m, d1);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d1), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d2), false);
	phi2_sew(m, d1, nv2.dart);
	phi2_sew(m, d2, nv1.dart);
	set_boundary(m, nv1.dart, is_boundary(m, d1));
	set_boundary(m, nv2.dart, is_boundary(m, d2));
	CMap2::Vertex v(nv1.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
			set_index(m, v, new_index<CMap2::Vertex>(m));
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			set_index(m, CMap2::HalfEdge(nv1.dart), new_index<CMap2::HalfEdge>(m));
			set_index(m, CMap2::HalfEdge(nv2.dart), new_index<CMap2::HalfEdge>(m));
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			copy_index<CMap2::Edge>(m, nv2.dart, d1);
			set_index(m, CMap2::Edge(nv1.dart), new_index<CMap2::Edge>(m));
		}
		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m, nv1.dart, d1);
			copy_index<CMap2::Face>(m, nv2.dart, d2);
		}
		if (is_indexed<CMap2::Volume>(m))
		{
			copy_index<CMap2::Volume>(m, nv1.dart, d1);
			copy_index<CMap2::Volume>(m, nv2.dart, d2);
		}
	}

	return v;
}

///////////
// CMap3 //
///////////

CMap3::Vertex cut_edge(CMap3& m, CMap3::Edge e, bool set_indices)
{
	Dart d0 = e.dart;
	Dart d23 = phi<2, 3>(m, d0);

	CMap3::Vertex v(cut_edge(static_cast<CMap2&>(m), CMap2::Edge(d0), false).dart);

	while (d23 != e.dart)
	{
		d0 = d23;
		d23 = phi<2, 3>(m, d23);

		cut_edge(static_cast<CMap2&>(m), CMap2::Edge(d0), false);

		const Dart d3 = phi3(m, d0);
		phi3_unsew(m, d0);

		phi3_sew(m, d0, phi1(m, d3));
		phi3_sew(m, d3, phi1(m, d0));
	}

	const Dart d3 = phi3(m, e.dart);
	phi3_unsew(m, e.dart);

	phi3_sew(m, e.dart, phi1(m, d3));
	phi3_sew(m, d3, phi1(m, e.dart));

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
			set_index(m, v, new_index<CMap3::Vertex>(m));
		if (is_indexed<CMap3::Vertex2>(m))
		{
			Dart d = v.dart;
			do
			{
				if (!is_boundary(m, d))
					set_index(m, CMap3::Vertex2(d), new_index<CMap3::Vertex2>(m));
				d = phi<2, 3>(m, d);
			} while (d != v.dart);
		}
		if (is_indexed<CMap3::Edge>(m))
		{
			set_index(m, CMap3::Edge(v.dart), new_index<CMap3::Edge>(m));
			set_index(m, e, index_of(m, e));
		}
		if (is_indexed<CMap3::Face2>(m))
		{
			Dart d = e.dart;
			do
			{
				// if (!is_boundary(m, d))
				copy_index<CMap3::Face2>(m, phi1(m, d), d);
				// if (!is_boundary(m, phi3(m, d)))
				copy_index<CMap3::Face2>(m, phi3(m, d), phi<1, 3>(m, d));
				d = phi<2, 3>(m, d);
			} while (d != e.dart);
		}
		if (is_indexed<CMap3::Face>(m))
		{
			Dart d = e.dart;
			do
			{
				copy_index<CMap3::Face>(m, phi1(m, d), d);
				copy_index<CMap3::Face>(m, phi3(m, d), d);
				// copy_index<CMap3::Face>(m, phi2(m, d), phi<1, 2>(m, d));
				d = phi<2, 3>(m, d);
			} while (d != e.dart);
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			foreach_dart_of_orbit(m, e, [&](Dart d) -> bool {
				if (is_boundary(m, d))
					return true;
				copy_index<CMap3::Volume>(m, phi1(m, d), d);
				copy_index<CMap3::Volume>(m, phi2(m, d), d);
				return true;
			});
		}
	}

	return v;
}

//////////
// CPH3 //
//////////

CPH3::CMAP::Vertex cut_edge(CPH3& m, CPH3::CMAP::Edge e, bool set_indices)
{
	CPH3::CMAP& map = static_cast<CPH3::CMAP&>(m);

	CPH3::CMAP::Vertex v = cut_edge(map, e, false);

	Dart d = e.dart;
	do
	{
		m.set_edge_id(phi1(map, d), m.edge_id(d));
		m.set_edge_id(phi3(map, d), m.edge_id(d));
		m.set_edge_id(phi2(map, d), m.edge_id(phi<1, 2>(map, d)));
		m.set_face_id(phi1(map, d), m.face_id(d));
		m.set_face_id(phi3(map, d), m.face_id(d));
		m.set_face_id(phi2(map, d), m.face_id(phi<1, 2>(map, d)));
		m.set_dart_level(phi1(map, d), m.current_level_);
		m.set_dart_level(phi2(map, d), m.current_level_);
		d = phi<2, 3>(map, d);
	} while (d != e.dart);
	if (set_indices)
	{
		if (is_indexed<CPH3::CMAP::Vertex>(m))
			set_index(map, v, new_index<CPH3::CMAP::Vertex>(m));
		if (is_indexed<CPH3::CMAP::Edge>(m))
		{
			uint32 ne = new_index<CPH3::CMAP::Edge>(m);
			foreach_dart_of_orbit(map, e, [&](Dart d) -> bool {
				if (m.dart_level(d) == m.current_level_)
					set_index<CPH3::CMAP::Edge>(m, d, ne);
				return true;
			});
			ne = new_index<CPH3::CMAP::Edge>(m);
			foreach_dart_of_orbit(map, CPH3::CMAP::Edge(phi1(map, e.dart)), [&](Dart d) -> bool {
				if (m.dart_level(d) == m.current_level_)
					set_index<CPH3::CMAP::Edge>(m, d, ne);
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Face>(m))
		{
			d = e.dart;
			do
			{
				Dart it = phi1(map, d);
				do
				{
					it = phi1(map, it);
				} while (m.dart_level(it) < m.current_level_ - 1 && it != d);

				copy_index<CPH3::CMAP::Face>(map, phi1(map, d), it);
				it = phi2(map, d);
				do
				{
					it = phi1(map, it);
				} while (m.dart_level(it) < m.current_level_ - 1 && it != phi2(map, phi1(map, d)));
				copy_index<CPH3::CMAP::Face>(map, phi2(map, d), it);
				d = phi<2, 3>(map, d);
			} while (d != e.dart);
		}
		if (is_indexed<CPH3::CMAP::Volume>(m))
		{
			d = e.dart;
			do
			{
				if (!is_boundary(m, d))
				{
					Dart it = phi1(map, d);
					do
					{
						it = phi1(map, it);
					} while (m.dart_level(it) < m.current_level_ - 1 && it != d);
					copy_index<CPH3::CMAP::Volume>(map, phi1(map, d), it);
					it = phi2(map, d);
					do
					{
						it = phi1(map, it);
					} while (m.dart_level(it) < m.current_level_ - 1 && it != phi2(map, phi1(map, d)));
					copy_index<CPH3::CMAP::Volume>(map, phi2(map, d), it);
				}
				d = phi<2, 3>(map, d);
			} while (d != e.dart);
		}
	}

	return v;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// collapse_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////
void remove_edge_stability(IncidenceGraph& ig, IncidenceGraph::Edge e)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;
	
	auto [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];
	
	remove_edge(ig, e);

	//remove isolated vertex
	if ((*ig.vertex_incident_edges_)[v1.index_].size() == 0)
	{
		remove_cell<Vertex>(ig,v1);
	}
	if ((*ig.vertex_incident_edges_)[v2.index_].size() == 0)
	{
		remove_cell<Vertex>(ig, v2);
	}

 }


 
std::pair<IncidenceGraph::Vertex, std::vector<IncidenceGraph::Edge>> collapse_edge_qmat(IncidenceGraph& ig,
	IncidenceGraph::Edge e)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;
	
	Vertex new_vertex = add_vertex(ig);

	std::vector<Edge> removed_edges;
	auto [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];
	std::cout << "delete edge: " << e.index_ << std::endl;
	remove_edge(ig, e);

	foreach_incident_edge(ig, v1, [&](Edge iedge1) {
		std::cout << "v1: incident edge: " << iedge1.index_ << std::endl;
		return true;
		});
	foreach_incident_edge(ig, v2, [&](Edge iedge2) {
		std::cout << "v2: incident edge: " << iedge2.index_ << std::endl;
		return true;		
		});

	std::cout << "-----------------------------------------------" << std::endl;
	// replace v1 by new vertex in incident edges of v1
	for (Edge iev1 : (*ig.vertex_incident_edges_)[v1.index_])
	{
		// add old edge to removed edges
		removed_edges.push_back(iev1);
		foreach_incident_face(ig, iev1, [&](Face f) {
			std::cout << "edge " << iev1.index_ << ": incident face: " << f.index_ << ", ";
			std::cout << ::std::endl;
			return true;
		});
		// replace v1 by new vertex in incident edges of v1
		replace_vertex_in_edge(ig, iev1, v1, new_vertex);
		(*ig.vertex_incident_edges_)[new_vertex.index_].push_back(iev1);
	}

	// replace v2 by new vertex in incident edges of v2
	for (Edge iev2 : (*ig.vertex_incident_edges_)[v2.index_])
	{
		// add old edge to removed edges
		removed_edges.push_back(iev2);
		foreach_incident_face(ig, iev2, [&](Face f) {
			std::cout << "edge " << iev2.index_ << ": incident face: " << f.index_ << ", ";
			std::cout << ::std::endl;
			return true;
		});
		// replace v2 by new vertex in incident edges of v2
		replace_vertex_in_edge(ig, iev2, v2, new_vertex);
		(*ig.vertex_incident_edges_)[new_vertex.index_].push_back(iev2);
		
	}

	foreach_incident_edge(ig, new_vertex, [&](Edge iie) {
		std::cout << "new vertex:" << new_vertex.index_ << " incident edge: " << iie.index_
				  << "vertex1: " << (*ig.edge_incident_vertices_)[iie.index_].first.index_<< ", vertex 2: " <<(*ig.edge_incident_vertices_)[iie.index_].second.index_ << std::endl;
		return true;
	});

	// Check duplicate edges incident to new vertex
	std::unordered_map<uint32, Edge> adjacent_vertex_new_vertex;//adjacent vertex index, incident edge
	std::vector<Edge> duplicate_edges_after_delete;
	std::cout << "start duplicate edge check" << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
	foreach_incident_edge(ig, new_vertex, [&](Edge ie) {
		std::cout << "loop incident edge: " << ie.index_ << std::endl;
 		auto [va, vb] = (*ig.edge_incident_vertices_)[ie.index_];
 		Vertex other_vertex = (va == new_vertex) ? vb : va;
 
 		auto it = adjacent_vertex_new_vertex.find(other_vertex.index_);
 		if (it == adjacent_vertex_new_vertex.end())
 		{	
 			adjacent_vertex_new_vertex.insert({other_vertex.index_, ie});
 			std::cout << "insert tuple: {" << ie.index_ << ", " << other_vertex.index_ << "}" << std::endl;
 		}
		else
		{
			// replace current edge by old edge
			foreach_incident_face(ig, ie, [&](Face f) {
				std::cout << "face: " << f.index_ << ", replace edge " << ie.index_ << " by edge " << it->second.index_
						  << std::endl;
				replace_edge_in_face(ig, f, ie, it->second);
				(*ig.edge_incident_faces_)[it->second.index_].push_back(f);
				return true;
			});
 			// remove duplicate edge
			remove_edge_in_vertex(ig, other_vertex, ie);
			duplicate_edges_after_delete.push_back(ie);
 			std::cout << "remove edge: " << ie.index_ << " in vertex " << other_vertex.index_ << std::endl;	
 			remove_cell<Edge>(ig, ie);
 		}
		std::cout << "-------------------------------------------------" << std::endl;
		return true;
	});
	for (Edge ie : duplicate_edges_after_delete)
	{
		remove_edge_in_vertex(ig, new_vertex, ie);
	}
	std::cout << "\n" << std::endl;
	remove_cell<Vertex>(ig, v1);
	remove_cell<Vertex>(ig, v2);
	return {new_vertex, removed_edges};
}

std::pair<IncidenceGraph::Vertex, std::vector<IncidenceGraph::Edge>> collapse_edge(IncidenceGraph& ig,
																				   IncidenceGraph::Edge e,
																				   bool /*set_indices*/)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	std::vector<Edge> removed_edges;
	auto [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];

	// remove e from its incident vertices
	remove_edge_in_vertex(ig, v1, e);
	remove_edge_in_vertex(ig, v2, e);

	// remove e from its incident faces
	for (Face iface : (*ig.edge_incident_faces_)[e.index_])
	{
		remove_edge_in_face(ig, iface, e);
		// remove degenerate faces
		if ((*ig.face_incident_edges_)[iface.index_].size() < 3)
			remove_face(ig, iface);
	}

	// replace v1 by v2 in incident edges of v1
	for (Edge iev1 : (*ig.vertex_incident_edges_)[v1.index_])
	{
		replace_vertex_in_edge(ig, iev1, v1, v2);
		// check for duplicate edges around v2
		Edge similar_edge_in_v2;
		for (uint32 i = 0; !similar_edge_in_v2.is_valid() && i < (*ig.vertex_incident_edges_)[v2.index_].size(); ++i)
		{
			Edge iev2 = (*ig.vertex_incident_edges_)[v2.index_][i];
			if (same_edge(ig, iev1, iev2))
				similar_edge_in_v2 = iev2;
		}
		//Didn't find a duplicate edge in v2
		if (!similar_edge_in_v2.is_valid())
			(*ig.vertex_incident_edges_)[v2.index_].push_back(iev1);
		else
		{
			// migrate faces of iev1 to the similar edge in v2
			for (Face iface : (*ig.edge_incident_faces_)[iev1.index_])
			{
				auto fit = std::find((*ig.edge_incident_faces_)[similar_edge_in_v2.index_].begin(),
									 (*ig.edge_incident_faces_)[similar_edge_in_v2.index_].end(), iface);
				if (fit == (*ig.edge_incident_faces_)[similar_edge_in_v2.index_].end())
				{
					replace_edge_in_face(ig, iface, iev1, similar_edge_in_v2);
					(*ig.edge_incident_faces_)[similar_edge_in_v2.index_].push_back(iface);
				}
			}
			// remove iev1 from its vertices
			auto [iev1v1, iev1v2] = (*ig.edge_incident_vertices_)[iev1.index_];
			remove_edge_in_vertex(ig, iev1v1, iev1);
			remove_edge_in_vertex(ig, iev1v2, iev1);
			// remove iev1
			remove_cell<Edge>(ig, iev1);
			removed_edges.push_back(iev1);
		}
	}

	// remove v1
	remove_cell<Vertex>(ig, v1);
	// remove e
	remove_cell<Edge>(ig, e);

	return {v2, removed_edges};
}

std::vector<IncidenceGraph::Edge> collapse_edge_with_fixed_vertices(
	IncidenceGraph& ig, IncidenceGraph::Edge e, std::vector<double> selected_points)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;
	std::vector<Edge> removed_edges;

	auto [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];
	if (selected_points[v1.index_] > 0 && selected_points[v2.index_] > 0)
	{
	}
	else
	{
		// if v1 is selected, do nothing
		// if v2 is selected, swap v1 and v2
		if (selected_points[index_of(ig, v1)])
		{
			Vertex tmp = v1;
			v1 = v2;
			v2 = tmp;
		}
		// if none is selected, do nothing
		// remove e from its incident vertices
		remove_edge_in_vertex(ig, v1, e);
		remove_edge_in_vertex(ig, v2, e);

		// remove e from its incident faces
		for (Face iface : (*ig.edge_incident_faces_)[e.index_])
		{
			remove_edge_in_face(ig, iface, e);
			// remove degenerate faces
			if ((*ig.face_incident_edges_)[iface.index_].size() < 3)
				remove_face(ig, iface);
		}

		// replace v1 by v2 in incident edges of v1
		for (Edge iev1 : (*ig.vertex_incident_edges_)[v1.index_])
		{
			replace_vertex_in_edge(ig, iev1, v1, v2);
			// check for duplicate edges around v2
			Edge similar_edge_in_v2;
			for (uint32 i = 0; !similar_edge_in_v2.is_valid() && i < (*ig.vertex_incident_edges_)[v2.index_].size();
				 ++i)
			{
				Edge iev2 = (*ig.vertex_incident_edges_)[v2.index_][i];
				if (same_edge(ig, iev1, iev2))
					similar_edge_in_v2 = iev2;
			}
			if (!similar_edge_in_v2.is_valid())
				(*ig.vertex_incident_edges_)[v2.index_].push_back(iev1);
			else
			{
				// migrate faces of iev1 to the similar edge in v2
				for (Face iface : (*ig.edge_incident_faces_)[iev1.index_])
				{
					auto fit = std::find((*ig.edge_incident_faces_)[similar_edge_in_v2.index_].begin(),
										 (*ig.edge_incident_faces_)[similar_edge_in_v2.index_].end(), iface);
					if (fit == (*ig.edge_incident_faces_)[similar_edge_in_v2.index_].end())
					{
						replace_edge_in_face(ig, iface, iev1, similar_edge_in_v2);
						(*ig.edge_incident_faces_)[similar_edge_in_v2.index_].push_back(iface);
					}
				}
				// remove iev1 from its vertices
				auto [iev1v1, iev1v2] = (*ig.edge_incident_vertices_)[iev1.index_];
				remove_edge_in_vertex(ig, iev1v1, iev1);
				remove_edge_in_vertex(ig, iev1v2, iev1);
				// remove iev1
				remove_cell<Edge>(ig, iev1);
				removed_edges.push_back(iev1);
			}
		}

		// remove v1
		remove_cell<Vertex>(ig, v1);
		// remove e
		remove_cell<Edge>(ig, e);
	}
	return removed_edges;
// 		collapse edge to v1
// 				std::unordered_set<Vertex*> adjacent_vertices_to_v2;
// 				// remove adjacent edges of v2 and store adjacent vertices of v2 which will serve to create new edges
// 				std::vector<Edge> adjacent_edges_to_v2 = (*ig.vertex_incident_edges_)[v2.index_];
// 				for (Edge& adjacent_edge : (*ig.vertex_incident_edges_)[v2.index_])
// 				{
// 					auto ev2 = (*ig.edge_incident_vertices_)[adjacent_edge.index_];
// 					if (ev2.first == v2 && ev2.second !=v1) 
// 					{
// 						adjacent_vertices_to_v2.insert(&ev2.second);
// 					}
// 					if (ev2.second == v2 && ev2.first != v1)
// 					{
// 						adjacent_vertices_to_v2.insert(&ev2.first);
// 					}
// 					// remove incident face of edge
// 					
// 					if (ev2.first != v1 && ev2.second != v1)
// 					{
// 						for (Face iface : (*ig.edge_incident_faces_)[adjacent_edge.index_])
// 						{
// 							remove_edge_in_face(ig, iface, adjacent_edge);
// 							// remove degenerate faces
// 							if ((*ig.face_incident_edges_)[iface.index_].size() < 3)
// 								remove_face(ig, iface);
// 						}
// 						// remove edge from vertices
// 						remove_edge_in_vertex(ig, ev2.first, adjacent_edge);
// 						remove_edge_in_vertex(ig, ev2.second, adjacent_edge);
// 						// remove edge from graph
// 						remove_cell<Edge>(ig, adjacent_edge);
// 						removed_edges.push_back(adjacent_edge);
// 					}
// 				}
// 				// test all adjacent vertices of v2 to see if they are already adjacent to v1, if not, create new edge
// 				for (Edge& adjacent_edge : (*ig.vertex_incident_edges_)[v1.index_])
// 				{
// 					auto ev1 = (*ig.edge_incident_vertices_)[adjacent_edge.index_];
// 					if (adjacent_vertices_to_v2.find(&ev1.first) != adjacent_vertices_to_v2.end())
// 					{
// 						adjacent_vertices_to_v2.erase(&ev1.first);
// 					}
// 					if (adjacent_vertices_to_v2.find(&ev1.second) != adjacent_vertices_to_v2.end())
// 					{
// 						adjacent_vertices_to_v2.erase(&ev1.second);
// 					}
// 				}
// 				// create new edges
// 				for (Vertex* va : adjacent_vertices_to_v2)
// 				{
// 					added_edges.push_back(add_edge(ig, v1, *va));
// 				}
// 				
// 				remove_edge(ig, e);
// 				// remove vertex
// 				remove_cell<Vertex>(ig, v2);
// 	}
// 	return std::make_pair(removed_edges, added_edges);
}
///////////
// Graph //
///////////

Graph::Vertex collapse_edge(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart d = e.dart;
	Dart d1 = alpha_1(g, d);
	Dart dd = alpha0(g, d);
	Dart dd1 = alpha_1(g, dd);

	Graph::Vertex v;

	if (dd1 != dd)
	{
		v.dart = dd1;
		alpha1_unsew(g, dd);
	}
	if (d1 != d)
	{
		v.dart = d1;
		alpha1_unsew(g, d);
	}
	if (d1 != d && dd1 != dd)
		alpha1_sew(g, d1, dd1);
	remove_dart(g, d);
	remove_dart(g, dd);

	if (set_indices && v.is_valid())
	{
		if (is_indexed<Graph::Vertex>(g))
			set_index(g, v, index_of(g, v));
	}

	return v;
}

///////////
// CMap1 //
///////////

CMap1::Vertex collapse_edge(CMap1& m, CMap1::Edge e, bool set_indices)
{
	Dart d = phi_1(m, e.dart);
	phi1_unsew(m, d);
	remove_dart(m, e.dart);
	CMap1::Vertex v(phi1(m, d));

	if (set_indices)
	{
	}

	return v;
}

///////////
// CMap2 //
///////////

CMap2::Vertex collapse_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart dd = e.dart;
	Dart dd_1 = phi_1(m, dd);
	Dart dd_12 = phi2(m, dd_1);
	Dart ee = phi2(m, dd);
	Dart ee_1 = phi_1(m, ee);
	Dart ee_12 = phi2(m, ee_1);

	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);

	if (codegree(m, CMap2::Face(dd_1)) == 2u)
	{
		Dart dd1 = phi1(m, dd_1);
		Dart dd12 = phi2(m, dd1);
		phi2_unsew(m, dd1);
		phi2_unsew(m, dd_1);
		phi2_sew(m, dd12, dd_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(dd1), false);
	}

	if (codegree(m, CMap2::Face(ee_1)) == 2u)
	{
		Dart ee1 = phi1(m, ee_1);
		Dart ee12 = phi2(m, ee1);
		phi2_unsew(m, ee1);
		phi2_unsew(m, ee_1);
		phi2_sew(m, ee12, ee_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(ee1), false);
	}

	CMap2::Vertex v(dd_12);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
			set_index(m, v, index_of(m, v));
		if (is_indexed<CMap2::Edge>(m))
		{
			copy_index<CMap2::Edge>(m, dd_12, phi2(m, dd_12));
			copy_index<CMap2::Edge>(m, ee_12, phi2(m, ee_12));
		}
	}

	return v;
}

/*****************************************************************************/

// template <typename MESH>
// bool
// flip_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

bool flip_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart d = e.dart;
	Dart dd = phi2(m, d);
	Dart d1 = phi1(m, d);
	Dart d_1 = phi_1(m, d);
	Dart dd1 = phi1(m, dd);
	Dart dd_1 = phi_1(m, dd);

	// // Cannot flip edge whose incident faces have co-degree 1
	// if (d == d1 || dd == dd1)
	// 	return false;

	// // Both vertices have degree 1 and thus nothing is done // TODO may return true ?
	// if (d == dd_1 && dd == d_1)
	// 	return false;

	// if (d != dd_1)
	// 	phi1_sew(m, d, dd_1); // Detach the edge from its
	// if (dd != d_1)
	// 	phi1_sew(m, dd, d_1); // two incident vertices

	// if (d != dd_1)
	// 	phi1_sew(m, d, d1); // Insert the first end in its new vertices
	// if (dd != d_1)
	// 	phi1_sew(m, dd, dd1); // Insert the second end in its new vertices

	phi1_sew(m, d, dd_1);
	phi1_sew(m, dd, d_1);
	phi1_sew(m, d, d1);
	phi1_sew(m, dd, dd1);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			copy_index<CMap2::Vertex>(m, d, phi1(m, dd));
			copy_index<CMap2::Vertex>(m, dd, phi1(m, d));
		}

		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m, phi_1(m, d), d);
			copy_index<CMap2::Face>(m, phi_1(m, dd), dd);
		}
	}

	return true;
}

} // namespace cgogn
