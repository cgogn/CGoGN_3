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

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>

#include <cgogn/core/types/cell_marker.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

IncidenceGraph::Vertex add_vertex(IncidenceGraph& ig)
{
	using Vertex = IncidenceGraph::Vertex;

	Vertex v = add_cell<Vertex>(ig);
	(*ig.vertex_incident_edges_)[v.index_].clear();
	return v;
}

IncidenceGraph::Edge connect_vertices(IncidenceGraph& ig, IncidenceGraph::Vertex v1, IncidenceGraph::Vertex v2)
{
	using Edge = IncidenceGraph::Edge;

	Edge e = add_cell<Edge>(ig);
	(*ig.edge_incident_vertices_)[e.index_] = {v1, v2};
	(*ig.edge_incident_faces_)[e.index_].clear();
	(*ig.vertex_incident_edges_)[v1.index_].push_back(e);
	(*ig.vertex_incident_edges_)[v2.index_].push_back(e);

	return e;
}

void remove_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v)
{
	using Vertex = IncidenceGraph::Vertex;

	while ((*ig.vertex_incident_edges_)[v.index_].size() > 0)
		remove_edge(ig, (*ig.vertex_incident_edges_)[v.index_].back());
	remove_cell<Vertex>(ig, v);
}

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

IncidenceGraph::Vertex cut_edge(IncidenceGraph& ig, IncidenceGraph::Edge e0, bool /*set_indices*/)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	auto [v0, v1] = (*ig.edge_incident_vertices_)[e0.index_];
	Vertex v = add_cell<Vertex>(ig);
	(*ig.edge_incident_vertices_)[e0.index_] = {v0, v};
	Edge e1 = add_edge(ig, v, v1);
	(*ig.vertex_incident_edges_)[v.index_].push_back(e0);
	for (Face f : (*ig.edge_incident_faces_)[e0.index_])
	{
		(*ig.edge_incident_faces_)[e1.index_].push_back(f);
		(*ig.face_incident_edges_)[f.index_].push_back(e1);
		sort_face_edges(ig, f); // TODO: could do more efficient (insert)
	}
	return v;
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

std::vector<IncidenceGraph::Edge> collapse_edge_qmat(IncidenceGraph& ig, IncidenceGraph::Edge e,
													 IncidenceGraph::Vertex& new_vertex)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	std::vector<Edge> removed_edges;
	auto [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];
	remove_edge(ig, e);

	// replace v1 by new vertex in incident edges of v1
	for (Edge iev1 : (*ig.vertex_incident_edges_)[v1.index_])
	{
		// add old edge to removed edges
		removed_edges.push_back(iev1);
		// replace v1 by new vertex in incident edges of v1
		replace_vertex_in_edge(ig, iev1, v1, new_vertex);
		(*ig.vertex_incident_edges_)[new_vertex.index_].push_back(iev1);
	}

	// replace v2 by new vertex in incident edges of v2
	for (Edge iev2 : (*ig.vertex_incident_edges_)[v2.index_])
	{
		// add old edge to removed edges
		removed_edges.push_back(iev2);
		// replace v2 by new vertex in incident edges of v2
		replace_vertex_in_edge(ig, iev2, v2, new_vertex);
		(*ig.vertex_incident_edges_)[new_vertex.index_].push_back(iev2);
	}

	// Check duplicate edges incident to new vertex
	std::unordered_map<uint32, Edge> adjacent_vertex_new_vertex; // adjacent vertex index, incident edge
	std::vector<Edge> duplicate_edges_after_delete;

	foreach_incident_edge(ig, new_vertex, [&](Edge ie) {
		auto [va, vb] = (*ig.edge_incident_vertices_)[ie.index_];
		Vertex other_vertex = (va == new_vertex) ? vb : va;

		auto it = adjacent_vertex_new_vertex.find(other_vertex.index_);
		if (it == adjacent_vertex_new_vertex.end())
		{
			adjacent_vertex_new_vertex.insert({other_vertex.index_, ie});
		}
		else
		{
			// replace current edge by old edge
			foreach_incident_face(ig, ie, [&](Face f) {
				replace_edge_in_face(ig, f, ie, it->second);
				sort_face_edges(ig, f);
				(*ig.edge_incident_faces_)[it->second.index_].push_back(f);
				return true;
			});
			// remove duplicate edge
			remove_edge_in_vertex(ig, other_vertex, ie);
			duplicate_edges_after_delete.push_back(ie);
			remove_cell<Edge>(ig, ie);
		}
		return true;
	});
	for (Edge ie : duplicate_edges_after_delete)
	{
		remove_edge_in_vertex(ig, new_vertex, ie);
	}
	remove_cell<Vertex>(ig, v1);
	remove_cell<Vertex>(ig, v2);
	return removed_edges;
}

std::pair<IncidenceGraph::Vertex, std::vector<IncidenceGraph::Edge>> collapse_edge_qmat_with_debug_info(
	IncidenceGraph& ig, IncidenceGraph::Edge e)
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
		// add old edge into sturcture to perform delay deletion
		removed_edges.push_back(iev1);
		foreach_incident_face(ig, iev1, [&](Face f) {
			std::cout << "edge " << iev1.index_ << ": incident face: " << f.index_ << ", "
					  << "incident vertex: " << (*ig.edge_incident_vertices_)[iev1.index_].first.index_ << ", "
					  << (*ig.edge_incident_vertices_)[iev1.index_].second.index_;
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
		// add old edge into sturcture to perform delay deletion
		removed_edges.push_back(iev2);
		foreach_incident_face(ig, iev2, [&](Face f) {
			std::cout << "edge " << iev2.index_ << ": incident face: " << f.index_ << ", "
					  << "incident vertex: " << (*ig.edge_incident_vertices_)[iev2.index_].first.index_ << ", "
					  << (*ig.edge_incident_vertices_)[iev2.index_].second.index_;
			std::cout << ::std::endl;
			return true;
		});
		// replace v2 by new vertex in incident edges of v2
		replace_vertex_in_edge(ig, iev2, v2, new_vertex);
		(*ig.vertex_incident_edges_)[new_vertex.index_].push_back(iev2);
	}

	foreach_incident_edge(ig, new_vertex, [&](Edge iie) {
		std::cout << "new vertex:" << new_vertex.index_ << " incident edge: " << iie.index_
				  << " vertex1: " << (*ig.edge_incident_vertices_)[iie.index_].first.index_
				  << ", vertex 2: " << (*ig.edge_incident_vertices_)[iie.index_].second.index_ << std::endl;
		return true;
	});

	// Check duplicate edges incident to new vertex
	std::unordered_map<uint32, Edge> adjacent_vertex_new_vertex; // adjacent vertex index, incident edge
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
				sort_face_edges(ig, f);
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
	for (Edge ie : (*ig.vertex_incident_edges_)[new_vertex.index_])
	{
		auto [iev1, iev2] = (*ig.edge_incident_vertices_)[ie.index_];
		std::cout << "edge " << ie.index_ << " vertex1 : " << iev1.index_ << ", vertex2 : " << iev2.index_
				  << " is incident to: \n";
		for (Face iface : (*ig.edge_incident_faces_)[ie.index_])
		{
			std::cout << "Face : " << iface.index_ << std::endl;
			std::cout << "Face: " << iface.index_ << " is incident to: \n";
			for (Edge iie : (*ig.face_incident_edges_)[iface.index_])
			{
				std::cout << "edge : " << iie.index_ << std::endl;
				std::cout << "edge : " << iie.index_ << " is incident to:\n ";
				auto iiv = (*ig.edge_incident_vertices_)[iie.index_];
				{
					std::cout << iiv.first.index_ << ", " << iiv.second.index_ << std::endl;
				}
			}
			std::cout << "Face: " << iface.index_ << " is incident to vertex: ";
			foreach_incident_vertex(ig, iface, [&](Vertex iiiv) {
				std::cout << iiiv.index_ << ", ";
				return true;
			});
			std::cout << std::endl;
		}
		std::cout << "--------------------------------------------------" << std::endl;
	}
	std::cout << "\n" << std::endl;
	remove_cell<Vertex>(ig, v1);
	remove_cell<Vertex>(ig, v2);
	return {new_vertex, removed_edges};
}

std::vector<IncidenceGraph::Edge> collapse_edge_with_fixed_vertices(IncidenceGraph& ig, IncidenceGraph::Edge e,
																	std::vector<double> selected_points)
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
}

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

IncidenceGraph::Face add_face(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge>& edges)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	Face f = add_cell<Face>(ig);
	(*ig.face_incident_edges_)[f.index_] = edges;
	if (sort_face_edges(ig, f))
	{
		for (Edge e : edges)
			(*ig.edge_incident_faces_)[e.index_].push_back(f);
		return f;
	}
	else
	{
		remove_cell<Face>(ig, f);
		return Face();
	}
}

void remove_face(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	for (Edge e : (*ig.face_incident_edges_)[f.index_])
		remove_face_in_edge(ig, e, f);
	remove_cell<Face>(ig, f);
}

IncidenceGraph::Edge cut_face(IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	// TODO: manage face_incident_edges_dir_ !!!

	// find common face
	std::vector<Face> faces0 = incident_faces(ig, v0);
	std::vector<Face> faces1 = incident_faces(ig, v1);

	Face face;
	for (uint32 i = 0; i < faces0.size(); ++i)
	{
		for (uint32 j = 0; j < faces1.size(); ++j)
		{
			if (faces0[i] == faces1[j])
			{
				face = faces0[i];
				break;
			}
		}
		if (face.is_valid())
			break;
	}

	if (!face.is_valid())
		return Edge();

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[face.index_];
	std::vector<Vertex> vertices = sorted_face_vertices(ig, face);

	std::vector<Edge> face_edge0;
	std::vector<Edge> face_edge1;

	bool inside = false;
	for (uint32 i = 0; i < edges.size(); ++i)
	{
		if (vertices[i] == v0 || vertices[i] == v1)
			inside = !inside;

		if (inside)
			face_edge1.push_back(edges[i]);
		else
			face_edge0.push_back(edges[i]);
	}

	remove_face(ig, face);
	Edge new_edge = add_edge(ig, v0, v1);
	face_edge0.push_back(new_edge);
	face_edge1.push_back(new_edge);
	add_face(ig, face_edge0);
	add_face(ig, face_edge1);

	return new_edge;
}

/*************************************************************************/
// Helper functions
/*************************************************************************/

bool same_edge(IncidenceGraph& ig, IncidenceGraph::Edge e1, IncidenceGraph::Edge e2)
{
	using Vertex = IncidenceGraph::Vertex;

	Vertex e1v1 = (*ig.edge_incident_vertices_)[e1.index_].first;
	Vertex e1v2 = (*ig.edge_incident_vertices_)[e1.index_].second;
	Vertex e2v1 = (*ig.edge_incident_vertices_)[e2.index_].first;
	Vertex e2v2 = (*ig.edge_incident_vertices_)[e2.index_].second;

	return (e1v1 == e2v1 && e1v2 == e2v2) || (e1v1 == e2v2 && e1v2 == e2v1);
}

IncidenceGraph::Vertex common_vertex(IncidenceGraph& ig, IncidenceGraph::Edge e0, IncidenceGraph::Edge e1)
{
	using Vertex = IncidenceGraph::Vertex;

	Vertex v0 = (*ig.edge_incident_vertices_)[e0.index_].first;
	Vertex v1 = (*ig.edge_incident_vertices_)[e0.index_].second;
	Vertex v2 = (*ig.edge_incident_vertices_)[e1.index_].first;
	Vertex v3 = (*ig.edge_incident_vertices_)[e1.index_].second;

	if (v0 == v2 || v0 == v3)
		return v0;

	if (v1 == v2 || v1 == v3)
		return v1;

	return Vertex();
}

void remove_edge_in_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v, IncidenceGraph::Edge edge_to_remove)
{
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.vertex_incident_edges_)[v.index_];
	auto eit = std::find(edges.begin(), edges.end(), edge_to_remove);
	if (eit != edges.end())
		edges.erase(eit);
}

void remove_face_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Face face_to_remove)
{
	using Face = IncidenceGraph::Face;

	std::vector<Face>& faces = (*ig.edge_incident_faces_)[e.index_];
	auto fit = std::find(faces.begin(), faces.end(), face_to_remove);
	if (fit != faces.end())
		faces.erase(fit);
}

void replace_vertex_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Vertex old_vertex,
							IncidenceGraph::Vertex new_vertex)
{
	if ((*ig.edge_incident_vertices_)[e.index_].first == old_vertex)
		(*ig.edge_incident_vertices_)[e.index_].first = new_vertex;
	if ((*ig.edge_incident_vertices_)[e.index_].second == old_vertex)
		(*ig.edge_incident_vertices_)[e.index_].second = new_vertex;
}

void remove_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge edge_to_remove)
{
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<uint8>& edges_dir = (*ig.face_incident_edges_dir_)[f.index_];
	auto eit = edges.begin();
	auto edit = edges_dir.begin();
	for (Edge e : edges)
	{
		if (e == edge_to_remove)
			break;
		eit++;
		edit++;
	}
	if (eit != edges.end())
	{
		edges.erase(eit);
		edges_dir.erase(edit);
	}
}

void replace_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge old_edge,
						  IncidenceGraph::Edge new_edge)
{
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	auto eit = std::find(edges.begin(), edges.end(), old_edge);
	if (eit != edges.end())
		*eit = new_edge;
}

std::vector<IncidenceGraph::Vertex> sorted_face_vertices(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<Vertex> sorted_vertices;
	for (uint32 i = 0; i < edges.size(); ++i)
		sorted_vertices.push_back(common_vertex(ig, edges[i], edges[(i + 1) % edges.size()]));
	sorted_vertices.insert(sorted_vertices.begin(), sorted_vertices.back());
	sorted_vertices.pop_back();
	return sorted_vertices;
}

bool sort_face_edges(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	std::vector<Edge>& edges = (*ig.face_incident_edges_)[f.index_];
	std::vector<uint8>& edges_dir = (*ig.face_incident_edges_dir_)[f.index_];
	edges_dir.clear();

	std::vector<Edge> unordered_edges;
	unordered_edges.swap(edges);

	edges.push_back(unordered_edges.front());
	edges_dir.push_back(0);

	unordered_edges.erase(unordered_edges.begin());
	uint32 vid0 = (*ig.edge_incident_vertices_)[edges.front().index_].first.index_;
	uint32 vid1 = (*ig.edge_incident_vertices_)[edges.front().index_].second.index_;
	bool broken = false;

	while (unordered_edges.size() > 0)
	{
		std::size_t i, end;
		for (i = 0, end = unordered_edges.size(); i < unordered_edges.size(); ++i)
		{
			Edge e = unordered_edges[i];
			std::pair<Vertex, Vertex>& evs = (*ig.edge_incident_vertices_)[e.index_];
			uint32 ev = (evs.first.index_ == vid1 ? evs.second.index_
												  : (evs.second.index_ == vid1 ? evs.first.index_ : INVALID_INDEX));

			if (ev != INVALID_INDEX)
			{
				vid1 = ev;
				edges.push_back(e);
				edges_dir.push_back(ev == evs.first.index_ ? 1 : 0);
				unordered_edges.erase(unordered_edges.begin() + i);
				break;
			}
		}

		broken = (vid1 == vid0) || (i == end);
		if (broken)
			break;
	}

	return (broken && unordered_edges.size() == 0);
}

} // namespace cgogn
