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
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>
#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// add_face(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

CMap1::Face add_face(CMap1& m, uint32 size, bool set_indices)
{
	Dart d = add_dart(m);
	for (uint32 i = 1u; i < size; ++i)
	{
		Dart e = add_dart(m);
		phi1_sew(m, d, e);
	}
	CMap1::Face f(d);

	if (set_indices)
	{
		if (is_indexed<CMap1::Vertex>(m))
		{
			foreach_incident_vertex(
				m, f,
				[&](CMap1::Vertex v) -> bool {
					set_index(m, v, new_index<CMap1::Vertex>(m));
					return true;
				},
				CMapBase::TraversalPolicy::DART_MARKING);
		}
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (is_indexed<CMap1::Face>(m))
			set_index(m, f, new_index<CMap1::Face>(m));
	}

	return f;
}

///////////
// CMap2 //
///////////

CMap2::Face add_face(CMap2& m, uint32 size, bool set_indices)
{
	CMap2::Face f = add_face(static_cast<CMap1&>(m), size, false);
	CMap2::Face b = add_face(static_cast<CMap1&>(m), size, false);
	Dart it = b.dart;
	foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
		set_boundary(m, it, true);
		phi2_sew(m, d, it);
		it = phi_1(m, it);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(
				m, f,
				[&](CMap2::Vertex v) -> bool {
					set_index(m, v, new_index<CMap2::Vertex>(m));
					return true;
				},
				CMapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(
				m, f,
				[&](CMap2::Edge e) -> bool {
					set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
					return true;
				},
				CMapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(
				m, f,
				[&](CMap2::Edge e) -> bool {
					set_index(m, e, new_index<CMap2::Edge>(m));
					return true;
				},
				CMapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Face>(m))
			set_index(m, f, new_index<CMap2::Face>(m));
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, CMap2::Volume(f.dart), new_index<CMap2::Volume>(m));
	}

	return f;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// add_face(MESH& m, std::vector<typename mesh_traits<MESH>::Edge edges);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

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

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// remove_face(MESH& m, typename mesh_traits<MESH>::Face f, bool set_indices = true);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

void remove_face(IncidenceGraph& ig, IncidenceGraph::Face f)
{
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;
	for (Edge e : (*ig.face_incident_edges_)[f.index_])
	{
		remove_face_in_edge(ig, e, f);
	}
	remove_cell<Face>(ig, f);
}

///////////
// CMap1 //
///////////

void remove_face(CMap1& m, CMap1::Face f, bool set_indices)
{
	Dart it = phi1(m, f.dart);
	while (it != f.dart)
	{
		Dart next = phi1(m, it);
		remove_dart(m, it);
		it = next;
	}
	remove_dart(m, f.dart);

	if (set_indices)
	{
	}
}

/*****************************************************************************/

// template <typename MESH>
// void
// merge_incident_faces(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

void CGOGN_CORE_EXPORT merge_incident_faces(CMap2& m, CMap2::Edge e, bool set_indices)
{
	if (is_incident_to_boundary(m, e))
		return;

	Dart d = e.dart;
	Dart d_1 = phi_1(m, d);
	Dart d2 = phi2(m, d);
	Dart d2_1 = phi_1(m, d2);

	phi1_sew(m, d_1, d2);
	phi1_sew(m, d2_1, d);

	remove_face(static_cast<CMap1&>(m), CMap1::Face(d), set_indices);

	if (set_indices)
	{
		if (is_indexed<CMap2::Face>(m))
			set_index(m, CMap2::Face(d_1), index_of(m, CMap2::Face(d_1)));
	}
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// cut_face(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_indices =
// true);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

IncidenceGraph::Edge CGOGN_CORE_EXPORT cut_face(IncidenceGraph& ig, IncidenceGraph::Vertex v0,
												IncidenceGraph::Vertex v1)
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

///////////
// CMap2 //
///////////

CMap2::Edge cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices)
{
	Dart dd = phi_1(m, v1.dart);
	Dart ee = phi_1(m, v2.dart);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);
	phi1_sew(m, nv1.dart, nv2.dart);
	phi2_sew(m, nv1.dart, nv2.dart);
	set_boundary(m, nv1.dart, is_boundary(m, dd));
	set_boundary(m, nv2.dart, is_boundary(m, ee));
	CMap2::Edge e(nv1.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			copy_index<CMap2::Vertex>(m, nv1.dart, v1.dart);
			copy_index<CMap2::Vertex>(m, nv2.dart, v2.dart);
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			set_index(m, CMap2::HalfEdge(nv1.dart), new_index<CMap2::HalfEdge>(m));
			set_index(m, CMap2::HalfEdge(nv2.dart), new_index<CMap2::HalfEdge>(m));
		}
		if (is_indexed<CMap2::Edge>(m))
			set_index(m, CMap2::Edge(nv1.dart), new_index<CMap2::Edge>(m));
		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m, nv2.dart, v1.dart);
			set_index(m, CMap2::Face(v2.dart), new_index<CMap2::Face>(m));
		}
		if (is_indexed<CMap2::Volume>(m))
		{
			copy_index<CMap2::Volume>(m, nv1.dart, v2.dart);
			copy_index<CMap2::Volume>(m, nv2.dart, v1.dart);
		}
	}

	return e;
}

///////////
// CMap3 //
///////////

CMap3::Edge cut_face(CMap3& m, CMap3::Vertex v1, CMap3::Vertex v2, bool set_indices)
{
	Dart d = v1.dart;
	Dart e = v2.dart;

	Dart dd = phi<3, 1>(m, d);
	Dart ee = phi<3, 1>(m, e);

	cut_face(static_cast<CMap2&>(m), CMap2::Vertex(v1.dart), CMap2::Vertex(e), false);
	cut_face(static_cast<CMap2&>(m), CMap2::Vertex(dd), CMap2::Vertex(ee), false);

	phi3_sew(m, phi_1(m, v1.dart), phi_1(m, ee));
	phi3_sew(m, phi_1(m, dd), phi_1(m, e));

	CMap3::Edge edge(phi_1(m, e));

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
		{
			copy_index<CMap3::Vertex>(m, phi_1(m, e), v1.dart);
			copy_index<CMap3::Vertex>(m, phi_1(m, ee), v1.dart);
			copy_index<CMap3::Vertex>(m, phi_1(m, d), e);
			copy_index<CMap3::Vertex>(m, phi_1(m, dd), e);
		}
		if (is_indexed<CMap3::Vertex2>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CMap3::Vertex2>(m, phi_1(m, d), e);
				copy_index<CMap3::Vertex2>(m, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CMap3::Vertex2>(m, phi_1(m, dd), ee);
				copy_index<CMap3::Vertex2>(m, phi_1(m, ee), dd);
			}
		}
		if (is_indexed<CMap3::Edge>(m))
			set_index(m, CMap3::Edge(phi_1(m, v1.dart)), new_index<CMap3::Edge>(m));
		if (is_indexed<CMap3::Face2>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CMap3::Face2>(m, phi_1(m, d), d);
				set_index(m, CMap3::Face2(e), new_index<CMap3::Face2>(m));
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CMap3::Face2>(m, phi_1(m, dd), dd);
				set_index(m, CMap3::Face2(ee), new_index<CMap3::Face2>(m));
			}
		}
		if (is_indexed<CMap3::Face>(m))
		{
			copy_index<CMap3::Face>(m, phi_1(m, ee), d);
			copy_index<CMap3::Face>(m, phi_1(m, d), d);
			set_index(m, CMap3::Face(e), new_index<CMap3::Face>(m));
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CMap3::Volume>(m, phi_1(m, d), d);
				copy_index<CMap3::Volume>(m, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CMap3::Volume>(m, phi_1(m, dd), dd);
				copy_index<CMap3::Volume>(m, phi_1(m, ee), dd);
			}
		}
	}

	return edge;
}

//////////
// CPH3 //
//////////

CPH3::CMAP::Edge cut_face(CPH3& m, CPH3::CMAP::Vertex v1, CPH3::CMAP::Vertex v2, bool set_indices)
{
	CPH3::CMAP& map = static_cast<CPH3::CMAP&>(m);

	Dart d = v1.dart;
	Dart e = v2.dart;

	Dart dd = phi<3, 1>(m, v1.dart);
	Dart ee = phi<3, 1>(m, e);

	CPH3::CMAP::Edge result = cut_face(map, v1, v2, false);

	uint32 eid = m.refinement_edge_id(v1.dart, v2.dart);

	foreach_dart_of_orbit(m, result, [&](Dart d) -> bool {
		m.set_edge_id(d, eid);
		m.set_face_id(d, m.face_id(v1.dart));
		m.set_dart_level(d, m.current_level_);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CPH3::CMAP::Vertex>(m))
		{
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, e), v1.dart);
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, ee), v1.dart);
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, d), e);
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, dd), e);
		}
		if (is_indexed<CPH3::CMAP::Edge>(m))
			set_index(m, CPH3::CMAP::Edge(phi_1(m, v1.dart)), new_index<CPH3::CMAP::Edge>(m));
		if (is_indexed<CPH3::CMAP::Face>(m))
		{
			uint32 nf1 = new_index<CPH3::CMAP::Face>(m);
			uint32 nf2 = new_index<CPH3::CMAP::Face>(m);
			foreach_dart_of_orbit(m, CPH3::CMAP::Face(d), [&](Dart df) -> bool {
				if (m.current_level_ == m.dart_level(df))
					set_index<CPH3::CMAP::Face>(m, df, nf1);
				return true;
			});
			foreach_dart_of_orbit(m, CPH3::CMAP::Face(e), [&](Dart df) -> bool {
				if (m.current_level_ == m.dart_level(df))
					set_index<CPH3::CMAP::Face>(m, df, nf2);
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Volume>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, d), d);
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, dd), dd);
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, ee), dd);
			}
		}
	}

	return result;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// close_hole(MESH& m, Dart d, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Face close_hole(CMap2& m, Dart d, bool set_indices)
{
	cgogn_message_assert(phi2(m, d) == d, "CMap2: close hole called on a dart that is not a phi2 fix point");

	Dart first = add_dart(m); // First edge of the face that will fill the hole
	phi2_sew(m, d, first);	  // 2-sew the new edge to the hole

	Dart d_next = d; // Turn around the hole
	Dart d_phi1;	 // to complete the face
	do
	{
		do
		{
			d_phi1 = phi1(m, d_next); // Search and put in d_next
			d_next = phi2(m, d_phi1); // the next dart of the hole
		} while (d_next != d_phi1 && d_phi1 != d);

		if (d_phi1 != d)
		{
			Dart next = add_dart(m); // Add a vertex into the built face
			phi1_sew(m, first, next);
			phi2_sew(m, d_next, next); // and 2-sew the face to the hole
		}
	} while (d_phi1 != d);

	CMap2::Face hole(first);

	if (set_indices)
	{
		foreach_dart_of_orbit(m, hole, [&](Dart hd) -> bool {
			Dart hd2 = phi2(m, hd);
			if (is_indexed<CMap2::Vertex>(m))
				copy_index<CMap2::Vertex>(m, hd, phi1(m, hd2));
			if (is_indexed<CMap2::Edge>(m))
				copy_index<CMap2::Edge>(m, hd, hd2);
			if (is_indexed<CMap2::Volume>(m))
				copy_index<CMap2::Volume>(m, hd, hd2);
			return true;
		});
	}

	return hole;
}

/*****************************************************************************/

// template <typename MESH>
// uint32
// close(MESH& m, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

uint32 close(CMap2& m, bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		if (phi2(m, d) == d)
			fix_point_darts.push_back(d);

	for (Dart d : fix_point_darts)
	{
		if (phi2(m, d) == d)
		{
			CMap2::Face h = close_hole(m, d, set_indices);
			foreach_dart_of_orbit(m, h, [&](Dart hd) -> bool {
				set_boundary(m, hd, true);
				return true;
			});
			++nb_holes;
		}
	}

	return nb_holes;
}

/*****************************************************************************/

// template <typename MESH>
// void
// reverse_orientation(MESH& m);

/*****************************************************************************/

///////////
// CMap2 //
///////////

void reverse_orientation(CMap2& m)
{
	if (is_indexed<CMap2::Vertex>(m))
	{
		auto new_vertex_indices = m.darts_.add_attribute<uint32>("__new_vertex_indices");
		new_vertex_indices->fill(INVALID_INDEX);
		for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
			(*new_vertex_indices)[d.index] = index_of(m, CMap2::Vertex(phi1(m, d)));
		m.cells_indices_[CMap2::Vertex::ORBIT]->swap(new_vertex_indices.get());
		m.darts_.remove_attribute(new_vertex_indices);
	}

	m.phi1_->swap(m.phi_1_.get());
}

} // namespace cgogn
