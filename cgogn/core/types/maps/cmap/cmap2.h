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

#ifndef CGOGN_CORE_TYPES_MAPS_CMAP_CMAP2_H_
#define CGOGN_CORE_TYPES_MAPS_CMAP_CMAP2_H_

#include <cgogn/core/types/maps/cmap/cmap1.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

struct CMap2 : public CMap1
{
	static const uint8 dimension = 2;

	using Parent = CMap1;

	using Vertex = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2>;
	using Face = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Volume;

	using BoundaryCell = Face;

	using Cells = std::tuple<Vertex, HalfEdge, Edge, Face, Volume>;

	std::shared_ptr<Attribute<Dart>> phi2_;

	CMap2() : CMap1()
	{
		phi2_ = add_relation("phi2");
	}
};

template <>
struct mesh_traits<CMap2>
{
	static constexpr const char* name = "CMap2";
	static constexpr const uint8 dimension = 2;

	using Parent = CMap2::Parent;

	using Vertex = CMap2::Vertex;
	using HalfEdge = CMap2::HalfEdge;
	using Edge = CMap2::Edge;
	using Face = CMap2::Face;
	using Volume = CMap2::Volume;

	using Cells = std::tuple<Vertex, HalfEdge, Edge, Face, Volume>;
	static constexpr const char* cell_names[] = {"Vertex", "HalfEdge", "Edge", "Face", "Volume"};

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

/*************************************************************************/
// Basic phi functions
/*************************************************************************/

inline Dart phi2(const CMap2& m, Dart d)
{
	return (*(m.phi2_))[d.index_];
}

inline void phi2_sew(CMap2& m, Dart d, Dart e)
{
	cgogn_assert(phi2(m, d) == d);
	cgogn_assert(phi2(m, e) == e);
	(*(m.phi2_))[d.index_] = e;
	(*(m.phi2_))[e.index_] = d;
}

inline void phi2_unsew(CMap2& m, Dart d)
{
	Dart e = phi2(m, d);
	(*(m.phi2_))[d.index_] = d;
	(*(m.phi2_))[e.index_] = e;
}

/*************************************************************************/
// Operators
/*************************************************************************/

CMap2::Vertex cut_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);
CMap2::Vertex collapse_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);
bool flip_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);

bool edge_can_collapse(const CMap2& m, CMap2::Edge e);
bool edge_can_flip(const CMap2& m, CMap2::Edge e);

CMap2::Face add_face(CMap2& m, uint32 size, bool set_indices = true);
void merge_incident_faces(CMap2& m, CMap2::Edge e, bool set_indices = true);
CMap2::Edge cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices = true);

CMap2::Volume add_pyramid(CMap2& m, uint32 size, bool set_indices = true);
CMap2::Volume add_prism(CMap2& m, uint32 size, bool set_indices = true);
void remove_volume(CMap2& m, CMap2::Volume v);

CMap2::Face close_hole(CMap2& m, Dart d, bool set_indices = true);
uint32 close(CMap2& m, bool set_indices = true);

void reverse_orientation(CMap2& m);

/*************************************************************************/
// High-level operators
/*************************************************************************/

CMap2::Vertex cut_edge_and_incident_triangles(CMap2& m, CMap2::Edge e);

/*************************************************************************/
// Specific accessors
/*************************************************************************/

// return the 4 vertices of a tetrahedron
std::array<CMap2::Vertex, 4> tet_vertices(CMap2& m, CMap2::Volume v);

// in a triangle mesh, return the opposite vertex of a halfedge
inline CMap2::Vertex opposite_vertex(CMap2& m, CMap2::HalfEdge he);

// in a triangle mesh, return the opposite vertices of an edge
inline std::vector<CMap2::Vertex> opposite_vertices(CMap2& m, CMap2::Edge e);

/*************************************************************************/
// Specific implementation of algorithms
/*************************************************************************/

// return the divergence at a vertex of a vector field defined on faces
geometry::Scalar face_vector_field_divergence(const CMap2& m, CMap2::Vertex v,
											  const CMap2::Attribute<geometry::Vec3>* face_gradient,
											  const CMap2::Attribute<geometry::Vec3>* vertex_position);

// return [horizon_halfedges, visible_faces]
// horizon_halfedges is a cycle of halfedges that are on the horizon
// visible_faces are the faces that are visible from the point, starting from the given face (that is supposed visible)
std::pair<std::vector<CMap2::HalfEdge>, std::vector<CMap2::Face>> build_horizon(
	CMap2& m, const CMap2::Attribute<geometry::Vec3>* vertex_position, const geometry::Vec3& point, CMap2::Face f);

// remove faces in visible_faces and fill the hole with a fan of triangles
// return the central vertex of the fan
CMap2::Vertex remove_faces_and_fill(CMap2& m, const std::vector<CMap2::HalfEdge>& area_boundary,
									const std::vector<CMap2::Face>& area_faces);

/*************************************************************************/
// Debugging helper functions
/*************************************************************************/

bool check_integrity(CMap2& m, bool verbose = true);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MAPS_CMAP_CMAP2_H_
