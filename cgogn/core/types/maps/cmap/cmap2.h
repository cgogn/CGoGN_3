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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP2_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP2_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/maps/cmap/cmap1.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap2 : public CMap1
{
	static const uint8 dimension = 2;

	using Parent = CMap1;

	using Vertex = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2>;
	using Face = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Volume;

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

CMap2::Vertex CGOGN_CORE_EXPORT cut_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);

CMap2::Vertex CGOGN_CORE_EXPORT collapse_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);

bool CGOGN_CORE_EXPORT flip_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);

CMap2::Face CGOGN_CORE_EXPORT add_face(CMap2& m, uint32 size, bool set_indices = true);

CMap2::Volume CGOGN_CORE_EXPORT add_pyramid(CMap2& m, uint32 size, bool set_indices = true);

CMap2::Volume CGOGN_CORE_EXPORT add_prism(CMap2& m, uint32 size, bool set_indices = true);

void CGOGN_CORE_EXPORT remove_volume(CMap2& m, CMap2::Volume v);

void CGOGN_CORE_EXPORT merge_incident_faces(CMap2& m, CMap2::Edge e, bool set_indices = true);

CMap2::Edge CGOGN_CORE_EXPORT cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices = true);

CMap2::Face close_hole(CMap2& m, Dart d, bool set_indices = true);

uint32 close(CMap2& m, bool set_indices = true);

void reverse_orientation(CMap2& m);


bool check_integrity(CMap2& m, bool verbose = true);

bool edge_can_collapse(const CMap2& m, CMap2::Edge e);

bool edge_can_flip(const CMap2& m, CMap2::Edge e);


inline Dart phi2(const CMap2& m, Dart d)
{
	return (*(m.phi2_))[d.index];
}

void phi2_sew(CMap2& m, Dart d, Dart e);

void phi2_unsew(CMap2& m, Dart d);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP2_H_
