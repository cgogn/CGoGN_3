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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP3_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP3_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/maps/cmap/cmap2.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap3 : public CMap2
{
	static const uint8 dimension = 3;

	using Parent = CMap2;

	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Cell<PHI2>;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Cell<PHI1_PHI2_PHI3>;

	using Cells = std::tuple<Vertex, Vertex2, HalfEdge, Edge, Edge2, Face, Face2, Volume>;

	std::shared_ptr<Attribute<Dart>> phi3_;

	CMap3() : CMap2()
	{
		phi3_ = add_relation("phi3");
	}
};

template <>
struct mesh_traits<CMap3>
{
	static constexpr const char* name = "CMap3";
	static constexpr const uint8 dimension = 3;

	using Parent = CMap3::Parent;

	using Vertex = CMap3::Vertex;
	using Vertex2 = CMap3::Vertex2;
	using HalfEdge = CMap3::HalfEdge;
	using Edge = CMap3::Edge;
	using Edge2 = CMap3::Edge2;
	using Face = CMap3::Face;
	using Face2 = CMap3::Face2;
	using Volume = CMap3::Volume;
	using CC = CMap3::CC;

	using Cells = std::tuple<Vertex, Vertex2, HalfEdge, Edge, Edge2, Face, Face2, Volume, CC>;
	static constexpr const char* cell_names[] = {"Vertex", "Vertex2", "HalfEdge", "Edge",
												 "Edge2",  "Face",	  "Face2",	  "Volume", "CC"};

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};


CMap3::Vertex CGOGN_CORE_EXPORT cut_edge(CMap3& m, CMap3::Edge e, bool set_indices = true);

CMap3::Edge CGOGN_CORE_EXPORT cut_face(CMap3& m, CMap3::Vertex v1, CMap3::Vertex v2, bool set_indices = true);

CMap3::Face CGOGN_CORE_EXPORT cut_volume(CMap3& m, const std::vector<Dart>& path, bool set_indices = true);

CMap3::Volume CGOGN_CORE_EXPORT close_hole(CMap3& m, Dart d, bool set_indices = true);

uint32 CGOGN_CORE_EXPORT close(CMap3& m, bool set_indices = true);

inline Dart phi3(const CMap3& m, Dart d)
{
	return (*(m.phi3_))[d.index];
}

void phi3_sew(CMap3& m, Dart d, Dart e);

void phi3_unsew(CMap3& m, Dart d);

bool check_integrity(CMap3& m, bool verbose = true);
} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP3_H_
