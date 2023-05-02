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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP1_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP1_H_

#include <cgogn/core/types/cmap/cmap0.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap1 : public CMap0
{
	static const uint8 dimension = 1;

	using Parent = CMap0;

	using Vertex = Cell<DART>;
	using Edge = Cell<DART>;
	using Face = Cell<PHI1>;

	using Cells = std::tuple<Vertex, Edge, Face>;

	std::shared_ptr<Attribute<Dart>> phi1_;
	std::shared_ptr<Attribute<Dart>> phi_1_;

	CMap1() : CMap0()
	{
		phi1_ = add_relation("phi1");
		phi_1_ = add_relation("phi_1");
	}
};

template <>
struct mesh_traits<CMap1>
{
	static constexpr const char* name = "CMap1";
	static constexpr const uint8 dimension = 1;

	using Parent = CMap1::Parent;

	using Vertex = CMap1::Vertex;
	using Edge = CMap1::Edge;
	using Face = CMap1::Face;

	using Cells = std::tuple<Vertex, Edge, Face>;
	static constexpr const char* cell_names[] = {"Vertex", "Edge", "Face"};

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

CMap1::Vertex CGOGN_CORE_EXPORT cut_edge(CMap1& m, CMap1::Edge e, bool set_indices = true);

CMap1::Vertex CGOGN_CORE_EXPORT collapse_edge(CMap1& m, CMap1::Edge e, bool set_indices = true);

CMap1::Face CGOGN_CORE_EXPORT add_face(CMap1& m, uint32 size, bool set_indices = true);

void CGOGN_CORE_EXPORT remove_face(CMap1& m, CMap1::Face f);

bool check_integrity(CMap1& m, bool verbose = true);

inline Dart phi1(const CMap1& m, Dart d)
{
	return (*(m.phi1_))[d.index];
}

inline Dart phi_1(const CMap1& m, Dart d)
{
	return (*(m.phi_1_))[d.index];
}

void phi1_sew(CMap1& m, Dart d, Dart e);

void phi1_unsew(CMap1& m, Dart d);


} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP1_H_
