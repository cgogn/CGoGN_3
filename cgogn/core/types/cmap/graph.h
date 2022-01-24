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

#ifndef CGOGN_CORE_TYPES_CMAP_GRAPH_H_
#define CGOGN_CORE_TYPES_CMAP_GRAPH_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cell.h>
#include <cgogn/core/types/cmap/cmap_base.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT Graph : public CMapBase
{
	static const uint8 dimension = 1;

	using Vertex = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2>;

	using Cells = std::tuple<Vertex, HalfEdge, Edge>;

	std::shared_ptr<Attribute<Dart>> alpha0_;
	std::shared_ptr<Attribute<Dart>> alpha1_;
	std::shared_ptr<Attribute<Dart>> alpha_1_;

	Graph()
	{
		alpha0_ = add_relation("alpha0");
		alpha1_ = add_relation("alpha1");
		alpha_1_ = add_relation("alpha_1");
	}
};

template <>
struct mesh_traits<Graph>
{
	static constexpr const char* name = "Graph";
	static constexpr const uint8 dimension = 1;

	using Vertex = Graph::Vertex;
	using HalfEdge = Graph::HalfEdge;
	using Edge = Graph::Edge;

	using Cells = std::tuple<Vertex, HalfEdge, Edge>;
	static constexpr const char* cell_names[] = {"Vertex", "HalfEdge", "Edge"};

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_GRAPH_H_
