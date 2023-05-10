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

Graph::Vertex CGOGN_CORE_EXPORT add_vertex(Graph& g, bool set_indices = true);

void CGOGN_CORE_EXPORT remove_vertex(Graph& g, Graph::Vertex v, bool set_indices = true);

Graph::Edge CGOGN_CORE_EXPORT connect_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices = true);

void CGOGN_CORE_EXPORT disconnect_vertices(Graph& g, Graph::Edge e, bool set_indices = true);

void CGOGN_CORE_EXPORT merge_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices = true);

Graph::Vertex CGOGN_CORE_EXPORT cut_edge(Graph& m, Graph::Edge e, bool set_indices = true);

Graph::Vertex CGOGN_CORE_EXPORT collapse_edge(Graph& g, Graph::Edge e, bool set_indices = true);


template <typename CELL, typename FUNC>
void foreach_dart_of_orbit(const Graph& m, CELL c, const FUNC& f)
{
	static_assert(is_in_tuple<CELL, typename Graph::Cells>::value, "Cell not supported in a Graph");
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	static const Orbit orbit = CELL::ORBIT;
	switch (orbit)
	{
	case DART:
		f(c.dart);
		break;
	case PHI2:
		foreach_dart_of_ALPHA0(m, c.dart, f);
		break;
	case PHI21:
		foreach_dart_of_ALPHA1(m, c.dart, f);
		break;
	default:
		break;
	}
}

inline Dart alpha0(const Graph& m, Dart d)
{
	return (*m.alpha0_)[d.index];
}

inline Dart alpha1(const Graph& m, Dart d)
{
	return (*m.alpha1_)[d.index];
}

inline Dart alpha_1(const Graph& m, Dart d)
{
	return (*m.alpha_1_)[d.index];
}

void alpha0_sew(Graph& m, Dart d, Dart e);

void alpha0_unsew(Graph& m, Dart d);

void alpha1_sew(Graph& m, Dart d, Dart e);

void alpha1_unsew(Graph& m, Dart d);

inline bool is_vertex_isolated(Graph& g, Graph::Vertex v)
{
	return alpha0(g, v.dart) == alpha1(g, v.dart); 
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_GRAPH_H_
