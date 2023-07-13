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

#ifndef CGOGN_CORE_TYPES_MAPS_CMAP_GRAPH_H_
#define CGOGN_CORE_TYPES_MAPS_CMAP_GRAPH_H_

#include <cgogn/core/types/maps/cmap/cmap_base.h>

namespace cgogn
{

struct Graph : public CMapBase
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

/*************************************************************************/
// Basic alpha functions
/*************************************************************************/

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

inline void alpha0_sew(Graph& m, Dart d, Dart e)
{
	(*m.alpha0_)[d.index] = e;
	(*m.alpha0_)[e.index] = d;
}

inline void alpha0_unsew(Graph& m, Dart d)
{
	Dart e = alpha0(m, d);
	(*m.alpha0_)[d.index] = d;
	(*m.alpha0_)[e.index] = e;
}

inline void alpha1_sew(Graph& m, Dart d, Dart e)
{
	Dart f = alpha1(m, d);
	Dart g = alpha1(m, e);
	(*m.alpha1_)[d.index] = g;
	(*m.alpha1_)[e.index] = f;
	(*m.alpha_1_)[g.index] = d;
	(*m.alpha_1_)[f.index] = e;
}

inline void alpha1_unsew(Graph& m, Dart d)
{
	Dart e = alpha1(m, d);
	Dart f = alpha_1(m, d);
	(*m.alpha1_)[f.index] = e;
	(*m.alpha1_)[d.index] = d;
	(*m.alpha_1_)[e.index] = f;
	(*m.alpha_1_)[d.index] = d;
}

/*************************************************************************/
// Orbit traversals
/*************************************************************************/

template <typename FUNC>
void foreach_dart_of_ALPHA0(const Graph& m, Dart d, const FUNC& f)
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (f(d))
		f(alpha0(m, d));
}

template <typename FUNC>
void foreach_dart_of_ALPHA1(const Graph& m, Dart d, const FUNC& f)
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			break;
		it = alpha1(m, it);
	} while (it != d);
}

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

/*************************************************************************/
// Operators
/*************************************************************************/

Graph::Vertex add_vertex(Graph& g, bool set_indices = true);
void remove_vertex(Graph& g, Graph::Vertex v, bool set_indices = true);
void merge_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices = true);
Graph::Edge connect_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices = true);
void disconnect_vertices(Graph& g, Graph::Edge e, bool set_indices = true);

inline bool is_vertex_isolated(Graph& g, Graph::Vertex v)
{
	return alpha0(g, v.dart) == alpha1(g, v.dart);
};

Graph::Vertex cut_edge(Graph& m, Graph::Edge e, bool set_indices = true);
Graph::Vertex collapse_edge(Graph& g, Graph::Edge e, bool set_indices = true);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MAPS_CMAP_GRAPH_H_
