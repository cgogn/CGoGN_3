/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_

#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Edge> incident_edges(MESH& m, CELL c);

/*****************************************************************************/

///////////
// CMap1 //
///////////

std::vector<CMap1::Vertex> incident_edges(const CMap1& m, CMap1::Face f)
{
	std::vector<CMap1::Edge> edges;
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool { edges.push_back(CMap1::Edge(d)); return true; });
	return edges;
}

///////////
// CMap2 //
///////////

std::vector<CMap2::Edge> incident_edges(const CMap2& m, CMap2::Vertex v)
{
	std::vector<CMap2::Edge> edges;
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool { edges.push_back(CMap2::Edge(d)); return true; });
	return edges;
}

std::vector<CMap2::Edge> incident_edges(const CMap2& m, CMap2::Face f)
{
	std::vector<CMap2::Edge> edges;
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool { edges.push_back(CMap2::Edge(d)); return true; });
	return edges;
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
std::vector<typename mesh_traits<MESH>::Vertex>
incident_edges(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return incident_edges(m.mesh(), c);
}

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_edge(MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////
// CMap1 //
///////////

template <typename FUNC>
void foreach_incident_edge(const CMap1& m, CMap1::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap1::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool { return func(CMap1::Edge(d)); });
}

///////////
// CMap2 //
///////////

template <typename FUNC>
void foreach_incident_edge(const CMap2& m, CMap2::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool { return func(CMap2::Edge(d)); });
}

template <typename FUNC>
void foreach_incident_edge(const CMap2& m, CMap2::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool { return func(CMap2::Edge(d)); });
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
foreach_incident_edge(const MESH& m, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_incident_edge(m.mesh(), c, func);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_
