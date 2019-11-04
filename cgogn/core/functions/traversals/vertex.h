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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::vector<typename mesh_traits<MESH>::Vertex> incident_vertices(const MESH& m, CELL c);

/*****************************************************************************/

///////////
// Graph //
///////////

std::vector<Graph::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const Graph& g, Graph::Edge e);

///////////
// CMap1 //
///////////

std::vector<CMap1::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap1& m, CMap1::Face f);

///////////
// CMap2 //
///////////

std::vector<CMap2::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap2& m, CMap2::Edge e);

std::vector<CMap2::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap2& m, CMap2::Face f);

std::vector<CMap2::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap2& m, CMap2::Volume v);

///////////
// CMap3 //
///////////

std::vector<CMap3::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap3& m, CMap3::Edge e);

std::vector<CMap3::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap3& m, CMap3::Face f);

std::vector<CMap3::Vertex>
CGOGN_CORE_EXPORT incident_vertices(const CMap3& m, CMap3::Volume v);

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
std::vector<typename mesh_traits<MESH>::Vertex>
incident_vertices(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return incident_vertices(m.mesh(), c);
}

/*****************************************************************************/

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_incident_vertex(const MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////
// Graph //
///////////

template <typename FUNC>
void foreach_incident_vertex(const Graph& g, Graph::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, Graph::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	g.foreach_dart_of_orbit(e, [&] (Dart d) -> bool { return func(Graph::Vertex(d)); });
}

///////////
// CMap1 //
///////////

template <typename FUNC>
void foreach_incident_vertex(const CMap1& m, CMap1::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap1::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool { return func(CMap1::Vertex(d)); });
}

///////////
// CMap2 //
///////////

template <typename FUNC>
void foreach_incident_vertex(const CMap2& m, CMap2::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(e, [&] (Dart d) -> bool { return func(CMap2::Vertex(d)); });
}

template <typename FUNC>
void foreach_incident_vertex(const CMap2& m, CMap2::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool { return func(CMap2::Vertex(d)); });
}

template <typename FUNC>
void foreach_incident_vertex(const CMap2& m, CMap2::Volume v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	DartMarkerStore marker(m);
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			m.foreach_dart_of_orbit(CMap2::Vertex(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap2::Vertex(d));
		}
		return true;
	});
}

///////////
// CMap3 //
///////////

template <typename FUNC>
void foreach_incident_vertex(const CMap3& m, CMap3::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	static_cast<const CMap2&>(m).foreach_dart_of_orbit(CMap2::Edge(e.dart), [&] (Dart d) -> bool { return func(CMap3::Vertex(d)); });
}

template <typename FUNC>
void foreach_incident_vertex(const CMap3& m, CMap3::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	static_cast<const CMap2&>(m).foreach_dart_of_orbit(CMap2::Face(f.dart), [&] (Dart d) -> bool { return func(CMap3::Vertex(d)); });
}

template <typename FUNC>
void foreach_incident_vertex(const CMap3& m, CMap3::Volume v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	DartMarkerStore marker(m);
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			static_cast<const CMap2&>(m).foreach_dart_of_orbit(CMap2::Vertex(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap3::Vertex(d));
		}
		return true;
	});
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
foreach_incident_vertex(const MESH& m, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_incident_vertex(m.mesh(), c, func);
}

/*****************************************************************************/

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_adjacent_vertex_through_edge(MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& f);

/*****************************************************************************/

///////////
// Graph //
///////////

template <typename FUNC>
void
foreach_adjacent_vertex_through_edge(const Graph& g, Graph::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, Graph::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	g.foreach_dart_of_orbit(v, [&] (Dart d) -> bool { return func(Graph::Vertex(g.alpha0(d))); });
}

///////////
// CMap2 //
///////////

template <typename FUNC>
void
foreach_adjacent_vertex_through_edge(const CMap2& m, CMap2::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool { return func(CMap2::Vertex(m.phi2(d))); });
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
foreach_adjacent_vertex_through_edge(const MESH& m, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_adjacent_vertex_through_edge(m.mesh(), c, func);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_VERTEX_H_
