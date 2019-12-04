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

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>
#include <cgogn/core/functions/traversals/dart.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Edge> incident_edges(MESH& m, CELL c);

/*****************************************************************************/

///////////
// Graph //
///////////

std::vector<Graph::Edge>
CGOGN_CORE_EXPORT incident_edges(const Graph& g, Graph::Vertex v);

///////////
// CMap1 //
///////////

std::vector<CMap1::Vertex>
CGOGN_CORE_EXPORT incident_edges(const CMap1& m, CMap1::Face f);

///////////
// CMap2 //
///////////

std::vector<CMap2::Edge>
CGOGN_CORE_EXPORT incident_edges(const CMap2& m, CMap2::Vertex v);

CMap2::Edge
CGOGN_CORE_EXPORT incident_edge(const CMap2& m, CMap2::HalfEdge h);

std::vector<CMap2::Edge>
CGOGN_CORE_EXPORT incident_edges(const CMap2& m, CMap2::Face f);

std::vector<CMap2::Edge>
CGOGN_CORE_EXPORT incident_edges(const CMap2& m, CMap2::Volume v);

///////////
// CMap3 //
///////////

std::vector<CMap3::Edge>
CGOGN_CORE_EXPORT incident_edges(const CMap3& m, CMap3::Vertex v);

std::vector<CMap3::Edge>
CGOGN_CORE_EXPORT incident_edges(const CMap3& m, CMap3::Face f);

std::vector<CMap3::Edge>
CGOGN_CORE_EXPORT incident_edges(const CMap3& m, CMap3::Volume v);

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
std::vector<typename mesh_traits<MESH>::Edge>
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
// Graph //
///////////

template <typename FUNC>
void foreach_incident_edge(const Graph& m, Graph::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, Graph::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool { return func(Graph::Edge(d)); });
}

///////////
// CMap1 //
///////////

template <typename FUNC>
void foreach_incident_edge(const CMap1& m, CMap1::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap1::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m,f, [&] (Dart d) -> bool { return func(CMap1::Edge(d)); });
}

///////////
// CMap2 //
///////////

template <typename FUNC>
void foreach_incident_edge(const CMap2& m, CMap2::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool { return func(CMap2::Edge(d)); });
}

template <typename FUNC>
void foreach_incident_edge(const CMap2& m, CMap2::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m,f, [&] (Dart d) -> bool { return func(CMap2::Edge(d)); });
}

template <typename FUNC>
void foreach_incident_edge(const CMap2& m, CMap2::Volume v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	DartMarkerStore marker(m);
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			foreach_dart_of_orbit(m,CMap2::Edge(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap2::Edge(d));
		}
		return true;
	});
}

///////////
// CMap3 //
///////////

template <typename FUNC>
void foreach_incident_edge(const CMap3& m, CMap3::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	DartMarkerStore marker(m);
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			foreach_dart_of_orbit(m,CMap3::Edge(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap3::Edge(d));
		}
		return true;
	});
}

template <typename FUNC>
void foreach_incident_edge(const CMap3& m, CMap3::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m,CMap3::Face2(f.dart), [&] (Dart d) -> bool { return func(CMap3::Edge(d)); });
}

template <typename FUNC>
void foreach_incident_edge(const CMap3& m, CMap3::Volume v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	DartMarkerStore marker(m);
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			// TODO: could mark only the darts of CMap2::Edge(d)
			foreach_dart_of_orbit(m,CMap3::Edge(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap3::Edge(d));
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
foreach_incident_edge(const MESH& m, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_incident_edge(m.mesh(), c, func);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_EDGE_H_
