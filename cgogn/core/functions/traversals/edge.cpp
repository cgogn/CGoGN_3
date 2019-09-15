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

#include <cgogn/core/functions/traversals/edge.h>

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

std::vector<CMap2::Edge> incident_edges(const CMap2& m, CMap2::Volume v)
{
	std::vector<CMap2::Edge> edges;
	foreach_incident_edge(m, v, [&] (CMap2::Edge e) -> bool { edges.push_back(e); return true; });
	return edges;
}

///////////
// CMap3 //
///////////

std::vector<CMap3::Edge> incident_edges(const CMap3& m, CMap3::Vertex v)
{
	std::vector<CMap3::Edge> edges;
	foreach_incident_edge(m, v, [&] (CMap3::Edge e) -> bool { edges.push_back(e); return true; });
	return edges;
}

std::vector<CMap3::Edge> incident_edges(const CMap3& m, CMap3::Face f)
{
	std::vector<CMap3::Edge> edges;
	static_cast<const CMap2&>(m).foreach_dart_of_orbit(CMap2::Face(f.dart), [&] (Dart d) -> bool { edges.push_back(CMap3::Edge(d)); return true; });
	return edges;
}

std::vector<CMap3::Edge> incident_edges(const CMap3& m, CMap3::Volume v)
{
	std::vector<CMap3::Edge> edges;
	foreach_incident_edge(m, v, [&] (CMap3::Edge e) -> bool { edges.push_back(e); return true; });
	return edges;
}

} // namespace cgogn
