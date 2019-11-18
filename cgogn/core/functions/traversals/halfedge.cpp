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

#include <cgogn/core/functions/traversals/halfedge.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::vector<typename mesh_traits<MESH>::HalfEdge> incident_halfedges(const MESH& m, CELL c);

/*****************************************************************************/

///////////
// Graph //
///////////

std::vector<Graph::HalfEdge> incident_halfedges(const Graph& g, Graph::Edge e)
{
	std::vector<Graph::HalfEdge> halfedges;
	halfedges.reserve(2u);
	halfedges.push_back(Graph::HalfEdge(e.dart));
	halfedges.push_back(Graph::HalfEdge(g.alpha0(e.dart)));
	return halfedges;
}

std::vector<Graph::HalfEdge> incident_halfedges(const Graph& g, Graph::Vertex v)
{
	std::vector<Graph::HalfEdge> halfedges;
	halfedges.reserve(8u);
	g.foreach_dart_of_orbit(v, [&] (Dart d) -> bool { halfedges.push_back(Graph::HalfEdge(d)); return true; });
	return halfedges;
}

///////////
// CMap2 //
///////////

std::vector<CMap2::HalfEdge>
CGOGN_CORE_EXPORT incident_halfedges(const CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::HalfEdge> edges1;
	edges1.reserve(2u);
	edges1.push_back(CMap2::HalfEdge(e.dart));
	edges1.push_back(CMap2::HalfEdge(m.phi2(e.dart)));
	return edges1;
}

} // namespace cgogn
