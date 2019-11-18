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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_HALFEDGE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_HALFEDGE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::HalfEdge> incident_halfedges(MESH& m, CELL c);

/*****************************************************************************/

///////////
// Graph //
///////////

std::vector<Graph::HalfEdge>
CGOGN_CORE_EXPORT incident_halfedges(const Graph& g, Graph::Edge e);

std::vector<Graph::HalfEdge>
CGOGN_CORE_EXPORT incident_halfedges(const Graph& g, Graph::Vertex v);

///////////
// CMap2 //
///////////

std::vector<CMap2::HalfEdge>
CGOGN_CORE_EXPORT incident_halfedges(const CMap2& m, CMap2::Edge e);

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
std::vector<typename mesh_traits<MESH>::HalfEdge>
incident_halfedges(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return incident_edges(m.mesh(), c);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_HALFEDGE_H_
