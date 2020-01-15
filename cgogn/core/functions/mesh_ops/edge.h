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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_OPS_EDGE_H_
#define CGOGN_CORE_FUNCTIONS_MESH_OPS_EDGE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// cut_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Vertex CGOGN_CORE_EXPORT cut_edge(Graph& m, Graph::Edge e, bool set_indices = true);

///////////
// CMap1 //
///////////

CMap1::Vertex CGOGN_CORE_EXPORT cut_edge(CMap1& m, CMap1::Edge e, bool set_indices = true);

///////////
// CMap2 //
///////////

CMap2::Vertex CGOGN_CORE_EXPORT cut_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);

///////////
// CMap3 //
///////////

CMap3::Vertex CGOGN_CORE_EXPORT cut_edge(CMap3& m, CMap3::Edge e, bool set_indices = true);

//////////
// CPH3 //
//////////

CPH3::CMAP::Vertex CGOGN_CORE_EXPORT cut_edge(CPH3& m, CPH3::CMAP::Edge e, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// collapse_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

CMap1::Vertex CGOGN_CORE_EXPORT collapse_edge(CMap1& m, CMap1::Edge e, bool set_indices = true);

///////////
// CMap2 //
///////////

CMap2::Vertex CGOGN_CORE_EXPORT collapse_edge(CMap2& m, CMap2::Edge e, bool set_indices = true);

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_EDGE_H_
