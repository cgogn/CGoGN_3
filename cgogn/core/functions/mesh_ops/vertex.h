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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_OPS_VERTEX_H_
#define CGOGN_CORE_FUNCTIONS_MESH_OPS_VERTEX_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// add_vertex(MESH& m, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Vertex CGOGN_CORE_EXPORT add_vertex(Graph& g, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// void
// remove_vertex(MESH& m, typename mesh_traits<MESH>::Vertex v, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

void CGOGN_CORE_EXPORT remove_vertex(Graph& g, Graph::Vertex v, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// connect_vertices(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool
// set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Edge CGOGN_CORE_EXPORT connect_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// void
// disconnect_vertices(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

void CGOGN_CORE_EXPORT disconnect_vertices(Graph& g, Graph::Edge e, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// merge_vertices(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool
// set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

void CGOGN_CORE_EXPORT merge_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices = true);

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_VERTEX_H_
