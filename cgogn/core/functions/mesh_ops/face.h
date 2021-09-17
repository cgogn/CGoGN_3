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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_OPS_FACE_H_
#define CGOGN_CORE_FUNCTIONS_MESH_OPS_FACE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// add_face(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

CMap1::Face CGOGN_CORE_EXPORT add_face(CMap1& m, uint32 size, bool set_indices = true);

///////////
// CMap2 //
///////////

CMap2::Face CGOGN_CORE_EXPORT add_face(CMap2& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// add_face(MESH& m, std::vector<typename mesh_traits<MESH>::Edge edges);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

IncidenceGraph::Face CGOGN_CORE_EXPORT add_face(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge>& edges);

/*****************************************************************************/

// template <typename MESH>
// void
// remove_face(MESH& m, typename mesh_traits<MESH>::Face f);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

void CGOGN_CORE_EXPORT remove_face(IncidenceGraph& ig, IncidenceGraph::Face f);

///////////
// CMap1 //
///////////

void CGOGN_CORE_EXPORT remove_face(CMap1& m, CMap1::Face f, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// void
// merge_incident_faces(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

void CGOGN_CORE_EXPORT merge_incident_faces(CMap2& m, CMap2::Edge e, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// cut_face(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_indices =
// true);

/*****************************************************************************/

////////////////////
// IncidenceGraph //
////////////////////

IncidenceGraph::Edge CGOGN_CORE_EXPORT cut_face(IncidenceGraph& m, IncidenceGraph::Vertex v1,
												IncidenceGraph::Vertex v2);

///////////
// CMap2 //
///////////

CMap2::Edge CGOGN_CORE_EXPORT cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices = true);

///////////
// CMap3 //
///////////

CMap3::Edge CGOGN_CORE_EXPORT cut_face(CMap3& m, CMap3::Vertex v1, CMap3::Vertex v2, bool set_indices = true);

//////////
// CPH3 //
//////////

CPH3::CMAP::Edge CGOGN_CORE_EXPORT cut_face(CPH3& m, CPH3::CMAP::Vertex v1, CPH3::CMAP::Vertex v2,
											bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// close_hole(MESH& m, Dart d, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Face close_hole(CMap2& m, Dart d, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// uint32
// close(MESH& m, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

uint32 close(CMap2& m, bool set_indices = true);

/*****************************************************************************/

// template <typename MESH>
// void
// reverse_orientation(MESH& m);

/*****************************************************************************/

///////////
// CMap2 //
///////////

void reverse_orientation(CMap2& m);

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_FACE_H_
