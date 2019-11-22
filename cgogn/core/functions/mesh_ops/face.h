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

CMap1::Face
CGOGN_CORE_EXPORT add_face(CMap1& m, uint32 size, bool set_indices = true);

///////////
// CMap2 //
///////////

CMap2::Face
CGOGN_CORE_EXPORT add_face(CMap2& m, uint32 size, bool set_indices = true);

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::Face
add_face(MESH& m, uint32 size, bool set_indices = true)
{
	return add_face(m.mesh(), size, set_indices);
}

/*****************************************************************************/

// template <typename MESH>
// void
// remove_face(MESH& m, typename mesh_traits<MESH>::Face f, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

void
CGOGN_CORE_EXPORT remove_face(CMap1& m, CMap1::Face f, bool set_indices = true);

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
remove_face(MESH& m, typename mesh_traits<MESH>::Face f, bool set_indices = true)
{
	return remove_face(m.mesh(), f, set_indices);
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// cut_face(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Edge
CGOGN_CORE_EXPORT cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices = true);

///////////
// CMap3 //
///////////

CMap3::Edge
CGOGN_CORE_EXPORT cut_face(CMap3& m, CMap3::Vertex v1, CMap3::Vertex v2, bool set_indices = true);

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::Edge
cut_face(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_indices = true)
{
	return cut_face(m.mesh(), v1, v2, set_indices);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_FACE_H_
