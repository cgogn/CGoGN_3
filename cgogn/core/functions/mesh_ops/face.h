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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/cmap/cmap_ops.h>

#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/mesh_ops/edge.h>

#include <string>

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
add_face(CMap1& m, uint32 size, bool set_indices = true)
{
	Dart d = m.add_dart();
	for (uint32 i = 1u; i < size; ++i)
	{
		Dart e = m.add_dart();
		m.phi1_sew(d, e);
	}
	CMap1::Face f(d);

	if (set_indices)
	{
		if (m.is_embedded<CMap1::Vertex>())
			foreach_incident_vertex(m, f, [&] (CMap1::Vertex v) -> bool { create_embedding(m, v); return true; });
		if (m.is_embedded<CMap1::Face>())
			create_embedding(m, f);
	}

	return f;
}

///////////
// CMap2 //
///////////

CMap2::Face
add_face(CMap2& m, uint32 size, bool set_indices = true)
{
	CMap2::Face f = add_face(static_cast<CMap1&>(m), size, false);
	CMap2::Face b = add_face(static_cast<CMap1&>(m), size, false);
	Dart it = b.dart;
	m.foreach_dart_of_orbit(f, [&] (Dart d) -> bool
	{
		m.set_boundary(it, true);
		m.phi2_sew(d, it);
		it = m.phi_1(it);
		return true;
	});

	if (set_indices)
	{
		if (m.is_embedded<CMap2::Vertex>())
			foreach_incident_vertex(m, f, [&] (CMap2::Vertex v) -> bool { create_embedding(m, v); return true; });
		if (m.is_embedded<CMap2::Edge>())
			foreach_incident_edge(m, f, [&] (CMap2::Edge e) -> bool { create_embedding(m, e); return true; });
		if (m.is_embedded<CMap2::Face>())
			create_embedding(m, f);
	}

	return f;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::Face
add_face(MESH& m, uint32 size, bool set_incides = true)
{
	return add_face(m.mesh(), size, set_incides);
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
cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices = true)
{
	Dart dd = m.phi_1(v1.dart);
	Dart ee = m.phi_1(v2.dart);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);
	m.phi1_sew(nv1.dart, nv2.dart);
	m.phi2_sew(nv1.dart, nv2.dart);
	m.set_boundary(nv1.dart, m.is_boundary(dd));
	m.set_boundary(nv2.dart, m.is_boundary(ee));
	CMap2::Edge e(nv1.dart);

	if (set_indices)
	{
		if (m.is_embedded<CMap2::Vertex>())
		{
			m.copy_embedding<CMap2::Vertex>(nv1.dart, v1.dart);
			m.copy_embedding<CMap2::Vertex>(nv2.dart, v2.dart);
		}
		if (m.is_embedded<CMap2::Edge>())
			create_embedding(m, CMap2::Edge(nv1.dart));
		if (m.is_embedded<CMap2::Face>())
		{
			m.copy_embedding<CMap2::Face>(nv2.dart, v1.dart);
			create_embedding(m, CMap2::Face(v2.dart));
		}
		if (m.is_embedded<CMap2::Volume>())
		{
			m.copy_embedding<CMap2::Volume>(nv1.dart, v2.dart);
			m.copy_embedding<CMap2::Volume>(nv2.dart, v1.dart);
		}
	}

	return e;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::Edge
cut_face(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_incides = true)
{
	return cut_face(m.mesh(), v1, v2, set_incides);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_H_
