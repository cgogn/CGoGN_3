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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_OPS_EDGE_H_
#define CGOGN_CORE_FUNCTIONS_MESH_OPS_EDGE_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/cmap/cmap_ops.h>

#include <string>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// cut_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

CMap1::Vertex
cut_edge(CMap1& m, CMap1::Edge e, bool set_indices = true)
{
	Dart d = m.add_dart();
	m.phi1_sew(e.dart, d);
	CMap1::Vertex v(d);

	if (set_indices)
	{
		if (m.is_embedded<CMap1::Vertex>())
			create_embedding(m, v);
		if (m.is_embedded<CMap1::Face>())
			m.copy_embedding<CMap1::Face>(d, e.dart);
	}

	return v;
}

///////////
// CMap2 //
///////////

CMap2::Vertex
cut_edge(CMap2& m, CMap2::Edge e, bool set_indices = true)
{
	Dart d1 = e.dart;
	Dart d2 = m.phi2(d1);
	m.phi2_unsew(d1);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d1), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d2), false);
	m.phi2_sew(d1, nv2.dart);
	m.phi2_sew(d2, nv1.dart);
	m.set_boundary(nv1.dart, m.is_boundary(d1));
	m.set_boundary(nv2.dart, m.is_boundary(d2));
	CMap2::Vertex v(nv1.dart);

	if (set_indices)
	{
		if (m.is_embedded<CMap2::Vertex>())
			create_embedding(m, v);
		if (m.is_embedded<CMap2::Edge>())
		{
			m.copy_embedding<CMap2::Edge>(nv2.dart, d1);
			create_embedding(m, CMap2::Edge(nv1.dart));
		}
		if (m.is_embedded<CMap2::Face>())
		{
			m.copy_embedding<CMap2::Face>(nv1.dart, d1);
			m.copy_embedding<CMap2::Face>(nv2.dart, d2);
		}
		if (m.is_embedded<CMap2::Volume>())
		{
			m.copy_embedding<CMap2::Volume>(nv1.dart, d1);
			m.copy_embedding<CMap2::Volume>(nv2.dart, d2);
		}
	}

	return v;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::Vertex
cut_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true)
{
	return cut_edge(m.mesh(), e, set_indices);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_EDGE_H_
