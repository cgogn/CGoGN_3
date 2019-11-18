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

#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// add_pyramid(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Volume
add_pyramid(CMap2& m, uint32 size, bool set_indices)
{
	CMap1::Face first = add_face(static_cast<CMap1&>(m), 3u, false); // First triangle
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next triangles
	{
		CMap1::Face next = add_face(static_cast<CMap1&>(m), 3u, false);
		m.phi2_sew(m.phi_1(current), m.phi1(next.dart));
		current = next.dart;
	}
	m.phi2_sew(m.phi_1(current), m.phi1(first.dart)); // Finish the umbrella
	CMap2::Face base = m.close_hole(first.dart, false); // Add the base face
	
	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (m.is_indexed<CMap2::Vertex>())
		{
			foreach_incident_vertex(m, vol, [&] (CMap2::Vertex v) -> bool
			{
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::HalfEdge>())
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				set_index(m, CMap2::HalfEdge(m.phi2(e.dart)), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::Edge>())
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::Face>())
		{
			foreach_incident_face(m, vol, [&] (CMap2::Face f) -> bool
			{
				set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::Volume>())
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

CMap2::Volume
add_prism(CMap2& m, uint32 size, bool set_indices)
{
	CMap1::Face first = add_face(static_cast<CMap1&>(m), 4u, false); // first quad
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next quads
	{
		CMap1::Face next = add_face(static_cast<CMap1&>(m), 4u, false);
		m.phi2_sew(m.phi_1(current), m.phi1(next.dart));
		current = next.dart;
	}
	m.phi2_sew(m.phi_1(current), m.phi1(first.dart)); // Finish the sides
	CMap2::Face base = m.close_hole(first.dart, false); // Add the base face
	CMap2::Face top = m.close_hole(m.phi<11>(first.dart), false); // Add the top face
	
	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (m.is_indexed<CMap2::Vertex>())
		{
			foreach_incident_vertex(m, vol, [&] (CMap2::Vertex v) -> bool
			{
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::HalfEdge>())
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				set_index(m, CMap2::HalfEdge(m.phi2(e.dart)), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::Edge>())
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::Face>())
		{
			foreach_incident_face(m, vol, [&] (CMap2::Face f) -> bool
			{
				set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			});
		}
		if (m.is_indexed<CMap2::Volume>())
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

} // namespace cgogn
