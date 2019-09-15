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
#include <cgogn/core/types/cmap/cmap_ops.h>
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
	Dart base = m.close_hole(first.dart); // Add the base face
	
	CMap2::Volume vol(base);

	if (set_indices)
	{
		if (m.is_embedded<CMap2::Vertex>())
			foreach_incident_vertex(m, vol, [&] (CMap2::Vertex v) -> bool { create_embedding(m, v); return true; });
		if (m.is_embedded<CMap2::Edge>())
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool { create_embedding(m, e); return true; });
		if (m.is_embedded<CMap2::Face>())
			foreach_incident_face(m, vol, [&] (CMap2::Face f) -> bool { create_embedding(m, f); return true; });
		if (m.is_embedded<CMap2::Volume>())
			create_embedding(m, vol);
	}

	return vol;
}

} // namespace cgogn
