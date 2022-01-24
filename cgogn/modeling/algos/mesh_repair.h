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

#ifndef CGOGN_MODELING_ALGOS_MESH_REPAIR_H_
#define CGOGN_MODELING_ALGOS_MESH_REPAIR_H_

#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

namespace cgogn
{

namespace modeling
{

///////////
// CMap2 //
///////////

inline void fill_holes(CMap2& m, bool set_indices = true)
{
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		if (is_boundary(m, d))
		{
			set_boundary(m, d, false);
			if (set_indices)
			{
				if (is_indexed<CMap2::Face>(m))
				{
					if (index_of(m, CMap2::Face(d)) == INVALID_INDEX)
						set_index(m, CMap2::Face(d), new_index<CMap2::Face>(m));
				}
			}
		}
	}
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void remove_small_components(MESH& m, uint32 min_vertices)
{
	if constexpr (mesh_traits<MESH>::dimension == 2)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Volume = typename mesh_traits<MESH>::Volume;
		foreach_cell(m, [&](Volume vol) -> bool {
			uint32 n = 0;
			foreach_incident_vertex(m, vol, [&](Vertex /*v*/) -> bool {
				++n;
				return true;
			});
			if (n < min_vertices)
				remove_volume(m, vol);
			return true;
		});
	}
	else if constexpr (mesh_traits<MESH>::dimension == 3)
	{
		// TODO
	}
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_MESH_REPAIR_H_
