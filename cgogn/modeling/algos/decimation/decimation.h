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

#ifndef CGOGN_GEOMETRY_ALGOS_DECIMATION_H_
#define CGOGN_GEOMETRY_ALGOS_DECIMATION_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/mesh_ops/edge.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/decimation/edge_approximator_mid_edge.h>
#include <cgogn/modeling/algos/decimation/edge_traversor_edge_length.h>

namespace cgogn
{

namespace modeling
{

/////////////
// GENERIC //
/////////////

template <typename VEC, typename MESH>
void decimate(MESH& m, typename mesh_traits<MESH>::template AttributePtr<VEC> vertex_position, uint32 nb_vertices_to_remove)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	EdgeApproximator_MidEdge<MESH, VEC> approximator(m, vertex_position);
	EdgeTraversor_EdgeLength<MESH, VEC> traversor(m, vertex_position);

	uint32 count = 0;
	for (auto it = traversor.begin(); it != traversor.end(); ++it)
	{
		VEC newpos = approximator(*it);

		traversor.pre_collapse(*it);

		Vertex v = collapse_edge(m, *it);
		value<VEC>(m, vertex_position, v) = newpos;

		traversor.post_collapse();

		++count;
		if (count >= nb_vertices_to_remove)
			break;
	}
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_DECIMATION_H_
