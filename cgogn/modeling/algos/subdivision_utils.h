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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_UTILS_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_UTILS_H_

#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{
//forward for SFINAE
struct CMap2;

namespace modeling
{

using Vec3 = geometry::Vec3;

///////////
// CMap2 //
///////////

template <typename MESH, typename std::enable_if_t<std::is_same_v<MESH&, CMap2&>>* = nullptr>
void hexagon_to_triangles(MESH& m, typename mesh_traits<MESH>::Face f)
{
	using Vertex = typename mesh_traits<MESH>:Vertex;
	cgogn_message_assert(codegree(m, f) == 6, "hexagon_to_triangles: given face should have 6 edges");
	Dart d0 = phi1(m, f.dart);
	Dart d1 = phi<1, 1>(m, d0);
	cut_face(m, Vertex(d0), Vertex(d1));
	Dart d2 = phi<1, 1>(m, d1);
	cut_face(m, Vertex(d1), Vertex(d2));
	Dart d3 = phi<1, 1>(m, d2);
	cut_face(m, CMap2::Vertex(d2), CMap2::Vertex(d3));
}

//////////////
// CMapBase //
//////////////

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>* = nullptr>
typename mesh_traits<MESH>::Vertex quadrangulate_face(MESH& m, typename mesh_traits<MESH>::Face f)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	cgogn_message_assert(codegree(m, f) % 2 == 0, "quadrangulate_face: given face should have a pair codegree");

	Dart d0 = phi1(m, f.dart);
	Dart d1 = phi<1, 1>(m, d0);

	cut_face(m, Vertex(d0), Vertex(d1));
	cut_edge(m, Edge(phi_1(m, d0)));

	Dart x = phi<-1, 2>(m, d0);
	Dart dd = phi<1, 1, 1, 1>(m, x);
	while (dd != x)
	{
		Dart next = phi<1, 1>(m, dd);
		cut_face(m, Vertex(dd), Vertex(phi1(m, x)));
		dd = next;
	}

	return Vertex(phi2(m, x));
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_UTILS_H_
