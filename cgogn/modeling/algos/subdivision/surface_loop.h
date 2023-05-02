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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_LOOP_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_LOOP_H_

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

using geometry::Vec3;
using geometry::Scalar;

///////////
// CMap2 //
///////////

//void subdivide_loop(CMap2& m, CMap2::Attribute<Vec3>* vertex_position);

template <typename MAP2, typename std::enable_if_t<std::is_convertible_v<MAP2&, MapBase&>>* = nullptr>
void subdivide_loop(MAP2& m, typename MAP2::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename MAP2::Vertex;
	using Edge = typename MAP2::Edge;
	using Face = typename MAP2::Face;

	CellCache<MAP2> cache_old(m);
	cache_old.template build<Vertex>();
	cache_old.template build<Edge>();
	cache_old.template build<Face>();

	CellCache<MAP2> cache_new(m);

	auto vertex_position_new = add_attribute<Vec3, Vertex>(m, "__vertex_position_loop");

	foreach_cell(cache_old, [&](Edge e) -> bool {
		cache_new.add(cut_edge(m, e));
		return true;
	});

	parallel_foreach_cell(cache_new, [&](Vertex v) -> bool {
		Vertex v1 = Vertex(phi_1(m, v.dart));
		Vertex v2 = Vertex(phi1(m, v.dart));
		if (is_incident_to_boundary(m, v))
		{
			value<Vec3>(m, vertex_position_new, v) =
				0.5 * (value<Vec3>(m, vertex_position, v1) + value<Vec3>(m, vertex_position, v2));
		}
		else
		{
			Vertex op1 = Vertex(phi<1, 1, 1>(m, v.dart));
			Vertex op2 = Vertex(phi<2, -1, -1>(m, v.dart));
			value<Vec3>(m, vertex_position_new, v) =
				3.0 / 8.0 * (value<Vec3>(m, vertex_position, v1) + value<Vec3>(m, vertex_position, v2)) +
				1.0 / 8.0 * (value<Vec3>(m, vertex_position, op1) + value<Vec3>(m, vertex_position, op2));
		}
		return true;
	});

	parallel_foreach_cell(cache_old, [&](Vertex v) -> bool {
		Vec3 sum{0, 0, 0};

		Dart boundary;
		uint32 nb_e = 0;
		foreach_incident_edge(m, v, [&](Edge e) -> bool {
			nb_e++;
			sum += value<Vec3>(m, vertex_position, CMap2::Vertex(phi<1, 1>(m, e.dart)));
			if (is_boundary(m, e.dart))
				boundary = e.dart;
			return true;
		});

		if (!boundary.is_nil()) // boundary case
		{
			Vertex e1(phi<1, 1>(m, boundary));
			Vertex e2(phi<-1, -1>(m, boundary));
			value<Vec3>(m, vertex_position_new, v) =
				Scalar(3.0 / 4.0) * value<Vec3>(m, vertex_position, v) +
				Scalar(1.0 / 8.0) * (value<Vec3>(m, vertex_position, e1) + value<Vec3>(m, vertex_position, e2));
		}
		else
		{
			Scalar beta = 3.0 / 16.0;
			if (nb_e > 3)
				beta = 3.0 / (8.0 * nb_e);
			value<Vec3>(m, vertex_position_new, v) =
				beta * sum + (1.0 - beta * nb_e) * value<Vec3>(m, vertex_position, v);
		}

		return true;
	});

	// add edges inside faces
	foreach_cell(cache_old, [&](Face f) -> bool {
		Dart d0 = phi1(m, f.dart);
		Dart d1 = phi<1, 1>(m, d0);
		cut_face(m, Vertex(d0), Vertex(d1));
		Dart d2 = phi<1, 1>(m, d1);
		cut_face(m, Vertex(d1), Vertex(d2));
		Dart d3 = phi<1, 1>(m, d2);
		cut_face(m, Vertex(d2), Vertex(d3));
		return true;
	});

	vertex_position->swap(vertex_position_new.get());
	remove_attribute<Vertex>(m, vertex_position_new);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_LOOP_H_
