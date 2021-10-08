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

#include <cgogn/modeling/algos/subdivision/surface_loop.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

using Scalar = geometry::Scalar;
using Vec3 = geometry::Vec3;

///////////
// CMap2 //
///////////

void subdivide_loop(CMap2& m, CMap2::Attribute<Vec3>* vertex_position)
{
	CellCache<CMap2> cache_old(m);
	cache_old.template build<CMap2::Vertex>();
	cache_old.template build<CMap2::Edge>();
	cache_old.template build<CMap2::Face>();

	CellCache<CMap2> cache_new(m);

	auto vertex_position_new = add_attribute<Vec3, CMap2::Vertex>(m, "__vertex_position_loop");

	foreach_cell(cache_old, [&](CMap2::Edge e) -> bool {
		cache_new.add(cut_edge(m, e));
		return true;
	});

	parallel_foreach_cell(cache_new, [&](CMap2::Vertex v) -> bool {
		CMap2::Vertex v1 = CMap2::Vertex(phi_1(m, v.dart));
		CMap2::Vertex v2 = CMap2::Vertex(phi1(m, v.dart));
		if (is_incident_to_boundary(m, v))
		{
			value<Vec3>(m, vertex_position_new, v) =
				0.5 * (value<Vec3>(m, vertex_position, v1) + value<Vec3>(m, vertex_position, v2));
		}
		else
		{
			CMap2::Vertex op1 = CMap2::Vertex(phi<111>(m, v.dart));
			CMap2::Vertex op2 = CMap2::Vertex(phi_1(m, phi_1(m, phi2(m, v.dart))));
			value<Vec3>(m, vertex_position_new, v) =
				3.0 / 8.0 * (value<Vec3>(m, vertex_position, v1) + value<Vec3>(m, vertex_position, v2)) +
				1.0 / 8.0 * (value<Vec3>(m, vertex_position, op1) + value<Vec3>(m, vertex_position, op2));
		}
		return true;
	});

	parallel_foreach_cell(cache_old, [&](CMap2::Vertex v) -> bool {
		Vec3 sum{0, 0, 0};

		Dart boundary;
		int nb_e = 0;
		foreach_incident_edge(m, v, [&](CMap2::Edge e) -> bool {
			nb_e++;
			sum += value<Vec3>(m, vertex_position, CMap2::Vertex(phi<11>(m, e.dart)));
			if (is_boundary(m, e.dart))
				boundary = e.dart;
			return true;
		});

		if (!boundary.is_nil()) // boundary case
		{
			CMap2::Vertex e1(phi<11>(m, boundary));
			CMap2::Vertex e2(phi_1(m, phi_1(m, boundary)));
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
	foreach_cell(cache_old, [&](CMap2::Face f) -> bool {
		Dart d0 = phi1(m, f.dart);
		Dart d1 = phi<11>(m, d0);
		cut_face(m, CMap2::Vertex(d0), CMap2::Vertex(d1));
		Dart d2 = phi<11>(m, d1);
		cut_face(m, CMap2::Vertex(d1), CMap2::Vertex(d2));
		Dart d3 = phi<11>(m, d2);
		cut_face(m, CMap2::Vertex(d2), CMap2::Vertex(d3));
		return true;
	});

	vertex_position->swap(vertex_position_new.get());
	remove_attribute<CMap2::Vertex>(m, vertex_position_new);
}

} // namespace modeling

} // namespace cgogn
