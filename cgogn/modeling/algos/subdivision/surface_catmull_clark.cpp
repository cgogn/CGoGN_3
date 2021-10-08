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

#include <cgogn/modeling/algos/subdivision/surface_catmull_clark.h>
#include <cgogn/modeling/algos/subdivision_utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/algos/centroid.h>

namespace cgogn
{

namespace modeling
{

using Scalar = geometry::Scalar;
using Vec3 = geometry::Vec3;

///////////
// CMap2 //
///////////

void subdivide_catmull_clark(CMap2& m, CMap2::Attribute<Vec3>* vertex_position)
{
	CellCache<CMap2> cache_old(m);
	cache_old.template build<CMap2::Vertex>();
	cache_old.template build<CMap2::Edge>();
	cache_old.template build<CMap2::Face>();

	CellCache<CMap2> cache_new(m);

	foreach_cell(cache_old, [&](CMap2::Edge e) -> bool {
		std::vector<CMap2::Vertex> iv = incident_vertices(m, e);
		CMap2::Vertex v = cut_edge(m, e);
		cache_new.add(v);
		value<Vec3>(m, vertex_position, v) =
			0.5 * (value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1]));
		return true;
	});

	foreach_cell(cache_old, [&](CMap2::Face f) -> bool {
		Vec3 center = geometry::centroid<Vec3>(m, f, vertex_position);
		CMap2::Vertex v = quadrangulate_face(m, f);
		value<Vec3>(m, vertex_position, v) = center;
		return true;
	});

	parallel_foreach_cell(cache_new, [&](CMap2::Vertex v) -> bool {
		if (!is_incident_to_boundary(m, v))
		{
			CMap2::Vertex f1(phi_1(m, v.dart));
			CMap2::Vertex f2(phi<211>(m, v.dart));
			value<Vec3>(m, vertex_position, v) =
				0.25 * (2.0 * value<Vec3>(m, vertex_position, v) + value<Vec3>(m, vertex_position, f1) +
						value<Vec3>(m, vertex_position, f2));
		}
		return true;
	});

	parallel_foreach_cell(cache_old, [&](CMap2::Vertex v) -> bool {
		Vec3 sum_F{0, 0, 0};
		Vec3 sum_E{0, 0, 0};

		Dart boundary;
		uint32 nb_f = 0;
		uint32 nb_e = 0;

		foreach_incident_edge(m, v, [&](CMap2::Edge e) -> bool {
			++nb_e;
			sum_E += 0.5 * (value<Vec3>(m, vertex_position, v) +
							value<Vec3>(m, vertex_position, CMap2::Vertex(phi<1211>(m, e.dart))));
			if (!is_boundary(m, e.dart))
			{
				++nb_f;
				sum_F += value<Vec3>(m, vertex_position, CMap2::Vertex(phi<11>(m, e.dart)));
			}
			else
				boundary = e.dart;
			return true;
		});

		sum_F /= nb_f;
		sum_E /= nb_e;

		if (!boundary.is_nil()) // boundary case
		{
			CMap2::Vertex e1(phi<11>(m, boundary));
			CMap2::Vertex e2(phi_1(m, phi_1(m, boundary)));
			value<Vec3>(m, vertex_position, v) =
				3.0 / 4.0 * value<Vec3>(m, vertex_position, v) +
				1.0 / 8.0 * (value<Vec3>(m, vertex_position, e1) + value<Vec3>(m, vertex_position, e2));
		}
		else
		{
			value<Vec3>(m, vertex_position, v) =
				1.0 / nb_e * (sum_F + 2.0 * sum_E + (nb_e - 3.0) * value<Vec3>(m, vertex_position, v));
		}

		return true;
	});
}

} // namespace modeling

} // namespace cgogn
