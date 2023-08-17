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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_CATMULL_CLARK_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_CATMULL_CLARK_H_

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/modeling/algos/subdivision_utils.h>

#include <cgogn/core/types/mesh_views/cell_cache.h>

namespace cgogn
{

struct MapBase;

namespace modeling
{

using geometry::Scalar;
using geometry::Vec3;

///////////////
// MapBase:2 //
///////////////

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&> &&
												   (mesh_traits<MESH>::dimension == 2)>* = nullptr>
void subdivide_catmull_clark(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	CellCache<MESH> cache_old(m);
	cache_old.template build<Vertex>();
	cache_old.template build<Edge>();
	cache_old.template build<Face>();

	CellCache<MESH> cache_new(m);

	foreach_cell(cache_old, [&](Edge e) -> bool {
		std::vector<Vertex> iv = incident_vertices(m, e);
		Vertex v = cut_edge(m, e);
		cache_new.add(v);
		value<Vec3>(m, vertex_position, v) =
			0.5 * (value<Vec3>(m, vertex_position, iv[0]) + value<Vec3>(m, vertex_position, iv[1]));
		return true;
	});

	foreach_cell(cache_old, [&](Face f) -> bool {
		Vec3 center = geometry::centroid<Vec3>(m, f, vertex_position);
		Vertex v = quadrangulate_face(m, f);
		value<Vec3>(m, vertex_position, v) = center;
		return true;
	});

	parallel_foreach_cell(cache_new, [&](Vertex v) -> bool {
		if (!is_incident_to_boundary(m, v))
		{
			Vertex f1(phi_1(m, v.dart_));
			Vertex f2(phi<2, 1, 1>(m, v.dart_));
			value<Vec3>(m, vertex_position, v) =
				0.25 * (2.0 * value<Vec3>(m, vertex_position, v) + value<Vec3>(m, vertex_position, f1) +
						value<Vec3>(m, vertex_position, f2));
		}
		return true;
	});

	parallel_foreach_cell(cache_old, [&](Vertex v) -> bool {
		Vec3 sum_F{0, 0, 0};
		Vec3 sum_E{0, 0, 0};

		Dart boundary;
		uint32 nb_f = 0;
		uint32 nb_e = 0;

		foreach_incident_edge(m, v, [&](Edge e) -> bool {
			++nb_e;
			sum_E += 0.5 * (value<Vec3>(m, vertex_position, v) +
							value<Vec3>(m, vertex_position, Vertex(phi<1, 2, 1, 1>(m, e.dart_))));
			if (!is_boundary(m, e.dart_))
			{
				++nb_f;
				sum_F += value<Vec3>(m, vertex_position, Vertex(phi<1, 1>(m, e.dart_)));
			}
			else
				boundary = e.dart_;
			return true;
		});

		sum_F /= nb_f;
		sum_E /= nb_e;

		if (!boundary.is_nil()) // boundary case
		{
			Vertex e1(phi<1, 1>(m, boundary));
			Vertex e2(phi<-1, -1>(m, boundary));
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

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_CATMULL_CLARK_H_
