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

#ifndef CGOGN_GEOMETRY_ALGOS_AREA_H_
#define CGOGN_GEOMETRY_ALGOS_AREA_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/face.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/functions/area.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
Scalar area(const MESH& m, typename mesh_traits<MESH>::Face f,
			const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, f);
	Scalar face_area{0};
	for (uint32 i = 1, size = uint32(vertices.size()); i < size - 1; ++i)
	{
		face_area += area(value<Vec3>(m, vertex_position, vertices[0]), value<Vec3>(m, vertex_position, vertices[i]),
						  value<Vec3>(m, vertex_position, vertices[(i + 1) % size]));
	}
	return face_area;
}

template <typename MESH>
Scalar area(const MESH& m, typename mesh_traits<MESH>::Vertex v,
			const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Face = typename mesh_traits<MESH>::Face;
	Scalar vertex_area{0};
	foreach_incident_face(m, v, [&](Face iface) -> bool {
		vertex_area += area(m, iface, vertex_position) / 3.0;
		return true;
	});
	return vertex_area;
}

template <typename CELL, typename MESH>
void compute_area(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
				  typename mesh_traits<MESH>::template Attribute<Scalar>* cell_area)
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "MESH dimension should be >= 2");

	parallel_foreach_cell(m, [&](CELL c) -> bool {
		value<Scalar>(m, cell_area, c) = area(m, c, vertex_position);
		return true;
	});
}

template <typename MESH>
Scalar area(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Face = typename mesh_traits<MESH>::Face;

	Scalar result = 0;
	foreach_cell(m, [&](Face f) -> bool {
		result += area(m, f, vertex_position);
		return true;
	});
	return result;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_AREA_H_
