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

#ifndef CGOGN_GEOMETRY_ALGOS_VOLUME_H_
#define CGOGN_GEOMETRY_ALGOS_VOLUME_H_

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/face.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
Scalar volume(const MESH& m, typename mesh_traits<MESH>::Volume v,
			  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Face = typename mesh_traits<MESH>::Face;
	using Face2 = typename mesh_traits<MESH>::Face2;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	Scalar result = 0.0f;
	foreach_incident_face(m, v, [&](Face f) -> bool {
		Vec3 n = normal(m, Face2(f.dart_), vertex_position);
		Scalar a = area(m, f, vertex_position);
		result += value<Vec3>(m, vertex_position, Vertex(f.dart_)).dot(n) * a;
		return true;
	});
	result = fabs(result) / 3.0f;
	return result;
}

template <typename MESH>
Scalar volume(const MESH& m, typename mesh_traits<MESH>::Volume v,
			  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
			  const typename mesh_traits<MESH>::template Attribute<double>* area_face)
{
	using Face = typename mesh_traits<MESH>::Face;
	using Face2 = typename mesh_traits<MESH>::Face2;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	Scalar result = 0.0f;
	foreach_incident_face(m, v, [&](Face f) -> bool {
		Vec3 n = normal(m, Face2(f.dart), vertex_position);
		double a = value<double>(m, area_face, f);
		result += value<Vec3>(m, vertex_position, Vertex(f.dart)).dot(n) * a;
		return true;
	});
	result = fabs(result) / 3.0f;
	return result;
}

template <typename MESH>
void compute_volume(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
					typename mesh_traits<MESH>::template Attribute<double>* computed_volume)
{
	using Volume = typename mesh_traits<MESH>::Volume;

	foreach_cell(m, [&](Volume v) -> bool {
		value<double>(m, computed_volume, v) = volume(m, v, vertex_position);
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_VOLUME_H_
