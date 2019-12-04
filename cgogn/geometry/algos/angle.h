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

#ifndef CGOGN_GEOMETRY_ALGOS_ANGLE_H_
#define CGOGN_GEOMETRY_ALGOS_ANGLE_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

std::vector<Scalar>
opposite_angles(
	const CMap2& m,
	typename CMap2::Edge e,
	const typename mesh_traits<CMap2>::template Attribute<Vec3>* vertex_position
)
{
	if (!is_incident_to_boundary(m, e))
	{
		const Vec3& p1 = value<Vec3>(m, vertex_position, CMap2::Vertex(e.dart));
		const Vec3& p2 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m,e.dart)));
		const Vec3& p3 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi_1(m,e.dart)));
		const Vec3& p4 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi_1(m,phi2(m,e.dart))));
		return { angle(p1 - p3, p2 - p3), angle(p2 - p4, p1 - p4) };
	}
	else
	{
		const Vec3& p1 = value<Vec3>(m, vertex_position, CMap2::Vertex(e.dart));
		const Vec3& p2 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m,e.dart)));
		const Vec3& p3 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi_1(m,e.dart)));
		return { angle(p1 - p3, p2 - p3) };
	}
}

template <typename MESH>
Scalar
angle(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position
)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	std::vector<Face> faces = incident_faces(m, e);
	if (faces.size() < 2)
		return 0;

	const Vec3 n1 = normal(m, faces[0], vertex_position);
	const Vec3 n2 = normal(m, faces[1], vertex_position);

	std::vector<Vertex> vertices = incident_vertices(m, e);
	Vec3 edge = value<Vec3>(m, vertex_position, vertices[1]) - value<Vec3>(m, vertex_position, vertices[0]);
	edge.normalize();
	Scalar s = edge.dot(n1.cross(n2));
	Scalar c = n1.dot(n2);
	Scalar a(0);

	// the following trick is useful to avoid NaNs (due to floating point errors)
	if (c > Scalar(0.5)) a = std::asin(s);
	else
	{
		if (c < -1) c = -1;
		if (s >= 0) a = std::acos(c);
		else a = -std::acos(c);
	}

	return a;
	
	// return angle(
	//	 normal(m, faces[0], vertex_position),
	//	 normal(m, faces[1], vertex_position)
	// );
}

template <typename MESH>
Scalar
angle(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* face_normal
)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Face = typename mesh_traits<MESH>::Face;
	std::vector<Face> faces = incident_faces(m, e);
	if (faces.size() < 2)
		return 0;
	return angle(
		value<Vec3>(m, face_normal, faces[0]),
		value<Vec3>(m, face_normal, faces[1])
	);
}

template <typename MESH>
void
compute_angle(
	const MESH& m,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	typename mesh_traits<MESH>::template Attribute<Scalar>* edge_angle
)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Edge = typename mesh_traits<MESH>::Edge;

	parallel_foreach_cell(m, [&] (Edge e) -> bool
	{
		value<Scalar>(m, edge_angle, e) = angle(m, e, vertex_position);
		return true;
	});
}

template <typename MESH>
void
compute_angle(
	const MESH& m,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* face_normal,
	typename mesh_traits<MESH>::template Attribute<Scalar>* edge_angle
)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Edge = typename mesh_traits<MESH>::Edge;
	parallel_foreach_cell(m, [&] (Edge e) -> bool
	{
		value<Scalar>(m, edge_angle, e) = angle(m, e, vertex_position, face_normal);
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_ANGLE_H_
