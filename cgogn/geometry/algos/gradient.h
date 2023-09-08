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

#ifndef CGOGN_GEOMETRY_ALGOS_GRADIENT_H_
#define CGOGN_GEOMETRY_ALGOS_GRADIENT_H_

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

/////////////
// GENERIC //
/////////////

template <typename MESH>
void compute_gradient_of_vertex_scalar_field(
	const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_function,
	typename mesh_traits<MESH>::template Attribute<Vec3>* face_gradient)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	parallel_foreach_cell(m, [&](Face f) -> bool {
		Vec3 g(0, 0, 0);
		Vec3 n = geometry::normal(m, f, vertex_position);
		Scalar a = geometry::area(m, f, vertex_position);
		std::vector<Vertex> vertices = incident_vertices(m, f);
		for (uint32 i = 0; i < vertices.size(); ++i)
		{
			Vec3 e = value<Vec3>(m, vertex_position, vertices[(i + 2) % vertices.size()]) -
					 value<Vec3>(m, vertex_position, vertices[(i + 1) % vertices.size()]);
			g += value<Scalar>(m, vertex_function, vertices[i]) * n.cross(e);
		}
		g /= 2 * a;
		value<Vec3>(m, face_gradient, f) = -1.0 * g.normalized();
		return true;
	});
}

///////////
// CMap2 //
///////////

inline Scalar vertex_gradient_divergence(const CMap2& m, CMap2::Vertex v, const CMap2::Attribute<Vec3>* face_gradient,
										 const CMap2::Attribute<Vec3>* vertex_position)
{
	Scalar div = 0.0;
	std::vector<CMap2::Edge> edges = incident_edges(m, v);
	for (uint32 i = 0; i < edges.size(); ++i)
	{
		CMap2::Edge e1 = edges[i];
		CMap2::Edge e2 = edges[(i + 1) % edges.size()];

		CMap2::Face f(e1.dart_);

		const Vec3& X = value<Vec3>(m, face_gradient, f);

		const Vec3& p0 = value<Vec3>(m, vertex_position, v);
		const Vec3& p1 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e1.dart_)));
		const Vec3& p2 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e2.dart_)));

		Vec3 vecR = p0 - p2;
		Vec3 vecL = p1 - p2;
		Scalar cotValue1 = vecR.dot(vecL) / vecR.cross(vecL).norm();

		vecR = p2 - p1;
		vecL = p0 - p1;
		Scalar cotValue2 = vecR.dot(vecL) / vecR.cross(vecL).norm();

		div += cotValue1 * (p1 - p0).dot(X) + cotValue2 * (p2 - p0).dot(X);
	}
	return div / 2.0;
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void compute_div_of_face_vector_field(const MESH& m,
									  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
									  const typename mesh_traits<MESH>::template Attribute<Vec3>* face_vector_field,
									  typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_divergence)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		Scalar d = vertex_gradient_divergence(m, v, face_vector_field, vertex_position);
		value<Scalar>(m, vertex_divergence, v) = d;
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_GRADIENT_H_
