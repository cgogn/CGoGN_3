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

// computes the gradient of a vertex scalar field at the faces of a mesh
template <typename MESH>
void compute_gradient_of_vertex_scalar_field(
	const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_scalar_field,
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
			g += value<Scalar>(m, vertex_scalar_field, vertices[i]) * n.cross(e);
		}
		g /= 2 * a;
		value<Vec3>(m, face_gradient, f) = -1.0 * g.normalized();
		return true;
	});
}

// computes the divergence of a face vector field at the vertices of a mesh
template <typename MESH>
void compute_div_of_face_vector_field(const MESH& m,
									  const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
									  const typename mesh_traits<MESH>::template Attribute<Vec3>* face_vector_field,
									  typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_divergence)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		Scalar d = face_vector_field_divergence(m, v, face_vector_field, vertex_position);
		value<Scalar>(m, vertex_divergence, v) = d;
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_GRADIENT_H_
