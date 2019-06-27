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

#ifndef CGOGN_GEOMETRY_ALGOS_NORMAL_H_
#define CGOGN_GEOMETRY_ALGOS_NORMAL_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/vector_ops.h>
#include <cgogn/geometry/functions/normal.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
Vec3
normal(
	const MESH& m,
	typename mesh_traits<MESH>::Face f,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3> vertex_position
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, f);
	if (vertices.size() == 3)
	{
		Vec3 n = normal(
			value<Vec3>(m, vertex_position, vertices[0]),
			value<Vec3>(m, vertex_position, vertices[1]),
			value<Vec3>(m, vertex_position, vertices[2])
		);
		n.normalize();
		return n;
	}
	else
	{
		Vec3 n{0.0, 0.0, 0.0};
		for (uint32 i = 0; i < vertices.size() - 1; ++i)
		{
			const Vec3& p = value<Vec3>(m, vertex_position, vertices[i]);
			const Vec3& q = value<Vec3>(m, vertex_position, vertices[i+1]);
			n[0] += (p[1] - q[1]) * (p[2] + q[2]);
			n[1] += (p[2] - q[2]) * (p[0] + q[0]);
			n[2] += (p[0] - q[0]) * (p[1] + q[1]);
		}
		n.normalize();
		return n;
	}
}

template <typename MESH>
Vec3
normal(
	const MESH& m,
	typename mesh_traits<MESH>::Vertex v,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3> vertex_position
)
{
	using Face = typename mesh_traits<MESH>::Face;
	Vec3 n{0.0, 0.0, 0.0};
	foreach_incident_face(m, v, [&] (Face f) -> bool
	{
		n += normal(m, f, vertex_position);
		return true;
	});
	n.normalize();
	return n;
}

template <typename MESH>
void
compute_normal(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3> vertex_position,
	typename mesh_traits<MESH>::template AttributePtr<Vec3> vertex_normal
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	foreach_cell(m, [&] (Vertex v) -> bool
	{
		value<Vec3>(m, vertex_normal, v) = normal(m, v, vertex_position);
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_NORMAL_H_
