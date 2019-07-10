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

template <typename MESH>
Scalar
angle(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3>& vertex_position
)
{
	using Face = typename mesh_traits<MESH>::Face;
	std::vector<Face> faces = incident_faces(m, e);
	if (faces.size() < 2)
		return 0;
	return angle(
        normal(m, faces[0], vertex_position),
        normal(m, faces[1], vertex_position)
    );
}

template <typename MESH>
Scalar
angle(
	const MESH& m,
	typename mesh_traits<MESH>::Edge e,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3>& vertex_position,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3>& face_normal
)
{
	using Face = typename mesh_traits<MESH>::Face;
	std::vector<Face> faces = incident_faces(m, e);
	return angle(
        value<Vec3>(m, faces[0], face_normal),
        value<Vec3>(m, faces[1], face_normal)
    );
}

template <typename MESH>
void
compute_angle(
	const MESH& m,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3>& vertex_position,
	const typename mesh_traits<MESH>::template AttributePtr<Scalar>& edge_angle
)
{
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
	const typename mesh_traits<MESH>::template AttributePtr<Vec3>& vertex_position,
	const typename mesh_traits<MESH>::template AttributePtr<Vec3>& face_normal,
	const typename mesh_traits<MESH>::template AttributePtr<Scalar>& edge_angle
)
{
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
