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

#ifndef CGOGN_GEOMETRY_ALGOS_GEODESIC_H_
#define CGOGN_GEOMETRY_ALGOS_GEODESIC_H_

#include <cgogn/core/functions/mesh_ops/face.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/inclusion.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/intrinsic_triangulation.h>

namespace cgogn
{

namespace geometry
{

/**
 * performs a flip out operation on a intrinsic triangulated mesh
 * given a flexible joint formed by halfedges a and b respectively incoming and outcoming of one vertex
 * @param intr an intrinsic triangulation, will be modified
 * @param a an incoming halfedge
 * @param b an outcoming halfedge
*/
template <typename MESH>
void flip_out(IntrinsicTriangulation& intr,
	 const typename mesh_traits<MESH>::HalfEdge a,
	 const typename mesh_traits<MESH>::HalfEdge b)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	//TODO
}

/**
 * compute geodesic path with flip out algo on an intrinsic triangulation
 * @param intr an intrinsic triangulation, will be modified according to the geodesic
 * @param path a connected path from the begin to the end of the vector, will be updated to a geodesic path
*/
template <typename MESH>
void geodesic_path(IntrinsicTriangulation& intr,
	 std::vector<const typename mesh_traits<MESH>::Vertex>& path)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	//TODO
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_EAR_TRIANGULATION_H_
