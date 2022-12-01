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

namespace // helper
{

inline auto indexOfBi(IntrinsicTriangulation& intr, std::vector<Dart> wedge)
{
    return std::find_if(begin(wedge)+1, end(wedge)-1, [&](Dart d) -> bool {
		Dart n = phi2(intr.getMesh(), d);
		Scalar a = intr.getAngle(phi<2, 1>(intr.getMesh(), n));
		Scalar b = intr.getAngle(phi<-1, 2>(intr.getMesh(), n));
		if (a>b)
			return 2*M_PI - abs(a-b) < M_PI;
		return abs(a-b) < M_PI;
		
	});
}

} // end helper

/**
 * performs a flip out operation on a intrinsic triangulated mesh
 * given a flexible joint formed by halfedges a and b respectively incoming and outcoming of one vertex
 * @param intr an intrinsic triangulation, will be modified
 * @param a the incoming halfedge
 * @param b the outcoming halfedge
 * @returns a shorter path between a and phi1(b)
*/
std::vector<typename mesh_traits<CMap2>::Edge> flip_out(
	IntrinsicTriangulation& intr,
	const Dart& a,
	const Dart& b)
{
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
	using Edge = typename mesh_traits<CMap2>::Edge;
	using Face = typename mesh_traits<CMap2>::Face;

	CMap2& mesh = intr.getMesh();
	std::vector<Dart> wedge;
	for(Dart d = phi2(mesh, a); d != b; d = phi<-1, 2>(mesh, d))
		wedge.push_back(d);
	wedge.push_back(b);

	auto index=indexOfBi(intr, wedge);
	while (*index != wedge.back())
	{
		intr.flip_edge(Edge(*index));
		wedge.erase(index);
		index=indexOfBi(intr, wedge);
	}

	// build shorter path from simplified wedge
	wedge.pop_back();
	std::vector<Edge> shorter;
	for(Dart d : wedge) {
		shorter.push_back(Edge(phi<2, -1>(mesh, d)));
	}
	return shorter;
}

/**
 * compute geodesic path with flip out algo on an intrinsic triangulation
 * @param intr an intrinsic triangulation, will be modified according to the geodesic
 * @param path a connected path from the begin to the end of the vector, will be updated to a geodesic path
*/
void geodesic_path(IntrinsicTriangulation& intr,
	 std::vector<typename mesh_traits<CMap2>::Edge>& path)
{
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
	using Edge = typename mesh_traits<CMap2>::Edge;
	using Face = typename mesh_traits<CMap2>::Face;
	//TODO
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_EAR_TRIANGULATION_H_
