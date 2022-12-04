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

inline bool isSameVertex(const CMap2& m, Dart a, Dart b)
{
	bool r = false;
	foreach_dart_of_PHI21(m, a, [&] (Dart d) -> bool {
		if (d == b) {
			r = true;
			return false;
		}
		return true;
	});
	return r;
}

struct Joint
{
	Dart a;
	Dart b;
	std::list<typename mesh_traits<CMap2>::Edge>::iterator edge_ref;
	Scalar angle;
	bool inverted = false;

	/**
	 * construct a joint A -> B from 2 adjacent edges
	 * @param intr intrinsic triangulation
	 * @param edge_iterator a list iterator pointing on an edge of a path, it have to be followed by at least one element
	*/
	Joint (IntrinsicTriangulation& intr,
		std::list<typename mesh_traits<CMap2>::Edge>::iterator edge_iterator)
	{
		using Vertex = typename mesh_traits<CMap2>::Vertex;
		CMap2& m = intr.getMesh();
		edge_ref = edge_iterator;	// keep a ref to edge A
		a = (*edge_iterator).dart;
		b = (*(++edge_iterator)).dart;

		// reorder as a -> b
		if (isSameVertex(m, a, b)) {
			a = phi2(m, a);
		}
		else if (isSameVertex(m, phi2(m, b), a)) {
			a = phi2(m, a);
			b = phi2(m, b);
		}
		else if (isSameVertex(m, phi2(m, a), phi2(m, b))) {
			b = phi2(m, b);
		}
		else {
			assert(isSameVertex(m, phi2(m, a), b)); // a_edge and b_edge are not adjacent
		}
		
		angle = intr.interiorAngle(phi2(m, a), b);
		
		if (angle > M_PI)
		{
			// inverse to get the lowest angle
			Dart t = phi2(m, b);
			b = phi2(m, a);
			a = t;
			angle = 2 * M_PI - angle;
			inverted = true;
		}
	}
};

inline bool isFlexible(IntrinsicTriangulation& intr,
	std::list<typename mesh_traits<CMap2>::Edge>& path,
	const Joint& joint)
{
	using Edge = typename mesh_traits<CMap2>::Edge;

	const CMap2& mesh = intr.getMesh();
	for(Dart d = phi<2, -1, 2>(mesh, joint.a); d != joint.b; d = phi<-1, 2>(mesh, d))
	{
		if (edge_can_flip(intr.getMesh(), Edge(d)))	// TODO check if wedge is crossing path
			return true;
	}
	
	return  false;
}

inline auto buildJointList(IntrinsicTriangulation& intr,
	std::list<typename mesh_traits<CMap2>::Edge>& path)
{
	std::list<Joint> joints;
	auto end_iter = end(path);
	end_iter--;
	for (auto e = begin(path); e != end_iter; e++)
	{
		Joint j {intr, e};
		if (isFlexible(intr, path, j))
			joints.push_back(j);
	}
	// TODO sort joints by argmin angle
	return joints;
}

inline auto indexOfBi(IntrinsicTriangulation& intr, const std::list<Dart>& wedge)
{
	return std::find_if(begin(wedge), end(wedge), [&](const Dart& d) -> bool {
		Dart a = phi1(intr.getMesh(), d);
		Dart b = phi<2, -1, 2>(intr.getMesh(), d);
		return intr.interiorAngle(a, b) < M_PI;
	});
}

} // end helper

/**
 * performs a flip out operation on a intrinsic triangulated mesh
 * given a flexible joint formed by halfedges a and b respectively incoming and outcoming of one vertex
 * @param intr an intrinsic triangulation, will be modified
 * @param a the incoming halfedge
 * @param b the outcoming halfedge
 * @returns a shorter path between a and phi1(b) vertices
*/
std::list<typename mesh_traits<CMap2>::Edge> flip_out(
	IntrinsicTriangulation& intr,
	Dart a,
	Dart b)
{
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using Edge = typename mesh_traits<CMap2>::Edge;
	
	CMap2& mesh = intr.getMesh();
	assert(isSameVertex(mesh, phi1(mesh, a), b));

	std::list<Dart> wedge;	// contains the adjacent edges of a and b on the side with the smaller angle
	for(Dart d = phi<2, -1, 2>(mesh, a); d != b; d = phi<-1, 2>(mesh, d))
		wedge.push_back(d);

	auto index = indexOfBi(intr, wedge);
	bool update = true;
	while (update)	// update until stability
	{
		update = false;
		while (index != end(wedge))
		{
			auto e = end(wedge);
			if(intr.flip_edge(Edge(*index)))
				update = true;
			wedge.erase(index);
			index = indexOfBi(intr, wedge);
		}
	}

	// build shorter path from simplified wedge
	std::list<Edge> shorter;
	shorter.push_back(Edge(phi<2, 1>(mesh, a)));
	for(Dart d : wedge) {
		shorter.push_back(Edge(phi1(mesh, d)));
	}
	return shorter;
}

/**
 * compute geodesic path with flip out algo on an intrinsic triangulation
 * @param intr an intrinsic triangulation, will be modified according to the geodesic
 * @param path a connected path from the begin to the end of the vector, will be updated to a geodesic path
 * @param iteration optional, the maximum number of flip out to do, -1 mean to iterate until the path is the shortest path
*/
void geodesic_path(IntrinsicTriangulation& intr,
	 std::list<typename mesh_traits<CMap2>::Edge>& path,
	 int iteration = -1)
{
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
	using Edge = typename mesh_traits<CMap2>::Edge;

	assert (iteration >= -1);
	const CMap2& mesh = intr.getMesh();

	auto joints_priority_queue = buildJointList(intr, path);
	while (joints_priority_queue.size() > 0 && iteration != 0)
	{
		Joint j = joints_priority_queue.back();
		auto shorter_path = flip_out(intr, j.a, j.b);

		// update
		if (j.inverted)
			shorter_path.reverse();
		std::list<Edge>::iterator it = j.edge_ref;
		++it;
		it = path.erase(j.edge_ref, ++it);
		path.splice(it, shorter_path);

		joints_priority_queue = buildJointList(intr, path);	// optimizable
		--iteration;
	}
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_EAR_TRIANGULATION_H_
