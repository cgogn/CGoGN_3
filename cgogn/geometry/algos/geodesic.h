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

using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
using Vertex = typename mesh_traits<CMap2>::Vertex;
using Edge = typename mesh_traits<CMap2>::Edge;

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
	std::list<Edge>::iterator edge_ref;
	Scalar angle;
	bool inverted = false;

	/**
	 * construct a joint A -> B from 2 consecutive edges
	 * @param intr intrinsic triangulation
	 * @param edge_iterator a list iterator pointing on an edge of a path, it has to be followed by at least one element
	*/
	Joint (IntrinsicTriangulation& intr,
		std::list<Edge>::iterator edge_iterator)
	{
		CMap2& m = intr.getMesh();
		edge_ref = edge_iterator;	// keep a ref to edge A
		a = (*edge_iterator).dart;
		b = (*(++edge_iterator)).dart;

		// reorder as a -> b
		if (isSameVertex(m, a, b)) {
			a = phi2(m, a);
		} else if (isSameVertex(m, phi2(m, b), a)) {
			a = phi2(m, a);
			b = phi2(m, b);
		} else if (isSameVertex(m, phi2(m, a), phi2(m, b))) {
			b = phi2(m, b);
		} else {
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

inline bool isInPath(const CMap2& mesh,
	std::list<Edge>& path,
	const Dart& d)
{
	for (auto e : path)
	{
		if (d == e.dart || d == phi2(mesh, e.dart))
			return true;
	}
	return false;
}

inline bool isFlexible(IntrinsicTriangulation& intr,
	std::list<Edge>& path,
	const Joint& joint)
{
	const CMap2& mesh = intr.getMesh();

	if (isInPath(mesh, path, phi <2, 1>(mesh, joint.a)))
		return false;
	for(Dart d = phi<2, -1, 2>(mesh, joint.a); d != joint.b; d = phi<-1, 2>(mesh, d))
	{
		if (isInPath(mesh, path, d))
			return false;
		if (isInPath(mesh, path, phi1(mesh, d)))
			return false;
	}

	for(Dart d = phi<2, -1, 2>(mesh, joint.a); d != joint.b; d = phi<-1, 2>(mesh, d))
	{
		if (edge_can_flip(intr.getMesh(), Edge(d)))	// TODO better and check if wedge is crossing path
			return true;
	}
	return  false;
}

inline std::list<Joint> buildJointList(IntrinsicTriangulation& intr,
	std::list<Edge>& path)
{
	std::list<Joint> joints;
	std::list<Edge>::iterator end_iter = end(path);
	end_iter--;
	for (std::list<Edge>::iterator e = begin(path); e != end_iter; e++)
	{
		Joint j {intr, e};
		if (isFlexible(intr, path, j))
			joints.push_back(j);
	}
	joints.sort([] (const Joint& a, const Joint& b) -> bool {
		return a.angle > b.angle;
	});
	return joints;
}

inline auto indexOfBi(IntrinsicTriangulation& intr, const std::list<Dart>& wedge, const Joint& j)
{
	return std::find_if(begin(wedge), end(wedge), [&](const Dart& d) -> bool {
		Dart a = phi1(intr.getMesh(), d);
		Dart b = phi<2, -1, 2>(intr.getMesh(), d);
		return intr.interiorAngle(a, b) < M_PI;
	});
}

/**
 * performs a flip out operation on a intrinsic triangulated mesh
 * given a flexible joint formed by halfedges a and b respectively incoming and outcoming of one vertex
 * @param intr an intrinsic triangulation, will be modified
 * @param j a flexible joint to shorten
 * @returns the shortest path between a and phi1(b) vertices
*/
std::list<Edge> flip_out(
	IntrinsicTriangulation& intr,
	const Joint& j)
{
	CMap2& mesh = intr.getMesh();
	assert(isSameVertex(mesh, phi1(mesh, j.a), j.b));

	std::list<Dart> wedge;	// contains the adjacent edges from a to b counterclockwise
	for(Dart d = phi<2, -1, 2>(mesh, j.a); d != j.b; d = phi<-1, 2>(mesh, d))
		wedge.push_back(d);

	bool update;
	do {
		update = false;
		auto index = indexOfBi(intr, wedge, j);
		while (index != end(wedge))
		{
			if(intr.flip_edge(Edge(*index)))
				update = true;
			wedge.erase(index);
			index = indexOfBi(intr, wedge, j);
		}
	} while (update);	// update until stability

	// build the shorter path
	std::list<Edge> shorter;
	shorter.push_back(Edge(phi<2, 1>(mesh, j.a)));
	for(Dart d : wedge)
		shorter.push_back(Edge(phi1(mesh, d)));
	
	return shorter;
}

} // end helper

/**
 * compute geodesic path with flip out algorithm on an intrinsic triangulation
 * @param intr an intrinsic triangulation, will be modified according to the geodesic
 * @param path an intrinsic connected path from the begin to the end of the vector, will be updated toward a geodesic
 * @param iteration optional, the maximum number of flip out to do, a negative value will loop until stability
*/
void geodesic_path(IntrinsicTriangulation& intr,
	 std::list<Edge>& path,
	 int iteration = -1)
{
	const CMap2& mesh = intr.getMesh();

	auto joints_priority_queue = buildJointList(intr, path);
	while (joints_priority_queue.size() > 0 && iteration != 0)
	{
		Joint j = joints_priority_queue.back();
		auto shorter_path = flip_out(intr, j);

		if (j.inverted)
		{
			shorter_path.reverse();
			for(auto e : shorter_path)
				e = Edge(phi2(mesh, e.dart));
		}
		// update
		std::list<Edge>::iterator it = j.edge_ref;
		++it;
		it = path.erase(j.edge_ref, ++it);
		path.splice(it, shorter_path);

		joints_priority_queue = buildJointList(intr, path);	// TODO optimizable
		--iteration;
	}
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_GEODESIC_H_
