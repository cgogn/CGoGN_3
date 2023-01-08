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

/**
 * algorithm from Sharp & Crane
 * "You can find geodesic paths in triangle meshes by just flipping edges"
*/

#ifndef CGOGN_GEOMETRY_ALGOS_GEODESIC_H_
#define CGOGN_GEOMETRY_ALGOS_GEODESIC_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/intrinsic_triangulation.h>
#include <cgogn/core/functions/traversals/vertex.h>

namespace cgogn
{

namespace geometry
{

namespace // helper
{

template <typename T>
using Attribute = typename mesh_traits<CMap2>::template Attribute<T>;
using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
using Vertex = typename mesh_traits<CMap2>::Vertex;
using Edge = typename mesh_traits<CMap2>::Edge;

constexpr Scalar epsilon = 0.00001;

// check if the dart a and b are on the same vertex
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
	Dart a;	// incoming halfedge
	Dart b;	// outcoming halfedge
	std::list<Edge>::iterator edge_ref_a;	// reference of the given edge a
	std::list<Edge>::iterator edge_ref_b;	// reference of the given edge b
	Scalar angle;		// extrinsic interior angle between a and b
	bool inverted = false;	// true if a and b has been inverted
	bool loop = false;	// true if a and b are disjoint so edge_ref_a is back and edge_ref_b is begin

	/**
	 * construct a joint A -> B from 2 consecutive edges
	 * @param intr intrinsic triangulation
	 * @param edge_iterator_a a list iterator pointing on an edge of a path
	 * @param edge_iterator_b a list iterator pointing on an edge of a path
	*/
	Joint (IntrinsicTriangulation& intr,
		std::list<Edge>::iterator edge_iterator_a,
		std::list<Edge>::iterator edge_iterator_b)
	{
		CMap2& m = intr.getMesh();
		edge_ref_a = edge_iterator_a;	// keep a ref to edge A
		edge_ref_b = edge_iterator_b;	// keep a ref to edge B
		a = (*edge_iterator_a).dart;
		b = (*edge_iterator_b).dart;

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
		
		angle = intr.intrinsicInteriorAngle(phi2(m, a), b);
		
		if (angle > M_PI)
		{
			// inverse to get a lower angle
			angle = intr.intrinsicInteriorAngle(b, phi2(m, a));
			Dart t = phi2(m, b);
			b = phi2(m, a);
			a = t;
			inverted = true;
		}
	}
};

// the priority queue of joints to flip out sorted by minimum angle
struct JointPriorityComparer {
	typedef Joint const& param_type;
	bool operator()(param_type a, param_type b) const {
		return a.angle > b.angle;
	};
};
using JointsPriorityQueue = std::priority_queue<Joint, std::vector<Joint>, JointPriorityComparer>;

/**
 * check if dart is in a list of edge
 * @param mesh a CMap2 data structure
 * @param path a list of connected edges
 * @param d the dart to check
 * @returns true if the dart is contained in any edge of the list
*/
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

/**
 * check if a joint is flexible, ie its wedges does not contain any path's edge
 * check also for a case where any wedge cannot flip and the inner and outer joints form a convex diamond
 * this last case will issue in the non convergence of the algorithm as it will loop over the two joints
 * @param intr an intrinsic triangulation
 * @param path a list of connected edges
 * @param joint the joint to check
 * @returns true if the joint can be used in flip out
*/
inline bool isFlexible(IntrinsicTriangulation& intr,
	std::list<Edge>& path,
	const Joint& joint)
{
	const CMap2& mesh = intr.getMesh();
	bool no_wedge_can_flip = true;

	if (isInPath(mesh, path, phi <2, 1>(mesh, joint.a)))
		return false;
	// iterate over wedge segments
	for(Dart d = phi<2, -1, 2>(mesh, joint.a); d != joint.b; d = phi<-1, 2>(mesh, d))
	{
		if (isInPath(mesh, path, d))
			return false;
		if (isInPath(mesh, path, phi1(mesh, d)))
			return false;
		if (no_wedge_can_flip && intr.can_be_flipped(Edge(d)))
			no_wedge_can_flip = false;
	}

	// check if the outer wedge is longer, if so return false
	if (no_wedge_can_flip)
	{
		Scalar old_length = intr.getLength(joint.a) + intr.getLength(joint.b);

		Scalar new_length = intr.getLength(phi<2, 1>(mesh, joint.a));
		for(Dart d = phi<2, -1, 2>(mesh, joint.a); d != joint.b; d = phi<-1, 2>(mesh, d))
			new_length += intr.getLength(phi1(mesh, d));
		
		if (old_length <= new_length)
			return false;
	}

	return true;
}

/**
 * construct a priority queue of joints that can be flipped out
 * @param intr an intrinsic triangulation
 * @param path a list of connected edges
 * @param loop true if the last nd first edge are connected
 * @returns queue of joints to flip out
*/
inline JointsPriorityQueue buildJointList(IntrinsicTriangulation& intr,
	std::list<Edge>& path,
	bool loop = false)
{
	JointsPriorityQueue joints;
	std::list<Edge>::iterator end_iter = end(path);
	end_iter--;
	for (std::list<Edge>::iterator e = begin(path); e != end_iter; e++)
	{
		std::list<Edge>::iterator b = e;
		++b;
		Joint j {intr, e, b};
		if (j.angle < M_PI && isFlexible(intr, path, j))
			joints.push(j);
	}
	if (loop)
	{
		Joint j {intr, end_iter, begin(path)};
		j.loop = true;
		if (j.angle < M_PI && isFlexible(intr, path, j))
			joints.push(j);
	}
	return joints;
}

/**
 * search the first wedge segment that can be flipped
 * @param intr an intrinsic triangulation
 * @param wedge a list of consecutive adjacent halfedges
 * @returns an iterator to the founded wedge or end(wedge) if not found
*/
inline std::list<Dart>::const_iterator indexOfBi(IntrinsicTriangulation& intr, const std::list<Dart>& wedge)
{
	auto chosen = end(wedge);
	auto it = begin(wedge);
	Scalar min_angle = M_PI;	// upper bound

	while (it != end(wedge))
	{
		it = std::find_if(begin(wedge), end(wedge), [&](const Dart& d) -> bool {
			if (!intr.can_be_flipped(Edge(d)))
				return false;

			Dart a = phi1(intr.getMesh(), d);
			Dart b = phi<2, -1, 2>(intr.getMesh(), d);
			Scalar angle = intr.intrinsicInteriorAngle(a, b);
			if (angle < min_angle)
			{
				min_angle = angle;
				return true;
			}
			return false;
		});
		if (it != end(wedge))
			chosen = it;
	}
	return chosen;

	/*	faster finding, the results looks the same but it is preferable to flip the biggest angle firstly
	return std::find_if(begin(wedge), end(wedge), [&](const Dart& d) -> bool {
		if (!intr.can_be_flipped(Edge(d)))
			return false;

		Dart a = phi1(intr.getMesh(), d);
		Dart b = phi<2, -1, 2>(intr.getMesh(), d);
		return (intr.intrinsicInteriorAngle(a, b) < M_PI);
	});
	*/
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

	auto index = indexOfBi(intr, wedge);
	while (index != end(wedge)) // flip until stability
	{
		intr.flip_edge(Edge(*index));
		wedge.erase(index);
		index = indexOfBi(intr, wedge);
	}

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
 * @param iteration optional, the maximum number of flip out to do, default 100
 * @param loop optional, true if the path is closed, default false
*/
void geodesic_path(IntrinsicTriangulation& intr, std::list<Edge>& path,
	 int iteration = 100, bool loop = false)
{
	const CMap2& mesh = intr.getMesh();
	
	auto joints_priority_queue = buildJointList(intr, path, loop);
	while (joints_priority_queue.size() > 0 && iteration > 0)
	{
		Joint j = joints_priority_queue.top();
		Scalar old_length = intr.getLength(j.a) + intr.getLength(j.b);
		auto shorter_path = flip_out(intr, j);

		// reverse the segment if needed to not break the path
		if (j.inverted)
		{
			shorter_path.reverse();
			for(auto e : shorter_path)
				e = Edge(phi2(mesh, e.dart));
		}

		// update
		if (j.loop)
		{
			path.pop_back();
			path.pop_front();
			path.splice(path.end(), shorter_path);
		}
		else
		{
			std::list<Edge>::iterator it = path.erase(j.edge_ref_a, ++j.edge_ref_b);
			path.splice(it, shorter_path);
		}

		// update the priority queue with the new shorter path, in fact rebuild the whole path as every joint flexibility need to be checked
		joints_priority_queue = buildJointList(intr, path, loop);
		--iteration;
	}
}

/**
 * find an arbitrary edge path from vertex a to b
 * TODO templated version
 * @param mesh a CMap2 mesh
 * @param a begin vertex
 * @param b end param
 * @returns a list of adjacent edges
*/
std::list<Edge> find_path(CMap2& mesh, Vertex a, Vertex b)
{
	std::list<Edge> path;
	
	// initialize every distance to inf
	std::shared_ptr<Attribute<Scalar>> attr_dist = add_attribute<Scalar, Vertex>(mesh, "__dist__");
	parallel_foreach_cell(mesh, [&](Vertex v) -> bool {
		value<Scalar>(mesh, attr_dist, v) = std::numeric_limits<Scalar>::max();
		return true;
	});
	value<Scalar>(mesh, attr_dist, b) = 0;

	Vertex v_current;
	std::stack<Vertex> stack;
	Scalar shoortest_path = std::numeric_limits<Scalar>::max();
	stack.emplace(b);

	// forward crossing setting the distance by neighborhood
	do {
		v_current = stack.top();
		stack.pop();
		Scalar current_dist = value<Scalar>(mesh, attr_dist, v_current);
		if (shoortest_path > current_dist)
			foreach_adjacent_vertex_through_edge(mesh, v_current, [&](Vertex v) -> bool {
				if (value<Scalar>(mesh, attr_dist, v) > current_dist + 1) {
					stack.emplace(v);
					value<Scalar>(mesh, attr_dist, v) = current_dist + 1;
				}
				if (isSameVertex(mesh, v.dart, a.dart)) {
					shoortest_path = current_dist+1;
					return false;
				}
				return true;
			});
	} while (!stack.empty());

	v_current = a;
	Scalar min_value = std::numeric_limits<Scalar>::max();

	// backward crossing looking at the minimal path
	do {
		Vertex nearest_vertex = v_current;
		Vertex next_vertex;
		Edge nearest_edge;
		foreach_incident_edge(mesh, v_current, [&](Edge e) -> bool {
			auto v_incidents = incident_vertices(mesh, e);
			next_vertex = v_incidents[0];
			if (isSameVertex(mesh, v_current.dart, next_vertex.dart))
				next_vertex = v_incidents[1];
			Scalar v_dist = value<Scalar>(mesh, attr_dist, next_vertex);
			if (min_value > v_dist) {
				min_value = v_dist;
				nearest_edge = e;
				nearest_vertex = next_vertex;
			}
			return 1;
		});
		path.push_back(nearest_edge);
		v_current = nearest_vertex;
	} while (min_value > 0);

	remove_attribute<Vertex>(mesh, attr_dist);
	return path;
}

/**
 * Reorder a given list of edge so it can be computed by geodesic_path
 * @param mesh a mesh
 * @param updated_path the path that will be updated
 * @returns false if the path could not be reordered when the path is not a connected component
*/
bool reorder_path(CMap2& mesh, std::list<Edge>& updated_path)
{
	// TODO
	return false;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_GEODESIC_H_
