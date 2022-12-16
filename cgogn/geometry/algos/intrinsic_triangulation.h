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

#ifndef CGOGN_CORE_TYPES_INTRINSIC_TRIANGULATION_H_
#define CGOGN_CORE_TYPES_INTRINSIC_TRIANGULATION_H_

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/angle.h>
#include <cgogn/core/functions/mesh_ops/edge.h>

namespace cgogn
{

namespace geometry
{
/**
 * an intrinsic data structure,
 * this class is a simplification of the paper "Navigating Intrinsic Triangulations" (Sharp, Crane, Soliman),
 * it use is only intended for the geodesic computation with the edgeflip function,
 * only manifold triangulated mesh are supported
*/
class IntrinsicTriangulation
{
	template <typename T>
	using Attribute = typename mesh_traits<CMap2>::template Attribute<T>;

	using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using Edge = typename mesh_traits<CMap2>::Edge;
	using Face = typename mesh_traits<CMap2>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	/*
	template <typename MESH>
	IntrinsicTriangulation(const MESH& m, const std::shared_ptr<Attribute<Vec3>> vertex_position)
	{
		// TODO
		// convert MESH to CMap2
		// IntrinsicTriangulation(m_converted, vertex_position):
	}
	*/
	
	/**
	 * construct an intrinsic triangulation
	 * @param m an extrinsic mesh
	 * @param vertex_position the position attribute related to m
	*/
	IntrinsicTriangulation(const CMap2& m, const std::shared_ptr<Attribute<Vec3>> vertex_position): extr_(m), vertex_position_(vertex_position)
	{
		// copy topology for cmap2
		copy(intr_, m);

		// compute edge length
		edge_length_ = add_attribute<Scalar, Edge>(intr_, "intr_length");
		parallel_foreach_cell(intr_, [&](Edge e) -> bool {
			const Vec3& a = value<Vec3>(intr_, vertex_position, Vertex(e.dart));
			const Vec3& b = value<Vec3>(intr_, vertex_position, Vertex(phi2(intr_, e.dart)));

			value<Scalar>(intr_, edge_length_, e) = (a - b).norm();
			return true;
		});

		// compute angle sum
		vertex_angle_sum_ = add_attribute<Scalar, Vertex>(intr_, "intr_angle_sum");
		vertex_ref_ = add_attribute<Dart, Vertex>(intr_, "intr_ref");
		halfedge_angle_ = add_attribute<Scalar, HalfEdge>(intr_, "intr_angle");
		parallel_foreach_cell(intr_, [&](Vertex v) -> bool {
			// sum interior angles
			value<Scalar>(intr_, vertex_angle_sum_, v) = vertex_angle_sum(extr_, vertex_position_.get(), v);
			// determine extrinsic reference
			value<Dart>(intr_, vertex_ref_, v) = v.dart;

			// compute angle of each halfedge around vertex
			Dart vdart = v.dart;
			Dart it = phi<-1, 2>(intr_, vdart);
			value<Scalar>(intr_, halfedge_angle_, HalfEdge(vdart)) = 0;	// reference direction
			while (it != vdart)
			{
				update_signpost(it);
				it = phi<-1, 2>(intr_, it);
			}
			return true;
		});
			
	}

	/**
	 * @returns the intrinsic mesh
	*/
	CMap2& getMesh()
	{
		return intr_;
	}

	/**
	 * @param d an intrinsic halfedge
	 * @returns the intrinsic angle of halfedge d
	*/
	Scalar getAngle(Dart d)
	{
		return value<Scalar>(intr_, halfedge_angle_, HalfEdge(d));
	}

	/**
	 * counterclockwise intrinsic interior angle from a to b
	 * @param a the first angle
	 * @param b the second angle
	 * @returns angle from a to b
	*/
	Scalar interiorAngle(Dart a, Dart b)
	{
		Scalar a_angle = getAngle(a);
		Scalar b_angle = getAngle(b);
		return (a_angle > b_angle) ? 2 * M_PI + b_angle - a_angle : b_angle - a_angle;
	}

	// from https://geometry-central.net/surface/intrinsic_triangulations/basics/
	inline Vec2 layoutTriangleVertex(Vec2& pA, Vec2& pB, Scalar& lBC, Scalar& lCA)
	{
		Scalar lAB = (pB - pA).norm();
		Scalar tArea = area_(lAB, lBC, lCA);
		Scalar h = 2. * tArea / lAB;
		Scalar w = std::sqrt(std::max(0., lCA * lCA - h * h));
		if (lBC * lBC > (lAB * lAB + lCA * lCA)) w *= -1.0;
		Vec2 vABn = (pB - pA) / lAB;
		Vec2 vABnPerp{-(vABn[1]), vABn[0]};
		Vec2 pC = pA + w * vABn + h * vABnPerp;
		return pC;
	}

	// from https://geometry-central.net/surface/intrinsic_triangulations/basics/
	inline std::array<Vec2, 4> layoutDiamond(Dart iHe)
	{
		Dart iHeA0 = iHe;
		Dart iHeA1 = phi1(intr_, iHeA0);
		Dart iHeA2 = phi1(intr_, iHeA1);
		Dart iHeB0 = phi2(intr_, iHe);
		Dart iHeB1 = phi1(intr_, iHeB0);
		Dart iHeB2 = phi1(intr_, iHeB1);

		double l01 = value<Scalar>(intr_, edge_length_, Edge(iHeA1));
		double l12 = value<Scalar>(intr_, edge_length_, Edge(iHeA2));
		double l23 = value<Scalar>(intr_, edge_length_, Edge(iHeB1));
		double l30 = value<Scalar>(intr_, edge_length_, Edge(iHeB2));
		double l02 = value<Scalar>(intr_, edge_length_, Edge(iHeA0));

		Vec2 p3{0., 0.};
		Vec2 p0{l30, 0.};
		Vec2 p2 = layoutTriangleVertex(p3, p0, l02, l23);
		Vec2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

		return {p0, p1, p2, p3};
	}

	/**
	 * check if an intrinsic edge can be flipped
	 * @param e the edge to flip
	 * @returns true if can be flipped, false otherwise
	*/
	bool can_be_flipped (Edge e)
	{
		if (!edge_can_flip(intr_, e))	// topology can flip
			return false;

		std::array<Vec2, 4> layoutPositions = layoutDiamond(e.dart);

		// Test if geometryically flippable flippable (both signed areas of new triangles are positive)
		Scalar A1 = cross_(layoutPositions[1] - layoutPositions[0], layoutPositions[3] - layoutPositions[0]);
		Scalar A2 = cross_(layoutPositions[3] - layoutPositions[2], layoutPositions[1] - layoutPositions[2]);
		if (A1 <= 0 || A2 <= 0) {
			return false;
		}
		return true;
	}

	/**
	 * flip an intrinsic edge, does not check if it is flippable
	 * @param e the edge to flip
	*/
	void flip_edge (Edge e)
	{
		Scalar future_length = flipped_length_(e);
		cgogn::flip_edge(intr_, e);	// flip topology
		value<Scalar>(intr_, edge_length_, e) = future_length;
		update_signpost(e.dart);
		update_signpost(phi2(intr_, e.dart));
	}

	/**
	 * compute the topologic path of a list of edge
	 * @param e_list a list of adjacent edges to trace
	 * @returns a list of pairs of positions
	*/
	std::vector<Vec3> edge_list_topology(std::list<Edge> e_list)
	{
		std::vector<Vec3> geodesic_segments;
		for (Edge e : e_list)
			for (Vertex v : incident_vertices(intr_, e))
				geodesic_segments.push_back(value<Vec3>(intr_, vertex_position_, v));
		
		return geodesic_segments;
	}

	/**
	 * compute the path of a list of edge to trace
	 * @param e_list a list of adjacent edges to trace
	 * @returns a list of pairs of positions
	*/
	std::vector<Vec3> edge_list_trace(std::list<Edge> e_list)
	{
		std::vector<Vec3> geodesic_segments;
		for (Edge e : e_list)
		{
			std::vector<Vec3> points = trace(e.dart);
			Vec3 last_p = points[0];
			for(Vec3 p : points)
			{
				if (p != last_p)
				{
					geodesic_segments.push_back(last_p);
					geodesic_segments.push_back(p);
					last_p = p;
				}
			}
		}
		return geodesic_segments;
	}

	/**
	 * compute the extrinsic representation of an intrinsic edge
	 * @param d the edge to trace
	 * @returns a list of extrinsic positions
	*/
	std::vector<Vec3> trace(Dart d, int max_intersection = 10)
	{
		std::vector<Vec3> positions;
		positions.push_back(value<Vec3>(intr_, vertex_position_, Vertex(d)));

		// "un-normalize" the angle to get the trace direction
		trace_parameters.crossed_dart = value<Dart>(intr_, vertex_ref_, Vertex(d));
		update_trace_parameter_through_vertex_(value<Scalar>(intr_, halfedge_angle_, HalfEdge(d)));

		// the traced curve length
		Scalar length_remaining = value<Scalar>(intr_, edge_length_, Edge(d));

		// compute first intersection from vertex
		trace_parameters.edge_crossed = false;
		trace_parameters.crossed_position = positions.back();

		int n_intersection = 0;
		// iterate over intersections
		while (n_intersection < max_intersection && length_remaining > epsilon)
		{
			if (trace_parameters.edge_crossed)
				trace_from_edge_();
			else
				trace_from_vertex_();
			length_remaining -= (positions.back() - trace_parameters.crossed_position).norm();
			positions.push_back(trace_parameters.crossed_position);
			++n_intersection;
		}

		if (length_remaining < -epsilon)
		{
			positions.pop_back();
			positions.push_back(value<Vec3>(intr_, vertex_position_, Vertex(phi2(intr_, d))));
		}

		return positions;
	}

private:
	void trace_from_edge_()
	{
		Scalar barycenter, beta, B_gamma, C_gamma;
		Vec3 B_inter, C_inter;
		bool C_side = false, B_side = false;
		
		Scalar alpha = trace_parameters.intersection_angle;
		Dart B_dart = trace_parameters.crossed_dart;
		Dart A_dart = phi_1(extr_, B_dart);
		Dart C_dart = phi1(extr_, B_dart);

		Vec3 X = trace_parameters.crossed_position;	// last_intersection
		Vec3 A = value<Vec3>(extr_, vertex_position_, Vertex(A_dart));
		Vec3 B = value<Vec3>(extr_, vertex_position_, Vertex(B_dart));
		Vec3 C = value<Vec3>(extr_, vertex_position_, Vertex(C_dart));

		beta = geometry::angle((B-C).normalized(), (A-C).normalized());
		C_side = beta + alpha < M_PI;	// if edge AC is intersected
		if(C_side)
		{
			Vec3 AC = A-C;
			C_gamma = M_PI - (alpha + geometry::angle(AC.normalized(), (B-C).normalized()));
			Scalar CX = (X-C).norm() * sin(alpha) / sin(C_gamma);	// distance between C_inter and C
			barycenter = CX / (AC).norm();
			if (barycenter > 1-epsilon) {
				C_side = false;
			} else {
				C_inter = C*(1-barycenter) + A*barycenter;
			}
		}

		beta = geometry::angle((C-B).normalized(), (A-B).normalized());
		alpha = M_PI - alpha;
		B_side = beta + alpha < M_PI ;	// if edge AB is intersected
		if(B_side)
		{
			Vec3 AB = A-B;
			B_gamma = M_PI - (alpha + geometry::angle(AB.normalized(), (C-B).normalized()));
			Scalar BX = (X-B).norm() * sin(alpha) / sin(B_gamma);	// distance between B_inter and C
			barycenter = BX / (AB).norm();
			if (barycenter > 1-epsilon){
				B_side = false;
			} else {
				B_inter = B*(1-barycenter) + A*barycenter;
			}
		}

		if (C_side && B_side)	// pick the nearest intersection
		{
			if ((X-C_inter).norm() < (X-B_inter).norm())
				B_side = false;
			else
				C_side = false;
		} 
		if (C_side)
		{
			trace_parameters.intersection_angle = M_PI - C_gamma;
			trace_parameters.crossed_position = C_inter;
			trace_parameters.crossed_dart = phi2(extr_, C_dart);
		}
		else if (B_side)
		{
			trace_parameters.intersection_angle = B_gamma;
			trace_parameters.crossed_position = B_inter;
			trace_parameters.crossed_dart = phi2(extr_, A_dart);
		} else {
			trace_parameters.edge_crossed = false;
			trace_parameters.crossed_position = A;
			trace_parameters.crossed_dart = phi<1,2>(extr_, B_dart);
			update_trace_parameter_through_vertex_(M_PI, geometry::angle((C-A).normalized(), (X-A).normalized()));
		}
	}

	void trace_from_vertex_()
	{
		Scalar angle = trace_parameters.intersection_angle;
		Dart A_dart = trace_parameters.crossed_dart;
		Dart C_dart = phi1(extr_, A_dart);
		Dart B_dart = phi2(extr_, C_dart);
		Vec3 A = trace_parameters.crossed_position;
		Vec3 B = value<Vec3>(extr_, vertex_position_, Vertex(B_dart));
		Vec3 C = value<Vec3>(extr_, vertex_position_, Vertex(C_dart));
		Vec3 CB = C-B;
		Vec3 AB = A-B;
		Scalar lambda = M_PI - (angle + geometry::angle(AB.normalized(), CB.normalized()));
		Scalar BX = AB.norm() * sin(angle) / sin(lambda);
		Scalar barycenter = BX / CB.norm();
		
		if (barycenter < epsilon) {	// intersect B
			trace_parameters.crossed_position = B;
			trace_parameters.crossed_dart = phi_1(extr_, A_dart);
			update_trace_parameter_through_vertex_(M_PI);
		}
		else if (barycenter > 1-epsilon) {	// intersect C
			trace_parameters.crossed_position = C;
			trace_parameters.crossed_dart = phi2(extr_, A_dart);
			update_trace_parameter_through_vertex_(M_PI);
		}
		else {	// intersect BC edge
			trace_parameters.crossed_position = B*(1-barycenter) + C*barycenter;
			trace_parameters.edge_crossed = true;
			trace_parameters.crossed_dart = B_dart;
			trace_parameters.intersection_angle = lambda;
		}
	}

	/**
	 * called when a vertex is intersected
	 * update trace_parameters next angle and dart when intersecting a vertex
	 * it will update the angle and dart such as trace_parameters.crossed_dart will be turned at angle_max
	 * @param angle_max normalized angle for rotation
	 * @param offset unnormalize angle, used when the initial direction is not a vertex
	*/
	void update_trace_parameter_through_vertex_(Scalar angle_max, Scalar offset = 0)
	{
		Dart it = trace_parameters.crossed_dart;
		Vec3 v_position = value<Vec3>(extr_, vertex_position_, Vertex(it));
		Vec3 current_dir = value<Vec3>(extr_, vertex_position_, Vertex(phi2(extr_,it))) - v_position;
		current_dir.normalize();
		Scalar angle_sum = vertex_angle_sum(extr_, vertex_position_.get(), Vertex(it));
		Scalar normalized_coeff = 2 * M_PI / angle_sum;
		Vec3 next_dir;
		Scalar cumulated_angle {offset};
		while (cumulated_angle < angle_max) {
			it = phi_1(extr_, it);
			next_dir = value<Vec3>(extr_, vertex_position_, Vertex(it)) - v_position;
			it = phi2(extr_, it);
			next_dir.normalize();
			cumulated_angle += normalized_coeff * geometry::angle(current_dir, next_dir);
			current_dir = next_dir;
		}
		trace_parameters.crossed_dart = phi<2,1>(extr_, it);
		trace_parameters.intersection_angle = angle_sum * (cumulated_angle - angle_max) / (M_PI * 2); // un-normalize
	}

	/**
	 * update the angle of dart ik with the 3 edges lengths of the triangle ijk
	 * @param jk the HalfEdge angle to update
	*/
	void update_signpost(Dart jk) {
		Dart ki = phi<2, 1>(intr_, jk);
		Scalar angle = value<Scalar>(intr_, halfedge_angle_, HalfEdge(ki));
		// normalize angle
		angle += 2 * M_PI * angle_(Edge(jk), Edge(ki), Edge(phi1(intr_, ki))) / value<Scalar>(intr_, vertex_angle_sum_, Vertex(jk));
		value<Scalar>(intr_, halfedge_angle_, HalfEdge(jk)) = std::fmod(angle, 2 * M_PI);
	}

	/**
	 * compute the interior angle of a triangle on point k from its lengths by Kahan's formula
	 * @param ij edge adjacent to jk and ki, opposite to k
	 * @param ki edge adjacent to ij and jk, incident to k
	 * @param jk edge adjacent to ij and ki, incident to k
	 * @returns the interior angle of ijk on point k incident to jk and ki
	*/
	inline Scalar angle_(Edge jk, Edge ki, Edge ij) {
		Scalar a = value<Scalar>(intr_, edge_length_, jk);
		Scalar b = value<Scalar>(intr_, edge_length_, ki);
		Scalar c = value<Scalar>(intr_, edge_length_, ij);
		Scalar mu;
		if (c >= 0 && b >= c)
			mu = c - a + b;
		else if (b >= 0 && c > b)
			mu = b - a + c;
		else
			assert(false);
		
		return 2 * atan(sqrt( ((a-b+c)*mu) / ((a+b+c)*(a-c+b)) ));
	}
	
	/**
	 * compute the area of a triangle face by Heron's formula
	 * @param f a triangle face
	 * @returns the area of the triangle
	*/
	inline Scalar area_(Face f) {
		Scalar a = value<Scalar>(intr_, edge_length_, Edge(f.dart));
		Scalar b = value<Scalar>(intr_, edge_length_, Edge(phi1(intr_, f.dart)));
		Scalar c = value<Scalar>(intr_, edge_length_, Edge(phi_1(intr_, f.dart)));
		Scalar s = (a+b+c) * Scalar(0.5);
		return sqrt(s * (s-a) * (s-b) * (s-c));
	}

	/**
	 * compute the area of a triangle face by Heron's formula
	 * @param a a triangle length
	 * @param b a triangle length
	 * @param c a triangle length
	 * @returns the area of the triangle
	*/
	inline Scalar area_(Scalar a, Scalar b, Scalar c) {
		Scalar s = (a+b+c) * Scalar(0.5);
		return sqrt(s * (s-a) * (s-b) * (s-c));
	}

	/**
	 * compute the future flipped length of an edge
	 * @param e an edge befor flipping
	 * @returns the length of this edge after flipping
	*/
	Scalar flipped_length_(Edge e) {
		std::vector<Face> faces = incident_faces(intr_, e);
		Scalar ij = value<Scalar>(intr_, edge_length_, e);
		ij*=ij;
		Dart d = e.dart;
		Dart opp = phi2(intr_, e.dart);
		Scalar im = value<Scalar>(intr_, edge_length_, Edge(phi1(intr_, d)));
		im*=im;
		Scalar jm = value<Scalar>(intr_, edge_length_, Edge(phi_1(intr_, d)));
		jm*=jm;
		Scalar jk = value<Scalar>(intr_, edge_length_, Edge(phi1(intr_, opp)));
		jk*=jk;
		Scalar ki = value<Scalar>(intr_, edge_length_, Edge(phi_1(intr_, opp)));
		ki*=ki;

		Scalar km = (im+jk+jm+ki-ij + (im-jm) * (jk-ki) / ij) * Scalar(0.5);
		km += (Scalar(8) * area_(faces[0]) * area_(faces[1])) / ij;
		return sqrt(km);
	}

	/**
	 * Compute the dual of the cross product
	*/
	inline Scalar cross_(Vec2 u, Vec2 v)
	{
		return u[0] * v[1] - u[1] * v[0];
	}

private:
	static constexpr Scalar epsilon = 0.00001;
	
	struct {
		Vec3 crossed_position;
		Dart crossed_dart;
		Scalar intersection_angle;
		bool edge_crossed = false;
	} trace_parameters;

	const cgogn::CMap2& extr_;									// extrinsic mesh
	const std::shared_ptr<Attribute<Vec3>> vertex_position_;	// extrinsic vertex position attribute

	cgogn::CMap2 intr_;											// intrinsic mesh
	std::shared_ptr<Attribute<Dart>> vertex_ref_;				// extrinsic link of intrinsic vertices
	std::shared_ptr<Attribute<Scalar>> edge_length_;			// intrinsic euclidian edge length
	std::shared_ptr<Attribute<Scalar>> halfedge_angle_;			// interior angle
	std::shared_ptr<Attribute<Scalar>> vertex_angle_sum_;		// vertex angle sum, avoid recomputing
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_INTRINSIC_TRIANGULATION_H_
