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
#include <cgogn/geometry/functions/angle.h>
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

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	template <typename MESH>
	IntrinsicTriangulation(const MESH& m, const std::shared_ptr<Attribute<Vec3>> vertex_position)
	{
		// TODO
		// convert MESH to CMap2
		// IntrinsicTriangulation(m_converted, vertex_position):
	}
	
	IntrinsicTriangulation(const CMap2& m, const std::shared_ptr<Attribute<Vec3>> vertex_position): vertex_position_(vertex_position)
	{
		// copy topology for cmap2
		copy(intr_, m);

		// compute edge length
		edge_length_ = add_attribute<Scalar, Edge>(intr_, "intr_length");
		parallel_foreach_cell(intr_, [&](Edge c) -> bool {
			std::vector<Vertex> vertices = incident_vertices(m, c);
			const Vec3& a = value<Vec3>(intr_, vertex_position, vertices[0]);
			const Vec3& b = value<Vec3>(intr_, vertex_position, vertices[1]);

			value<Scalar>(intr_, edge_length_, c) = (a - b).squaredNorm();
			return true;
		});

		// compute angle sum
		vertex_angle_sum_ = add_attribute<Scalar, Vertex>(intr_, "intr_length");
		vertex_ref_dart_ = add_attribute<Dart, Vertex>(intr_, "intr_ref_dir");
		halfedge_angle_ = add_attribute<Scalar, HalfEdge>(intr_, "intr_angle");
		parallel_foreach_cell(intr_, [&](Vertex v) -> bool {
			
			// sum extrinsic incident directions angles
			Scalar angle_sum{0};
			Vec3 v_position = value<Vec3>(intr_, vertex_position, v);
			std::vector<Vertex> vertices = adjacent_vertices_through_edge(m, v);
			for (uint32 i = 0, size = uint32(vertices.size()); i < size; ++i)
			{
				Vec3 current_vertex = value<Vec3>(m, vertex_position, vertices[(i+1)%size]);
				Vec3 next_vertex = value<Vec3>(m, vertex_position, vertices[i]);
				angle_sum += angle(current_vertex - v_position, next_vertex - v_position);
			}
			value<Scalar>(intr_, vertex_angle_sum_, v) = angle_sum;

			// determine extrinsic direction reference
			const Vec3& ref_pos = value<Vec3>(intr_, vertex_position_, Vertex(phi1(intr_, v.dart)));
			value<Vec3>(intr_, vertex_ref_dir_, v) = ref_pos - v_position;

			// compute angle of each halfedge around vertex
			foreach_dart_of_PHI21(intr_, v.dart, [&](Dart d) -> bool {
				HalfEdge h = HalfEdge(d);
				if (d == Vertex(d).dart)
					value<Scalar>(intr_, halfedge_angle_, h) = 0;	// init to 0
				else
					update_signpost(h);
				return true;
			});

			return true;
		});
	}

	CMap2* getMesh()
	{
		return &intr_;
	}

	void flip_out()
	{
		// TODO
	}

	void flip_edge (Edge e)
	{
		assert(edge_can_flip(intr_, e));
		Scalar future_length = flipped_length(e);
		cgogn::flip_edge(intr_, e);	// flip topology
		value<Scalar>(intr_, edge_length_, e) = future_length;
		update_signpost(e.dart);
		update_signpost(phi2(intr_, e.dart));
	}

private:
	/**
	 * update the angle of dart ik with the 3 edges of the triangle ijk
	 * @param ik the HalfEdge angle to update
	*/
	void update_signpost(HalfEdge h) {
		Dart ik = h.dart;
		Dart ij = phi<1, 2>(intr_, ik);
		Scalar angle = value<Scalar>(intr_, halfedge_angle_, HalfEdge(ij));
		angle += normalize_angle_(Vertex(ik), angle_(Edge(ik), Edge(ij), Edge(phi1(intr_, ij))));
		value<Scalar>(intr_, halfedge_angle_, h) = angle;
	}

	/**
	 * compute the interior angle on point i in the triangle ijk
	 * @param ij, jk, ik edges of the triangle
	 * @returns the interior angle of ijk on point i
	*/
	inline Scalar angle_(Edge ij, Edge ik, Edge jk) {
		Scalar a = value<Scalar>(intr_, edge_length_, ij);
		Scalar b = value<Scalar>(intr_, edge_length_, ik);
		Scalar c = value<Scalar>(intr_, edge_length_, jk);
		return acos((a*a + c*c - b*b) / (2*a*c));
	}
	
	/**
	 * compute the area of a triangle face by Heron's formula
	 * @param ik a triangle face
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
	 * compute the future flipped length of an edge
	 * @param e an edge befor flipping
	 * @returns the length of this edge after flipping
	*/
	Scalar flipped_length(Edge e) {
		std::vector<Face> faces = incident_faces(intr_, e);
		assert(faces.size == 2);
		Scalar ij = value<Scalar>(intr_, edge_length_, e);
		ij*=ij;
		// TODO (see intr_tri_course)
		Dart d = e.dart;
		Dart opp = phi2(intr_, e.dart);
		Scalar im = value<Scalar>(intr_, edge_length_, phi1(d));
		im*=im;
		Scalar jm = value<Scalar>(intr_, edge_length_, phi_1(d));
		jm*=jm;
		Scalar jk = value<Scalar>(intr_, edge_length_, phi1(opp));
		jk*=jk;
		Scalar ki = value<Scalar>(intr_, edge_length_, phi_1(opp));
		ki*=ki;

		Scalar km = (im+jk+jm+ki-ij + (im-jm) * (jk-ki) / ij) * Scalar(0.5);
		km += Scalar(8) * area_(faces[0]) * area_(faces[1]) / ij;
		return sqrt(km);
	}

	/**
	 * normalize an angle so the sum of all incident angles is 2*pi
	 * @param v the incident vertex
	 * @param angle the angle value
	 * @returns the normalized angle
	*/
	inline Scalar normalize_angle_(Vertex v, Scalar angle) {
		return 2 * M_PI * angle / value<Scalar>(intr_, vertex_angle_sum_, v);	// could be optimized with 2*M_PI/sum direct value
	}

private:
	cgogn::CMap2 intr_;											// intrinsic mesh
	const std::shared_ptr<Attribute<Vec3>> vertex_position_;	// extrinsic vertex position attribute
	std::shared_ptr<Attribute<Scalar>> edge_length_;			// intrinsic euclidian edge length, L1 norm is discutable
	std::shared_ptr<Attribute<Vec3>> vertex_ref_dir_;			// reference direction for angle
	std::shared_ptr<Attribute<Scalar>> halfedge_angle_;			// interior angle
	std::shared_ptr<Attribute<Scalar>> vertex_angle_sum_;		// vertex angle sum, avoid recomputing
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_INTRINSIC_TRIANGULATION_H_
