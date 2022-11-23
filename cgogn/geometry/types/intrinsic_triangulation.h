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

namespace cgogn
{

namespace geometry
{
/**
 * an intrinsic data structure,
 * this class is a simplification of the paper "Navigating Intrinsic Triangulations" (Sharp, Crane, Soliman),
 * it use is only intended for the geodesic computation with the edgeflip function,
 * only manifold mesh are supported
*/
class IntrinsicTriangulation
{
	template <typename T>
	using Attribute = typename mesh_traits<CMap2>::template Attribute<T>;

	using HalfEdge = typename mesh_traits<CMap2>::HalfEdge;
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using Edge = typename mesh_traits<CMap2>::Edge;

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
		parallel_foreach_cell(intr_, [&](Vertex v) -> bool {
			Scalar angle_sum{0};
			Vec3 vertex = value<Vec3>(intr_, vertex_position, v);
			std::vector<Vertex> vertices = adjacent_vertices_through_edge(m, v);
			
			// sum incident directions angles
			for (uint32 i = 0, size = uint32(vertices.size()); i < size; ++i)
			{
				Vec3 current_vertex = value<Vec3>(m, vertex_position, vertices[(i+1)%size]);
				Vec3 next_vertex = value<Vec3>(m, vertex_position, vertices[i]);
				angle_sum += angle(current_vertex - vertex, next_vertex - vertex);
			}
			value<Scalar>(intr_, vertex_angle_sum_, v) = angle_sum;
			return true;
		});

		// determine direction reference
		vertex_ref_dir_ = add_attribute<Vec3, Vertex>(intr_, "intr_ref_dir");
		parallel_foreach_cell(intr_, [&](Vertex v) -> bool {
			std::vector<Vertex> vertices = incident_vertices(m, v);
			const Vec3& d = value<Vec3>(intr_, vertex_position, vertices[0]);
			const Vec3& o = value<Vec3>(intr_, vertex_position, v);

			value<Vec3>(intr_, vertex_ref_dir_, v) = (d - o);
			return true;
		});

		// compute angle
		halfedge_angle_ = add_attribute<Scalar, HalfEdge>(intr_, "intr_angle");
		parallel_foreach_cell(intr_, [&](HalfEdge h) -> bool {
			value<Scalar>(intr_, halfedge_angle_, h) = angle_(h);
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
		// TODO
	}

private:
	inline Scalar angle_(HalfEdge h) {
		Dart& dart = h.dart;
		const Vec3& a = value<Vec3>(intr_, vertex_position_, Vertex(dart));
		const Vec3& b = value<Vec3>(intr_, vertex_position_, Vertex(phi1(intr_, dart)));
		const Vec3& c = value<Vec3>(intr_, vertex_ref_dir_, Vertex(dart));
		return normalized_angle(Vertex(dart), angle(b-a, c));
	}

	inline Scalar normalized_angle(Vertex v, Scalar angle) {
		return 2 * M_PI * angle / value<Scalar>(intr_, vertex_angle_sum_, v);
	}

private:
	cgogn::CMap2 intr_;											// intrinsic mesh
	const std::shared_ptr<Attribute<Vec3>> vertex_position_;	// extrinsic vertex position attribute
	std::shared_ptr<Attribute<Scalar>> edge_length_;			// intrinsic euclidian edge length, L1 norm is discutable
	std::shared_ptr<Attribute<Vec3>> vertex_ref_dir_;			// direction relative to the extrinsic mesh
	std::shared_ptr<Attribute<Scalar>> halfedge_angle_;			// angle relative to vertex_dir_
	std::shared_ptr<Attribute<Scalar>> vertex_angle_sum_;		// vertex angle sum, avoid recomputing
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_INTRINSIC_TRIANGULATION_H_
