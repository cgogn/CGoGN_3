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

namespace cgogn
{

namespace geometry
{
/**
 * an intrinsic data structure,
 * this class is a simplification of the paper "Navigating Intrinsic Triangulations" (Sharp, Crane, Soliman),
 * it use is only intended for the geodesic computation with the edgeflip function
*/
template <typename MESH,
		  typename std::enable_if_t<std::is_convertible_v<MESH&, cgogn::CMap2&>>* = nullptr>
class IntrinsicTriangulation
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	IntrinsicTriangulation(const CMap2& m, const std::shared_ptr<Attribute<Vec3>> vertex_position): extr_(m), vertex_position_(vertex_position)
	{
		// copy topology (even attributes but unused)
		copy(intr_, m);

		// compute edge length
		edge_length_ = get_or_add_attribute<Scalar, Edge>(intr_, "area");
		parallel_foreach_cell(intr_, [&](Edge c) -> bool {
			std::vector<Vertex> vertices = incident_vertices(m, c);
			const Vec3& a = value<Vec3>(intr_, vertex_position, vertices[0]);
			const Vec3& b = value<Vec3>(intr_, vertex_position, vertices[1]);

			value<Scalar>(intr_, edge_length_, c) = (a - b).squaredNorm();
			return true;
		});
	}

private:
	const MESH& extr_;
	cgogn::CMap2 intr_;
	const std::shared_ptr<Attribute<Vec3>> vertex_position_;	// extrinsic vertex position attribute
	std::shared_ptr<Attribute<Scalar>> edge_length_;	// intrinsic euclidian edge length, L1 norm is discutable
	std::shared_ptr<Attribute<Scalar>> halfedge_angle_;
	std::shared_ptr<Attribute<Dart>> vertex_ref_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_INTRINSIC_TRIANGULATION_H_
