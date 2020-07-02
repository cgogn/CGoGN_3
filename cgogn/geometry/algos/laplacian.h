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

#ifndef CGOGN_GEOMETRY_ALGOS_LAPLACIAN_H_
#define CGOGN_GEOMETRY_ALGOS_LAPLACIAN_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Sparse>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
void compute_laplacian(MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
					   typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_laplacian)
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "MESH dimension should be >= 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	auto vertex_index = get_attribute<uint32, Vertex>(m, "vertex_index");
	if (!vertex_index)
		vertex_index = add_attribute<uint32, Vertex>(m, "vertex_index");

	auto edge_weight = get_attribute<Scalar, Edge>(m, "edge_weight");
	if (!edge_weight)
		edge_weight = add_attribute<Scalar, Edge>(m, "edge_weight");

	// compute edges weight
	parallel_foreach_cell(m, [&](Edge e) -> bool {
		std::vector<Scalar> angles = geometry::opposite_angles(m, e, vertex_position);
		Scalar& weight = value<Scalar>(m, edge_weight, e);
		for (Scalar a : angles)
			weight += std::tan(M_PI_2 - a);
		weight /= uint32(angles.size());
		return true;
	});

	// compute vertices laplacian
	uint32 nb_vertices = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, vertex_index, v) = nb_vertices++;
		return true;
	});

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL(nb_vertices, nb_vertices);
	std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
	LAPLcoeffs.reserve(nb_vertices * 10);
	foreach_cell(m, [&](Edge e) -> bool {
		Scalar w = value<Scalar>(m, edge_weight, e);
		auto vertices = incident_vertices(m, e);
		uint32 vidx1 = value<uint32>(m, vertex_index, vertices[0]);
		uint32 vidx2 = value<uint32>(m, vertex_index, vertices[1]);
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w));
		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w));
		return true;
	});
	LAPL.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());

	Eigen::MatrixXd vpos(nb_vertices, 3);
	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		const Vec3& pv = value<Vec3>(m, vertex_position, v);
		uint32 vidx = value<uint32>(m, vertex_index, v);
		vpos(vidx, 0) = pv[0];
		vpos(vidx, 1) = pv[1];
		vpos(vidx, 2) = pv[2];
		return true;
	});

	Eigen::MatrixXd lapl(nb_vertices, 3);
	lapl = LAPL * vpos;

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		Vec3& dcv = value<Vec3>(m, vertex_laplacian, v);
		uint32 vidx = value<uint32>(m, vertex_index, v);
		dcv[0] = lapl(vidx, 0);
		dcv[1] = lapl(vidx, 1);
		dcv[2] = lapl(vidx, 2);
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_LAPLACIAN_H_
