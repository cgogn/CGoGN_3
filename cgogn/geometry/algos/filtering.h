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

#ifndef CGOGN_GEOMETRY_ALGOS_FILTERING_H_
#define CGOGN_GEOMETRY_ALGOS_FILTERING_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Sparse>

namespace cgogn
{

namespace geometry
{

template <typename T, typename MESH>
void filter_average(const MESH& m, const typename mesh_traits<MESH>::template Attribute<T>* attribute_in,
					typename mesh_traits<MESH>::template Attribute<T>* attribute_out)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		T sum;
		sum.setZero();
		uint32 count = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
			sum += value<T>(m, attribute_in, av);
			++count;
			return true;
		});
		value<T>(m, attribute_out, v) = sum / count;
		return true;
	});
}

template <typename MESH>
void filter_regularize(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");
	uint32 nb_vertices = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, vertex_index, v) = nb_vertices++;
		return true;
	});

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL_COT =
		geometry::cotan_laplacian_matrix(m, vertex_index.get(), vertex_position);

	Eigen::MatrixXd vpos(nb_vertices, 3);
	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		const Vec3& pv = value<Vec3>(m, vertex_position, v);
		uint32 vidx = value<uint32>(m, vertex_index, v);
		vpos(vidx, 0) = pv[0];
		vpos(vidx, 1) = pv[1];
		vpos(vidx, 2) = pv[2];
		return true;
	});
	Eigen::MatrixXd poslapl(nb_vertices, 3);
	poslapl = LAPL_COT * vpos;

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * nb_vertices, nb_vertices);
	std::vector<Eigen::Triplet<Scalar>> Acoeffs;
	Acoeffs.reserve(nb_vertices * 10);
	Eigen::MatrixXd b(2 * nb_vertices, 3);

	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		uint32 nbv = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
			uint32 avidx = value<uint32>(m, vertex_index, av);
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
			++nbv;
			return true;
		});
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(nbv)));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), 10.0));
		const Vec3& pv = value<Vec3>(m, vertex_position, v);
		b(vidx, 0) = poslapl(vidx, 0);
		b(vidx, 1) = poslapl(vidx, 1);
		b(vidx, 2) = poslapl(vidx, 2);
		b(nb_vertices + vidx, 0) = 10.0 * pv[0];
		b(nb_vertices + vidx, 1) = 10.0 * pv[1];
		b(nb_vertices + vidx, 2) = 10.0 * pv[2];
		return true;
	});
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At = A.transpose();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(At * A);
	vpos = solver.solve(At * b);

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		Vec3& pos = value<Vec3>(m, vertex_position, v);
		pos[0] = vpos(vidx, 0);
		pos[1] = vpos(vidx, 1);
		pos[2] = vpos(vidx, 2);
		return true;
	});

	remove_attribute<Vertex>(m, vertex_index);
}

// template <typename MAP, typename MASK, typename VERTEX_ATTR>
// void filter_bilateral(
//	const MAP& map,
//	const MASK& mask,
//	const VERTEX_ATTR& position_in,
//	VERTEX_ATTR& position_out,
//	const VERTEX_ATTR& normal
//)
//{
//	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position_in, position_out & normal must be a
// vertex attribute");

//	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
//	using Scalar = ScalarOf<VEC3>;
//	using Vertex = typename MAP::Vertex;
//	using Edge = typename MAP::Edge;

//	Scalar length_sum = 0;
//	Scalar angle_sum = 0;
//	uint32 nb_edges = 0;

//	map.foreach_cell([&] (Edge e)
//	{
//		std::pair<Vertex, Vertex> v = map.vertices(e);
//		VEC3 edge = position_in[v.first] - position_in[v.second];
//		length_sum += edge.norm();
//		angle_sum += angle(normal[v.first], normal[v.second]);
//		++nb_edges;
//	},
//	mask);

//	Scalar sigmaC = 1.0 * (length_sum / Scalar(nb_edges));
//	Scalar sigmaS = 2.5 * (angle_sum / Scalar(nb_edges));

//	map.parallel_foreach_cell([&] (Vertex v)
//	{
//		const VEC3& n = normal[v];

//		Scalar sum = 0, normalizer = 0;
//		map.foreach_adjacent_vertex_through_edge(v, [&] (Vertex av)
//		{
//			VEC3 edge = position_in[av] - position_in[v];
//			Scalar t = edge.norm();
//			Scalar h = n.dot(edge);
//			Scalar wcs = std::exp((-1.0 * (t * t) / (2.0 * sigmaC * sigmaC)) + (-1.0 * (h * h) / (2.0 * sigmaS *
// sigmaS))); 			sum += wcs * h; 			normalizer += wcs;
//		});

//		position_out[v] = position_in[v] + ((sum / normalizer) * n);
//	},
//	mask);
//}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_FILTERING_H_
