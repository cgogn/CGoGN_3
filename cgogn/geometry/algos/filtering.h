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

#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Sparse>

namespace cgogn
{

namespace geometry
{

template <typename T, typename MESH>
void filter_average(const MESH& m, const typename mesh_traits<MESH>::template Attribute<T>* vertex_attribute_in,
					typename mesh_traits<MESH>::template Attribute<T>* vertex_attribute_out)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		T sum;
		sum.setZero();
		uint32 count = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
			sum += value<T>(m, vertex_attribute_in, av);
			++count;
			return true;
		});
		value<T>(m, vertex_attribute_out, v) = sum / count;
		return true;
	});
}

template <typename MESH>
void filter_regularize(MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
					   typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_attribute, Scalar fit_to_data)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");
	uint32 nb_vertices = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, vertex_index, v) = nb_vertices++;
		return true;
	});

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL_COT =
		geometry::cotan_operator_matrix(m, vertex_index.get(), vertex_position);

	Eigen::MatrixXd vattr(nb_vertices, 3);
	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		const Vec3& val = value<Vec3>(m, vertex_attribute, v);
		uint32 vidx = value<uint32>(m, vertex_index, v);
		vattr(vidx, 0) = val[0];
		vattr(vidx, 1) = val[1];
		vattr(vidx, 2) = val[2];
		return true;
	});
	Eigen::MatrixXd attr_lapl(nb_vertices, 3);
	attr_lapl = LAPL_COT * vattr;

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * nb_vertices, nb_vertices);
	std::vector<Eigen::Triplet<Scalar>> Acoeffs;
	Acoeffs.reserve(nb_vertices * 10);
	Eigen::MatrixXd b(2 * nb_vertices, 3);

	foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		if (!is_incident_to_boundary(m, v))
		{
			uint32 nbv = 0;
			foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
				uint32 avidx = value<uint32>(m, vertex_index, av);
				Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1.0));
				++nbv;
				return true;
			});
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), (-1.0 * Scalar(nbv))));
			b(vidx, 0) = attr_lapl(vidx, 0);
			b(vidx, 1) = attr_lapl(vidx, 1);
			b(vidx, 2) = attr_lapl(vidx, 2);
		}
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), fit_to_data));
		const Vec3& val = value<Vec3>(m, vertex_attribute, v);
		b(nb_vertices + vidx, 0) = fit_to_data * val[0];
		b(nb_vertices + vidx, 1) = fit_to_data * val[1];
		b(nb_vertices + vidx, 2) = fit_to_data * val[2];
		return true;
	});
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At = A.transpose();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(At * A);
	vattr = solver.solve(At * b);

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(m, vertex_index, v);
		Vec3& pos = value<Vec3>(m, vertex_attribute, v);
		pos[0] = vattr(vidx, 0);
		pos[1] = vattr(vidx, 1);
		pos[2] = vattr(vidx, 2);
		return true;
	});

	remove_attribute<Vertex>(m, vertex_index);
}

template <typename MESH>
void filter_bilateral(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position_in,
					  typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position_out)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	Scalar length_sum = 0;
	Scalar angle_sum = 0;
	uint32 nb_edges = 0;

	foreach_cell(m, [&](Edge e) -> bool {
		length_sum += length(m, e, vertex_position_in);
		angle_sum += angle(m, e, vertex_position_in);
		++nb_edges;
		return true;
	});

	Scalar sigmaC = 1.0 * (length_sum / Scalar(nb_edges));
	Scalar sigmaS = 2.5 * (angle_sum / Scalar(nb_edges));

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		Vec3 n = normal(m, v, vertex_position_in);
		Scalar sum = 0, normalizer = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
			Vec3 edge = value<Vec3>(m, vertex_position_in, av) - value<Vec3>(m, vertex_position_in, v);
			Scalar t = edge.norm();
			Scalar h = n.dot(edge);
			Scalar wcs =
				std::exp((-1.0 * (t * t) / (2.0 * sigmaC * sigmaC)) + (-1.0 * (h * h) / (2.0 * sigmaS * sigmaS)));
			sum += wcs * h;
			normalizer += wcs;
			return true;
		});
		value<Vec3>(m, vertex_position_out, v) = value<Vec3>(m, vertex_position_in, v) + ((sum / normalizer) * n);
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_FILTERING_H_
