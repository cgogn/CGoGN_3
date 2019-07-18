/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/IterativeLinearSolvers>

namespace cgogn
{

namespace geometry
{

template <typename T, typename MESH>
void filter_average(
	const MESH& m,
	const typename mesh_traits<MESH>::template Attribute<T>* attribute_in,
	typename mesh_traits<MESH>::template Attribute<T>* attribute_out
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	parallel_foreach_cell(m, [&] (Vertex v) -> bool
	{
		T sum;
		sum.setZero();
		uint32 count = 0;
		foreach_adjacent_vertex_through_edge(m, v, [&] (Vertex av) -> bool
		{
			sum += value<T>(m, attribute_in, av);
			++count;
			return true;
		});
		value<T>(m, attribute_out, v) = sum / count;
		return true;
	});
}

//template <typename MAP, typename MASK, typename VERTEX_ATTR>
//void filter_bilateral(
//	const MAP& map,
//	const MASK& mask,
//	const VERTEX_ATTR& position_in,
//	VERTEX_ATTR& position_out,
//	const VERTEX_ATTR& normal
//)
//{
//	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position_in, position_out & normal must be a vertex attribute");

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
//			Scalar wcs = std::exp((-1.0 * (t * t) / (2.0 * sigmaC * sigmaC)) + (-1.0 * (h * h) / (2.0 * sigmaS * sigmaS)));
//			sum += wcs * h;
//			normalizer += wcs;
//		});

//		position_out[v] = position_in[v] + ((sum / normalizer) * n);
//	},
//	mask);
//}

//template <typename MAP, typename MASK, typename VERTEX_ATTR>
//void filter_taubin(
//	const MAP& map,
//	const MASK& mask,
//	VERTEX_ATTR& position,
//	VERTEX_ATTR& position_tmp)
//{
//	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position & position_tmp must be a vertex attribute");

//	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
//	using Scalar = ScalarOf<VEC3>;
//	using Vertex = typename MAP::Vertex;

//	const Scalar lambda = 0.6307;
//	const Scalar mu = -0.6532;

//	map.parallel_foreach_cell([&] (Vertex v)
//	{
//		VEC3 avg;
//		set_zero(avg);
//		uint32 count = 0;
//		map.foreach_adjacent_vertex_through_edge(v, [&] (Vertex av)
//		{
//			avg += position[av];
//			++count;
//		});
//		avg /= Scalar(count);
//		const VEC3& p = position[v];
//		position_tmp[v] = p + ((avg - p) * lambda);
//	},
//	mask);

//	map.parallel_foreach_cell([&] (Vertex v)
//	{
//		VEC3 avg;
//		set_zero(avg);
//		uint32 count = 0;
//		map.foreach_adjacent_vertex_through_edge(v, [&] (Vertex av)
//		{
//			avg += position_tmp[av];
//			++count;
//		});
//		avg /= Scalar(count);
//		const VEC3& p = position_tmp[v];
//		position[v] = p + ((avg - p) * mu);
//	},
//	mask);
//}

//template <typename MAP, typename MASK, typename VERTEX_ATTR>
//void filter_laplacian(
//	MAP& map,
//	const MASK& mask,
//	VERTEX_ATTR& position_in,
//	VERTEX_ATTR& position_out
//)
//{
//	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position_in & position_out must be a vertex attribute");

//	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
//	using Scalar = ScalarOf<VEC3>;
//	using CDart = typename MAP::CDart;
//	using Vertex = typename MAP::Vertex;
//	using Edge = typename MAP::Edge;

//	typename MAP::template EdgeAttribute<Scalar> edge_weight = map.template add_attribute<Scalar, Edge::ORBIT>("__edge_weight");
//	typename MAP::template VertexAttribute<uint32> vertex_index = map.template add_attribute<uint32, Vertex::ORBIT>("__vertex_index");
//	typename MAP::template VertexAttribute<VEC3> vertex_lapl = map.template add_attribute<VEC3, Vertex::ORBIT>("__vertex_lapl");

//	// compute edge weights
//	map.parallel_foreach_cell([&] (Edge e)
//	{
//		if (!map.is_incident_to_boundary(e))
//		{
//			edge_weight[e] = (
//				std::tan(M_PI_2 - angle(map, CDart(map.phi_1(e.dart)), position_in)) +
//				std::tan(M_PI_2 - angle(map, CDart(map.phi_1(map.phi2(e.dart))), position_in))
//			) / 2.0;
//		}
//		else
//		{
//			cgogn::Dart d = map.boundary_dart(e);
//			edge_weight[e] = std::tan(M_PI_2 - angle(map, CDart(map.phi_1(map.phi2(d))), position_in));
//		}
//	},
//	mask);

//	// compute vertices laplacian
//	uint32 nb_vertices = 0;
//	map.foreach_cell([&] (Vertex v) { vertex_index[v] = nb_vertices++; }, mask);
//	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL(nb_vertices, nb_vertices);
//	std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
//	LAPLcoeffs.reserve(nb_vertices * 10);
//	map.foreach_cell([&] (Vertex v)
//	{
//		uint32 vidx = vertex_index[v];
//		Scalar wsum = 0;
//		Scalar a = 1. / (4. * area(map, v, position_in));
//		map.foreach_incident_edge(v, [&] (Edge e) { wsum += a * edge_weight[e]; });
//		map.foreach_adjacent_vertex_through_edge(v, [&] (Vertex av)
//		{
//			LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(vidx, vertex_index[av], (a * edge_weight[Edge(av.dart)]) / wsum));
//		});
//		LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(vidx, vidx, -1.));
//	},
//	mask);
//	LAPL.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());

//	Eigen::MatrixXd vpos(nb_vertices, 3);
//	map.parallel_foreach_cell([&] (Vertex v)
//	{
//		const VEC3& pv = position_in[v];
//		uint32 vidx = vertex_index[v];
//		vpos(vidx, 0) = pv[0];
//		vpos(vidx, 1) = pv[1];
//		vpos(vidx, 2) = pv[2];
//	},
//	mask);

//	Eigen::MatrixXd lapl(nb_vertices, 3);
//	lapl = LAPL * vpos;
//	map.parallel_foreach_cell([&] (Vertex v)
//	{
//		VEC3& vl = vertex_lapl[v];
//		uint32 vidx = vertex_index[v];
//		vl[0] = lapl(vidx, 0);
//		vl[1] = lapl(vidx, 1);
//		vl[2] = lapl(vidx, 2);
//	},
//	mask);

//	// vertices displacement
//	map.parallel_foreach_cell([&] (Vertex v)
//	{
//		position_out[v] = position_in[v] + 0.1 * vertex_lapl[v];
//	},
//	mask);

//	map.remove_attribute(edge_weight);
//	map.remove_attribute(vertex_index);
//	map.remove_attribute(vertex_lapl);
//}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_FILTERING_H_
