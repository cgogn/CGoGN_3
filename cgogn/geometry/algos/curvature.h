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

#ifndef CGOGN_GEOMETRY_ALGOS_CURVATURE_H_
#define CGOGN_GEOMETRY_ALGOS_CURVATURE_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/functions/vector_ops.h>

namespace cgogn
{

namespace geometry
{

template <typename VEC3, typename MESH>
std::tuple<vector_traits<VEC3>::Scalar, vector_traits<VEC3>::Scalar, VEC3, VEC3, VEC3>
curvature(
	const MESH& m,
	typename mesh_traits<MESH>::Vertex v,
	typename vector_traits<VEC3>::Scalar radius,
	const typename mesh_traits<MESH>::template AttributePtr<VEC3> vertex_position,
	const typename mesh_traits<MESH>::template AttributePtr<VEC3> vertex_normal,
	const typename mesh_traits<MESH>::template AttributePtr<typename vector_traits<VEC3>::Scalar>>& edge_angle,
	const typename mesh_traits<MESH>::template AttributePtr<typename vector_traits<VEC3>::Scalar>>& edge_area
)
{
	using Scalar = typename vector_traits<VEC3>::Scalar;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	// collect the normal cycle tensor
	geometry::Collector_WithinSphere<VEC3, MAP> neighborhood(map, radius, position);
	neighborhood.collect(v);

	Eigen::Matrix3d tensor;
	tensor.setZero();

	map.foreach_cell(
		[&] (Edge e)
		{
			std::vector<Vertex> vv = incident_vertices(m, e);
			const VEC3& p1 = value<VEC3>(m, vertex_position, vv[0]);
			const VEC3& p2 = value<VEC3>(m, vertex_position, vv[1]);
			VEC3 ev = p2 - p1;
			tensor += (ev * transpose(ev)) * value<Scalar>(m, edge_angle, e) * (Scalar(1) / norm(ev));
		},
		neighborhood
	);

	neighborhood.foreach_border([&] (Dart d)
	{
		std::pair<Vertex2, Vertex2> vv = map.vertices(Edge2(d));
		const VEC3& p1 = position[vv.first];
		const VEC3& p2 = position[vv.second];
		Eigen::Vector3d ev = Eigen::Vector3d(p2[0], p2[1], p2[2]) - Eigen::Vector3d(p1[0], p1[1], p1[2]);
		Scalar alpha;
		geometry::intersection_sphere_segment<VEC3>(position[v], radius, position[Vertex2(d)], position[Vertex2(map.phi1(d))], alpha);
		tensor += (ev * ev.transpose()) * edge_angle[Edge2(d)] * (Scalar(1) / ev.norm()) * alpha;
	});

	tensor /= neighborhood.area(position);

	const VEC3& normal_v = normal[v];
	Eigen::Vector3d e_normal_v(normal_v[0], normal_v[1], normal_v[2]);

	// project the tensor
	Eigen::Matrix3d proj;
	proj.setIdentity();
	proj -= e_normal_v * e_normal_v.transpose();
	tensor = proj * tensor * proj;

	// solve eigen problem
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(tensor);
	const Eigen::Vector3d& ev = solver.eigenvalues();
	const Eigen::Matrix3d& evec = solver.eigenvectors();

	// sort eigen components : ev[inormal] has minimal absolute value ; kmin = ev[imin] <= ev[imax] = kmax
	uint32 inormal = 0, imin, imax;
	if (fabs(ev[1]) < fabs(ev[inormal])) inormal = 1;
	if (fabs(ev[2]) < fabs(ev[inormal])) inormal = 2;
	imin = (inormal + 1) % 3;
	imax = (inormal + 2) % 3;
	if (ev[imax] < ev[imin]) { std::swap(imin, imax); }

	// set curvatures from sorted eigen components
	// warning : Kmin and Kmax are switched w.r.t. kmin and kmax

	// normal direction : minimal absolute eigen value
	VEC3 Knormal_v { evec(0, inormal), evec(1, inormal), evec(2, inormal) };
	if (Knormal_v.dot(normal_v) < 0)
		Knormal_v *= Scalar(-1); // change orientation

	// min curvature
	Scalar kmin = ev[imin];
	VEC3 Kmin_v { evec(0, imax), evec(1, imax), evec(2, imax) };

	// max curvature
	Scalar kmax = ev[imax];
	VEC3 Kmax_v { evec(0, imin), evec(1, imin), evec(2, imin) };

	return { kmax, kmin, Kmax, Kmin, Knormal };
}

template <typename VEC3, typename MESH>
void
compute_curvature(
	const MESH& m,
	typename vector_traits<VEC3>::Scalar radius,
	const typename mesh_traits<MESH>::template AttributePtr<VEC3> vertex_position,
	const typename mesh_traits<MESH>::template AttributePtr<VEC3> vertex_normal,
	const typename mesh_traits<MESH>::template AttributePtr<typename vector_traits<VEC3>::Scalar>>& edge_angle,
	const typename mesh_traits<MESH>::template AttributePtr<typename vector_traits<VEC3>::Scalar>>& edge_area,
	typename mesh_traits<MESH>::template AttributePtr<typename vector_traits<VEC3>::Scalar>>& vertex_kmax,
	typename mesh_traits<MESH>::template AttributePtr<typename vector_traits<VEC3>::Scalar>>& vertex_kmin,
	typename mesh_traits<MESH>::template AttributePtr<VEC3>& vertex_Kmax,
	typename mesh_traits<MESH>::template AttributePtr<VEC3>& vertex_Kmin,
	typename mesh_traits<MESH>::template AttributePtr<VEC3>& vertex_Knormal
)
{
	using Scalar = typename vector_traits<VEC3>::Scalar;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	foreach_cell(m, [&] (Vertex v) -> bool
	{
		const auto& [kmax, kmin, Kmax, Kmin, Knormal] = curvature<VEC3>(m, v, vertex_position, vertex_normal, edge_angle, edge_area);
		value<Scalar>(m, vertex_kmax, v) = kmax;
		value<Scalar>(m, vertex_kmin, v) = kmin;
		value<VEC3>(m, vertex_Kmax, v) = Kmax;
		value<VEC3>(m, vertex_Kmin, v) = Kmin;
		value<VEC3>(m, vertex_Knormal, v) = Knormal;
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_CURVATURE_H_
