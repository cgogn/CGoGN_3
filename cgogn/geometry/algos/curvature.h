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
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/algos/selection.h>
#include <cgogn/geometry/algos/area.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
std::tuple<Scalar, Scalar, Vec3, Vec3, Vec3>
curvature(
	const MESH& m,
	typename mesh_traits<MESH>::Vertex v,
	Scalar radius,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	const typename mesh_traits<MESH>::template Attribute<Scalar>* edge_angle
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	CellCache<MESH> neighborhood = within_sphere(m, v, radius, vertex_position);

	Mat3 tensor;
	tensor.setZero();

	foreach_cell(neighborhood, [&] (Edge e) -> bool
	{
		std::vector<Vertex> vv = incident_vertices(m, e);
		Vec3 ev = value<Vec3>(m, vertex_position, vv[0]) - value<Vec3>(m, vertex_position, vv[1]);
		tensor += (ev * ev.transpose()) * value<Scalar>(m, edge_angle, e) * (Scalar(1) / ev.norm());
		return true;
	});

	tensor /= area(neighborhood, vertex_position);

	const Vec3& vnormal = value<Vec3>(m, vertex_normal, v);

	// project the tensor
	Mat3 proj;
	proj.setIdentity();
	proj -= vnormal * vnormal.transpose();
	tensor = proj * tensor * proj;

	// solve eigen problem
	Eigen::SelfAdjointEigenSolver<Mat3> solver(tensor);
	const Vec3& ev = solver.eigenvalues();
	const Mat3& evec = solver.eigenvectors();

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
	Vec3 Knormal;
	Knormal[0] = evec(0, inormal);
	Knormal[1] = evec(1, inormal);
	Knormal[2] = evec(2, inormal);
	if (Knormal.dot(vnormal) < 0)
		Knormal *= Scalar(-1); // change orientation

	// min curvature
	Scalar kmin = ev[imin];
	Vec3 Kmin;
	Kmin[0] = evec(0, imax);
	Kmin[1] = evec(1, imax);
	Kmin[2] = evec(2, imax);

	// max curvature
	Scalar kmax = ev[imax];
	Vec3 Kmax;
	Kmax[0] = evec(0, imin);
	Kmax[1] = evec(1, imin);
	Kmax[2] = evec(2, imin);

	return { kmax, kmin, Kmax, Kmin, Knormal };
}

template <typename MESH>
void
compute_curvature(
	const MESH& m,
	Scalar radius,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	const typename mesh_traits<MESH>::template Attribute<Scalar>* edge_angle,
	typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_kmax,
	typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_kmin,
	typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_Kmax,
	typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_Kmin,
	typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_Knormal
)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	parallel_foreach_cell(m, [&] (Vertex v) -> bool
	{
		const auto& [kmax, kmin, Kmax, Kmin, Knormal] = curvature(m, v, radius, vertex_position, vertex_normal, edge_angle);
		value<Scalar>(m, vertex_kmax, v) = kmax;
		value<Scalar>(m, vertex_kmin, v) = kmin;
		value<Vec3>(m, vertex_Kmax, v) = Kmax;
		value<Vec3>(m, vertex_Kmin, v) = Kmin;
		value<Vec3>(m, vertex_Knormal, v) = Knormal;
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_CURVATURE_H_
