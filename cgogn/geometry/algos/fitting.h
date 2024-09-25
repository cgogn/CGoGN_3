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

#ifndef CGOGN_GEOMETRY_ALGOS_FITTING_H_
#define CGOGN_GEOMETRY_ALGOS_FITTING_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
std::pair<Vec3, Scalar> sphere_fitting(const MESH& m, const std::vector<typename mesh_traits<MESH>::Vertex>& vertices,
									   const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	Eigen::MatrixXd A(vertices.size(), 4);
	Eigen::VectorXd b(vertices.size());
	uint32 idx = 0;
	for (Vertex v : vertices)
	{
		const Vec3& pos = value<Vec3>(m, vertex_position, v);
		A.row(idx) = Eigen::Vector4d(-2.0 * pos[0], -2.0 * pos[1], -2.0 * pos[2], 1.0);
		b(idx) = -(pos[0] * pos[0]) - (pos[1] * pos[1]) - (pos[2] * pos[2]);
		++idx;
	};
	Eigen::LDLT<Eigen::MatrixXd> solver(A.transpose() * A);
	Eigen::MatrixXd s1 = solver.solve(A.transpose() * b);
	s1(3) = std::sqrt(s1(0) * s1(0) + s1(1) * s1(1) + s1(2) * s1(2) - s1(3));

	Vec3 s1c = Vec3(s1(0), s1(1), s1(2));
	Scalar s1r = s1(3);

	Eigen::MatrixXd J(vertices.size(), 4);
	Eigen::VectorXd r(vertices.size());
	Eigen::VectorXd s2(4);
	s2 << s1(0), s1(1), s1(2), s1(3);
	for (uint32 i = 0; i < 5; ++i) // TODO: check number of iterations
	{
		idx = 0;
		for (Vertex v : vertices)
		{
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			Vec3 d = pos - Vec3(s2(0), s2(1), s2(2));
			Scalar l = d.norm();
			J.row(idx) = Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0);
			r(idx) = -(l - s2(3));
			++idx;
		}
		Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
		s2 += solver.solve(J.transpose() * r);
	}

	Vec3 s2c = Vec3(s2(0), s2(1), s2(2));
	Scalar s2r = s2(3);

	return {s2c, s2r};
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_FITTING_H_
