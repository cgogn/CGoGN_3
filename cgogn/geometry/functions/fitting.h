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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_FITTING_H_
#define CGOGN_GEOMETRY_FUNCTIONS_FITTING_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

namespace cgogn
{

namespace geometry
{

inline std::pair<Vec3, Scalar> plane_fitting(const std::vector<Vec3>& points)
{
	uint32 nb_points = points.size();
	Eigen::Matrix3Xd m_points(3, nb_points);
	for (uint32 i = 0; i < nb_points; ++i)
		m_points.col(i) = points[i];

	Vec3 mean = m_points.rowwise().mean();
	const Eigen::Matrix3Xd m_points_centered = m_points.colwise() - mean;

	Eigen::JacobiSVD<Eigen::Matrix3Xd> svd = m_points_centered.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

	Vec3 normal = svd.matrixU().col(2);
	Scalar d = normal.dot(mean);

	normal.normalize();

	return {normal, -d};
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_FITTING_H_
