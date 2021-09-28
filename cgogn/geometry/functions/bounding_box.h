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

#ifndef CGOGN_GEOMETRY_FUNCTIONS_BOUNDING_BOX_H_
#define CGOGN_GEOMETRY_FUNCTIONS_BOUNDING_BOX_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <template <typename VEC> typename CONTAINER, typename VEC>
std::pair<VEC, VEC> bounding_box(const CONTAINER<VEC>& container)
{
	using Scalar = typename vector_traits<VEC>::Scalar;
	const std::size_t dimension = vector_traits<VEC>::SIZE;

	VEC bb_min, bb_max;
	for (std::size_t i = 0; i < dimension; ++i)
	{
		bb_min[i] = std::numeric_limits<Scalar>::max();
		bb_max[i] = std::numeric_limits<Scalar>::lowest();
	}
	for (const VEC& v : container)
	{
		for (std::size_t i = 0; i < dimension; ++i)
		{
			if (v[i] < bb_min[i])
				bb_min[i] = v[i];
			if (v[i] > bb_max[i])
				bb_max[i] = v[i];
		}
	}
	return {bb_min, bb_max};
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_FUNCTIONS_BOUNDING_BOX_H_
