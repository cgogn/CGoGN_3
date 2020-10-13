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

#ifndef CGOGN_MODELING_ALGOS_HOLE_FILLING_H_
#define CGOGN_MODELING_ALGOS_HOLE_FILLING_H_

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/functions/orientation.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <array>
#include <vector>

namespace cgogn
{

namespace modeling
{

void fill_holes(CMap2& m)
{
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		if (is_boundary(m, d))
			set_boundary(m, d, false);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_HOLE_FILLING_H_
