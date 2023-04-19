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

#ifndef CGOGN_CORE_CMAP_CMAP_INFO_H_
#define CGOGN_CORE_CMAP_CMAP_INFO_H_

#include <cgogn/core/types/cmap/cmap_base.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

#include <iomanip>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename CMAP>
// uint32 nb_darts_of_orbit(const CMAP& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename CMAP>
uint32 nb_darts_of_orbit(const CMAP& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename CMAP::Cells>::value, "CELL not supported in this CMAP");
	uint32 result = 0;
	foreach_dart_of_orbit(m, c, [&](Dart) -> bool {
		++result;
		return true;
	});
	return result;
}

/*****************************************************************************/

// template <typename CMAP>
// void is_boundary(const CMAP& m, Dart d)

/*****************************************************************************/

//////////////
// CMapBase //
//////////////


} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_INFO_H_
