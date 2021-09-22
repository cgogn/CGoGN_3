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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_OPS_GLOBAL_H_
#define CGOGN_CORE_FUNCTIONS_MESH_OPS_GLOBAL_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// void
// clear(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

void clear(CMapBase& m)
{
	m.darts_.clear_attributes();
	for (CMapBase::AttributeContainer& container : m.attribute_containers_)
		container.clear_attributes();
}

/*****************************************************************************/

// template <typename MESH>
// void
// copy(MESH& dst, const MESH& src);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>* = nullptr>
void copy(MESH& dst, const MESH& src)
{
	// TODO
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_GLOBAL_H_
