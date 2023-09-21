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

#ifndef CGOGN_CORE_TYPES_MAPS_CMAP_CMAP0_H_
#define CGOGN_CORE_TYPES_MAPS_CMAP_CMAP0_H_

#include <cgogn/core/types/maps/cmap/cmap_base.h>

namespace cgogn
{

struct CMap0 : public CMapBase
{
	static const uint8 dimension = 0;

	using Vertex = Cell<DART>;

	using Cells = std::tuple<Vertex>;

	CMap0()
	{
	}
};

template <>
struct mesh_traits<CMap0>
{
	static constexpr const char* name = "CMap0";
	static constexpr const uint8 dimension = 0;

	using Vertex = typename CMap0::Vertex;

	using Cells = std::tuple<Vertex>;
	static constexpr const char* cell_names[] = {"Vertex"};

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

/*************************************************************************/
// Operators
/*************************************************************************/

CMap0::Vertex add_vertex(CMap0& m, bool set_indices = true);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MAPS_CMAP_CMAP0_H_
