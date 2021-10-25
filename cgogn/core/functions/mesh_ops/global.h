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

#include <cgogn/core/types/cmap/cmap_base.h>

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

void clear(CMapBase& m, bool keep_attributes = true)
{
	// clear darts and keep attributes (phi relations)
	m.darts_.clear_attributes();
	if (!keep_attributes)
	{
		// remove cells indices attributes
		for (uint32 orbit = 0; orbit < NB_ORBITS; ++orbit)
		{
			if (m.cells_indices_[orbit] != nullptr)
			{
				m.darts_.remove_attribute(m.cells_indices_[orbit]);
				m.cells_indices_[orbit].reset();
			}
		}
	}

	// clear all cell attributes
	for (CMapBase::AttributeContainer& container : m.attribute_containers_)
	{
		if (keep_attributes)
			container.clear_attributes();
		else
		{
			container.clear_attributes();
			// if there are still shared_ptr somewhere, some attributes may not be removed
			container.remove_attributes();
		}
	}
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
	clear(dst, false);
	for (uint32 orbit = 0; orbit < NB_ORBITS; ++orbit)
	{
		if (src.cells_indices_[orbit] != nullptr)
			init_cells_indexing(dst, Orbit(orbit));
	}
	dst.darts_.copy(src.darts_);
	for (uint32 i = 0; i < NB_ORBITS; ++i)
		dst.attribute_containers_[i].copy(src.attribute_containers_[i]);
	dst.boundary_marker_ = dst.darts_.get_mark_attribute();
	dst.boundary_marker_->copy(*src.boundary_marker_);
}

////////////////////
// IncidenceGraph //
////////////////////

void copy(IncidenceGraph& /*dst*/, const IncidenceGraph& /*src*/)
{
	// TODO
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_OPS_GLOBAL_H_
