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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_HALFEDGE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_HALFEDGE_H_

// #include <cgogn/core/cgogn_core_export.h>

// #include <cgogn/core/utils/tuples.h>
// #include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_halfedge(MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////


/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::HalfEdge> incident_halfedges(MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::HalfEdge> incident_halfedges(const MESH& m, CELL c)
{
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
	std::vector<HalfEdge> halfedges;
	halfedges.reserve(32u);
	foreach_incident_halfedge(m, c, [&](HalfEdge e) -> bool {
		halfedges.push_back(e);
		return true;
	});
	return halfedges;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_HALFEDGE_H_
