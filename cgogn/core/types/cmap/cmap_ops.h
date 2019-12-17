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

#ifndef CGOGN_CORE_CMAP_CMAP_OPS_H_
#define CGOGN_CORE_CMAP_CMAP_OPS_H_

#include <cgogn/core/types/cmap/cmap_base.h>
#include <cgogn/core/types/cmap/cph3.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CMAP>
// Dart add_dart(CMAP& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

inline Dart add_dart(CMapBase& m)
{
	uint32 index = m.darts_.new_index();
	Dart d(index);
	for (auto rel : m.relations_)
		(*rel)[d.index] = d;
	for (auto emb : m.cells_indices_)
		if (emb)
			(*emb)[d.index] = INVALID_INDEX;
	return d;
}

//////////
// CPH3 //
//////////

// inline Dart add_dart(CPH3& m)
// {
// 	Dart d = add_dart(static_cast<CPH3::CMAP&>(m));
// 	m.nb_darts_per_level_[m.current_level_]++;
// 	m.set_edge_id(d, 0u);
// 	m.set_face_id(d, 0u);
// 	m.set_dart_level(d, m.current_level_);

// 	// update max level if needed
// 	if (m.current_level_ > m.maximum_level_)
// 		m.maximum_level_ = m.current_level_;

// 	return d;
// }

/*****************************************************************************/

// template <typename CMAP>
// void remove_dart(CMAP& m, Dart d);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

inline void remove_dart(CMapBase& m, Dart d)
{
	for (uint32 orbit = 0; orbit < m.attribute_containers_.size(); ++orbit)
	{
		if (m.cells_indices_[orbit])
		{
			uint32 index = (*m.cells_indices_[orbit])[d.index];
			if (index != INVALID_INDEX)
				m.attribute_containers_[orbit].unref_index(index);
		}
	}
	m.darts_.release_index(d.index);
}

/*****************************************************************************/

// template <typename CMAP>
// void set_boundary(const CMAP& m, Dart d, bool b)

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

inline void set_boundary(const CMapBase& m, Dart d, bool b)
{
	(*m.boundary_marker_)[d.index] = b ? 1u : 0u;
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void set_index(MESH& m, Dart d, uint32 index);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
void set_index(CMapBase& m, Dart d, uint32 index)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	const uint32 old = (*m.cells_indices_[orbit])[d.index];
	// ref_index() is done before unref_index() to avoid deleting the index if old == index
	m.attribute_containers_[orbit].ref_index(index); // ref the new index
	if (old != INVALID_INDEX)
		m.attribute_containers_[orbit].unref_index(old); // unref the old index
	(*m.cells_indices_[orbit])[d.index] = index;		 // affect the index to the dart
}

/*****************************************************************************/

// template <typename CELL, typename CMAP>
// void set_index(CMAP& m, CELL c, uint32 index);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename CMAP>
void set_index(CMAP& m, CELL c, uint32 index)
{
	static_assert(is_in_tuple<CELL, typename CMAP::Cells>::value, "CELL not supported in this CMAP");
	foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
		set_index<CELL>(m, d, index);
		return true;
	});
}

/*****************************************************************************/

// template <typename CELL, typename CMAP>
// void copy_index(CMAP& m, Dart dest, Dart src);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename CMAP>
inline void copy_index(CMAP& m, Dart dest, Dart src)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	set_index<CELL>(m, dest, index_of(m, CELL(src)));
}

} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_OPS_H_
