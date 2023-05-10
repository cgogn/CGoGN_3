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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_IND_H_
#define CGOGN_CORE_CMAP_CMAP_BASE_IND_H_


#include <sstream>
#include <iostream>
#include <iomanip>

namespace cgogn
{



template <typename CELL>
bool is_indexed(const MapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}


inline bool is_indexed(const MapBase& m, Orbit orbit)
{
	cgogn_message_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}


template <typename CELL>
uint32 maximum_index(const MapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access a cell index of an unindexed cell type");
	return m.attribute_containers_[CELL::ORBIT].maximum_index();
}



template <typename CELL>
uint32 index_of(const MapBase& m, CELL c)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	return (*m.cells_indices_[orbit])[c.dart.index];
}

inline bool is_boundary(const MapBase& m, Dart d)
{
	return (*m.boundary_marker_)[d.index] != 0u;
}


template <typename CELL>
CELL of_index(const MapBase& m, uint32 i)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		if (!is_boundary(m, d))
		{
			const CELL c(d);
			if (index_of(m, c) == i)
				return c;
		}
	}
	return CELL();
}



template <typename CELL>
uint32 new_index(const MapBase& m)
{
	return m.attribute_containers_[CELL::ORBIT].new_index();
}


template <typename CELL>
void set_index(MapBase& m, Dart d, uint32 index)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");

	const uint32 old = (*m.cells_indices_[orbit])[d.index];
	// ref_index() is done before unref_index() to avoid deleting the index if old == index
	if (index != INVALID_INDEX)
		m.attribute_containers_[orbit].ref_index(index); // ref the new index
	if (old != INVALID_INDEX)
		m.attribute_containers_[orbit].unref_index(old); // unref the old index
	(*m.cells_indices_[orbit])[d.index] = index;		 // affect the index to the dart
}

template <typename CELL>
void init_cells_indexing(MapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	if (!is_indexed<CELL>(m))
	{
		std::ostringstream oss;
		oss << "__index_" << orbit_name(orbit);
		m.cells_indices_[orbit] = m.darts_.add_attribute<uint32>(oss.str());
		m.cells_indices_[orbit]->fill(INVALID_INDEX);
	}
}


inline void init_cells_indexing(MapBase& m, Orbit orbit)
{
	cgogn_message_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	if (!is_indexed(m, orbit))
	{
		std::ostringstream oss;
		oss << "__index_" << orbit_name(orbit);
		m.cells_indices_[orbit] = m.darts_.add_attribute<uint32>(oss.str());
		m.cells_indices_[orbit]->fill(INVALID_INDEX);
	}
}



template <typename CELL, typename MESH>
auto set_index(MESH& m, CELL c, uint32 index) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
		set_index<CELL>(m, d, index);
		return true;
	});
}


template <typename CELL, typename MESH>
auto copy_index(MESH& m, Dart dest, Dart src) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	set_index<CELL>(m, dest, index_of(m, CELL(src)));
}






template <typename CELL, typename MESH>
auto index_cells(MESH& m) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		init_cells_indexing<CELL>(m);

	MapBase& base = static_cast<MapBase&>(m);
	DartMarker dm(m);
	for (Dart d = base.begin(), end = base.end(); d != end; d = base.next(d))
	{
		if (!is_boundary(m, d) && !dm.is_marked(d))
		{
			const CELL c(d);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				dm.mark(d);
				return true;
			});

			if (index_of(m, c) == INVALID_INDEX)
				set_index(m, c, new_index<CELL>(m));
		}
	}
}


} // namespace cgogn


#endif // CGOGN_CORE_CMAP_CMAP_BASE_IND_H_
