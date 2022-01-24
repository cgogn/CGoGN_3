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

#ifndef CGOGN_CORE_FUNCTIONS_CELLS_H_
#define CGOGN_CORE_FUNCTIONS_CELLS_H_

#include <cgogn/core/types/cmap/dart_marker.h>

#include <cgogn/core/types/cmap/cmap_base.h>
#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>

#include <sstream>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL>
// bool is_indexed(CMapBase& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
bool is_indexed(const CMapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}

inline bool is_indexed(const CMapBase& m, Orbit orbit)
{
	cgogn_message_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}

////////////////////
// IncidenceGraph //
////////////////////

template <typename CELL>
bool is_indexed(const IncidenceGraph& /*m*/)
{
	return true;
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 maximum_index(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
uint32 maximum_index(const CMapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access a cell index of an unindexed cell type");
	return m.attribute_containers_[CELL::ORBIT].maximum_index();
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 index_of(MESH& m, CELL c);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
uint32 index_of(const CMapBase& m, CELL c)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	return (*m.cells_indices_[orbit])[c.dart.index];
}

////////////////////
// IncidenceGraph //
////////////////////

template <typename CELL>
uint32 index_of(const IncidenceGraph& /*m*/, CELL c)
{
	return c.index_;
}

//////////
// CPH3 //
//////////

template <typename MRMAP, typename CELL>
inline auto index_of(const MRMAP& m, CELL c) -> std::enable_if_t<std::is_convertible_v<MRMAP&, CPH3&>, uint32>
{
	static const Orbit orbit = CELL::ORBIT;

	if constexpr (orbit == CPH3::CMAP::Edge::ORBIT)
		c.dart = m.edge_youngest_dart(c.dart);
	if constexpr (orbit == CPH3::CMAP::Face::ORBIT)
		c.dart = m.face_youngest_dart(c.dart);
	if constexpr (orbit == CPH3::CMAP::Volume::ORBIT)
		c.dart = m.volume_youngest_dart(c.dart);

	return index_of(static_cast<const CPH3::CMAP&>(m), c);
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// CELL of_index(MESH& m, uint32 i);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
CELL of_index(const CMapBase& m, uint32 i)
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

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 new_index(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
uint32 new_index(const CMapBase& m)
{
	return m.attribute_containers_[CELL::ORBIT].new_index();
}

////////////////////
// IncidenceGraph //
////////////////////

template <typename CELL>
uint32 new_index(const IncidenceGraph& ig)
{
	uint32 id = ig.attribute_containers_[CELL::CELL_INDEX].new_index();
	// (*ig.cells_indices_[CELL::CELL_INDEX])[id] = id;
	return id;
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// bool init_cells_indexing(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
void init_cells_indexing(CMapBase& m)
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

inline void init_cells_indexing(CMapBase& m, Orbit orbit)
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

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 set_index(MESH& m, CELL c);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH>
auto set_index(MESH& m, CELL c, uint32 index) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
		set_index<CELL>(m, d, index);
		return true;
	});
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void index_cells(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH>
auto index_cells(MESH& m) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		init_cells_indexing<CELL>(m);

	CMapBase& base = static_cast<CMapBase&>(m);
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

/////////////
// GENERIC //
/////////////

// template <typename CELL, typename MESH>
// void index_cells(MESH& m)
// {
// 	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
// 	if (!is_indexed<CELL>(m))
// 		init_cells_indexing<CELL>(m);

// 	foreach_cell(m, [&](CELL c) -> bool {
// 		if (index_of(m, c) == INVALID_INDEX)
// 			set_index(m, c, new_index<CELL>(m));
// 		return true;
// 	});
// }

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_CELLS_H_
