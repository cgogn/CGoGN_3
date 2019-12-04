/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/dart.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/functions/cmapbase_infos.h>
#include <functional>

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
inline
bool is_indexed(const CMapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
bool is_indexed(const MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return is_indexed<CELL>(m.mesh());
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 index_of(MESH& m, CELL c);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
inline uint32 index_of(const CMapBase& m,CELL c)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	return (*m.cells_indices_[orbit])[c.dart.index];
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
uint32 index_of(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return index_of(m.mesh(), c);
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 new_index(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
inline
auto new_index(const CMapBase& m)
{
	return m.attribute_containers_[CELL::ORBIT].new_index();
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
uint32 new_index(const MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return new_index<CELL>(m.mesh());
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void set_index(MESH& m, CELL c, uint32 index);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
inline void set_index(CMapBase& m, Dart d, uint32 index)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	const uint32 old = (*m.cells_indices_[orbit])[d.index];
	// ref_index() is done before unref_index() to avoid deleting the index if old == index
	m.attribute_containers_[orbit].ref_index(index);		// ref the new index
	if (old != INVALID_INDEX)
		m.attribute_containers_[orbit].unref_index(old);	// unref the old index
	(*m.cells_indices_[orbit])[d.index] = index;			// affect the index to the dart
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
void set_index(MESH& m, Dart d, uint32 index)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	set_index<CELL>(m.mesh(),d, index);
}


template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
void set_index(MESH& m, CELL c, uint32 index)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_dart_of_orbit(m,c, [&] (Dart d) -> bool
	{
		set_index<CELL>(m.mesh(),d, index);
		return true;
	});
}

/*****************************************************************************/

// template <typename CELL>
// bool is_indexed(CMapBase& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
inline
void init_cells_indexing(CMapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	if (!is_indexed<CELL>(m))
	{
		std::ostringstream oss;
		oss << "__index_" << orbit_name(orbit);
		m.cells_indices_[orbit] = m.topology_.add_attribute<uint32>(oss.str());
		m.cells_indices_[orbit]->fill(INVALID_INDEX);
	}
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
void init_cells_indexing(MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	init_cells_indexing<CELL>(m.mesh());
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void index_cells(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////


//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
void index_cells(MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if(!is_indexed<CELL>(m))
	{
		init_cells_indexing<CELL>(m);
	}
	
	DartMarker dm(m);
	foreach_dart(m,[&] (Dart d) -> bool
	{
		if (!is_boundary(m,d) && !dm.is_marked(d))
		{
			const CELL c(d);
			foreach_dart_of_orbit(m,c, [&] (Dart d) -> bool { dm.mark(d); return true; });
			if (index_of(m, c) == INVALID_INDEX)
				set_index<CELL>(m, c, new_index<CELL>(m));
			return true;
		}
		return true;
	});
}


//////////////
// CMapBase //
//////////////

template <typename CELL>
inline void copy_index(CMapBase& m,Dart dest, Dart src)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	set_index<CELL>(m,dest, index_of(m,CELL(src)));
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL,typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
void copy_index(MESH& m,Dart dest, Dart src)
{
	return copy_index<CELL>(m.mesh(),dest,src);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_CELLS_H_
