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

#ifndef CGOGN_CORE_TYPES_MARKER_H_
#define CGOGN_CORE_TYPES_MARKER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/functions/cells.h>
#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// typename mesh_traits<MESH>::MarkAttribute* get_mark_attribute(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH>
auto get_mark_attribute(const MESH& m)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, typename mesh_traits<MESH>::MarkAttribute*>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		index_cells<CELL>(const_cast<MESH&>(m));
	const CMapBase& mb = static_cast<const CMapBase&>(m);
	return mb.attribute_containers_[CELL::ORBIT].get_mark_attribute();
}

////////////////////
// IncidenceGraph //
////////////////////
template <typename CELL>
auto get_mark_attribute(const IncidenceGraph& ig)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this MESH");
	return ig.attribute_containers_[CELL::CELL_INDEX].get_mark_attribute();
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void release_mark_attribute(const MESH& m, typename mesh_traits<MESH>::MarkAttribute* attribute);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
void release_mark_attribute(const CMapBase& m, CMapBase::MarkAttribute* attribute)
{
	return m.attribute_containers_[CELL::ORBIT].release_mark_attribute(attribute);
}

////////////////////
// IncidenceGraph //
////////////////////

template <typename CELL>
void release_mark_attribute(const IncidenceGraph& ig, IncidenceGraph::MarkAttribute* attribute)
{
	return ig.attribute_containers_[CELL::CELL_INDEX].release_mark_attribute(attribute);
}

/*****************************************************************************/

template <typename MESH, typename CELL>
class CellMarker
{
private:
	const MESH& mesh_;
	typename mesh_traits<MESH>::MarkAttribute* mark_attribute_;

public:
	CellMarker(const MESH& mesh) : mesh_(mesh)
	{
		mark_attribute_ = get_mark_attribute<CELL>(mesh_);
	}

	~CellMarker()
	{
		unmark_all();
		release_mark_attribute<CELL>(mesh_, mark_attribute_);
	}

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(CellMarker);

	inline void mark(CELL c)
	{
		(*mark_attribute_)[index_of(mesh_, c)] = 1u;
	}
	inline void unmark(CELL c)
	{
		(*mark_attribute_)[index_of(mesh_, c)] = 0u;
	}

	inline bool is_marked(CELL c) const
	{
		return (*mark_attribute_)[index_of(mesh_, c)] != 0u;
	}

	inline void unmark_all()
	{
		mark_attribute_->fill(0u);
	}
};

template <typename MESH, typename CELL>
class CellMarkerStore
{
private:
	const MESH& mesh_;
	typename mesh_traits<MESH>::MarkAttribute* mark_attribute_;
	std::vector<uint32> marked_cells_;

public:
	inline CellMarkerStore(const MESH& mesh) : mesh_(mesh)
	{
		mark_attribute_ = get_mark_attribute<CELL>(mesh_);
		marked_cells_.reserve(512u);
	}

	~CellMarkerStore()
	{
		unmark_all();
		release_mark_attribute<CELL>(mesh_, mark_attribute_);
	}

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(CellMarkerStore);

	inline void mark(CELL c)
	{
		if (!is_marked(c))
		{
			uint32 index = index_of(mesh_, c);
			(*mark_attribute_)[index] = 1u;
			marked_cells_.push_back(index);
		}
	}

	inline void unmark(CELL c)
	{
		uint32 index = index_of(mesh_, c);
		auto it = std::find(marked_cells_.begin(), marked_cells_.end(), index);
		if (it != marked_cells_.end())
		{
			(*mark_attribute_)[index] = 0u;
			std::swap(*it, marked_cells_.back());
			marked_cells_.pop_back();
		}
	}

	inline bool is_marked(CELL c) const
	{
		return (*mark_attribute_)[index_of(mesh_, c)] != 0u;
	}

	inline void unmark_all()
	{
		for (uint32 i : marked_cells_)
			(*mark_attribute_)[i] = 0u;
		marked_cells_.clear();
	}

	inline const std::vector<uint32>& marked_cells() const
	{
		return marked_cells_;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MARKER_H_
