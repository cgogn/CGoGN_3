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

#ifndef CGOGN_CORE_TYPES_CMAP_CELL_MARKER_H_
#define CGOGN_CORE_TYPES_CMAP_CELL_MARKER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap_base.h>

namespace cgogn
{

template <typename CELL>
class CellMarker
{
	static const Orbit orbit = CELL::ORBIT;

protected:

	const CMapBase& map_;
	CMapBase::Attribute<uint8>* mark_attribute_;

public:

	CellMarker(const CMapBase& map) : map_(map)
	{
		mark_attribute_ = map.attribute_containers_[orbit].add_mark_attribute();
	}

	virtual ~CellMarker()
	{
		delete mark_attribute_;
	}

	inline void mark(CELL c) { (*mark_attribute_)[map_.embedding(c)] = 1u; }
	inline void unmark(CELL c) { (*mark_attribute_)[map_.embedding(c)] = 0u; }

	inline bool is_marked(CELL c) const
	{
		return (*mark_attribute_)[map_.embedding(c)] != 0u;
	}

	inline void unmark_all()
	{
		std::fill(mark_attribute_->begin(), mark_attribute_->end(), 0u);
	}
};

template <typename CELL>
class CellMarkerStore : public CellMarker<CELL>
{
	std::vector<uint32> marked_cells_;

public:

	inline CellMarkerStore(const CMapBase& map) : CellMarker<CELL>(map)
	{}

	~CellMarkerStore() override
	{}

	inline void mark(CELL c)
	{
		if (!is_marked(c))
		{
			CellMarker<CELL>::mark(c);
			marked_cells_.push_back(this->map_.embedding(c));
		}
	}

	inline void unmark(CELL c)
	{
		auto it = std::find(marked_cells_.begin(), marked_cells_.end(), this->map_.embedding(c));
		if (it !=  marked_cells_.end())
		{
			CellMarker<CELL>::unmark(c);
			std::swap(*it, marked_cells_.back());
			marked_cells_.pop_back();
		}
	}

	inline void unmark_all()
	{
		for (uint32 i : marked_cells_)
			return (*this->mark_attribute_)[i] != 0u;
		marked_cells_.clear();
	}

	inline const std::vector<uint32>& marked_cells() const
	{
		return marked_cells_;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CELL_MARKER_H_
