/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_CORE_TYPES_CELLS_SET_H_
#define CGOGN_CORE_TYPES_CELLS_SET_H_

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/utils/tuples.h>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

template <typename MESH, typename CELL>
struct CellsSet
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");

public:
	CellsSet(const MESH& m, const std::string& name) : m_(m), marker_(m), name_(name)
	{
	}

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(CellsSet);

	const std::string& name() const
	{
		return name_;
	}

	uint32 size() const
	{
		return uint32(cells_.size());
	}

	inline void select(CELL c)
	{
		if (!marker_.is_marked(c))
		{
			marker_.mark(c);
			cells_.emplace(index_of(m_, c), c);
		}
	}

	inline void unselect(CELL c)
	{
		if (marker_.is_marked(c))
		{
			auto it = cells_.find(index_of(m_, c));
			if (it != cells_.end())
			{
				marker_.unmark(it->second);
				cells_.erase(it);
			}
		}
	}

	inline bool contains(CELL c) const
	{
		return marker_.is_marked(c);
	}

	inline void rebuild()
	{
		cells_.clear();
		cgogn::foreach_cell(m_, [&](CELL c) -> bool {
			if (marker_.is_marked(c))
				cells_.emplace(index_of(m_, c), c);
			return true;
		});
	}

	inline void clear()
	{
		cells_.clear();
		marker_.unmark_all();
	}

	template <typename FUNC>
	void foreach_cell(const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, CELL>::value, "Wrong function parameter type");
		for (auto& [index, cell] : cells_)
			f(cell);
	}

	template <typename FUNC>
	void foreach_cell_index(const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, uint32>::value, "Wrong function parameter type");
		for (auto& [index, cell] : cells_)
			f(index);
	}

	class const_iterator
	{
		const CellsSet<MESH, CELL>* cs_;
		typename std::unordered_map<uint32, CELL>::const_iterator map_it_;

	public:
		inline const_iterator(const CellsSet<MESH, CELL>* cs,
							  typename std::unordered_map<uint32, CELL>::const_iterator&& map_it)
			: cs_(cs), map_it_(map_it)
		{
		}
		inline const_iterator(const const_iterator& it) : cs_(it.cs_), map_it_(it.map_it_)
		{
		}
		inline const_iterator& operator=(const const_iterator& it)
		{
			cs_ = it.cs_;
			map_it_ = it.map_it_;
			return *this;
		}
		inline bool operator!=(const_iterator it) const
		{
			cgogn_assert(cs_ == it.cs_);
			return map_it_ != it.map_it_;
		}
		inline const_iterator& operator++()
		{
			map_it_++;
			return *this;
		}
		inline CELL operator*() const
		{
			return map_it_->second;
		}
	};
	inline const_iterator begin() const
	{
		return const_iterator(this, cells_.begin());
	}
	inline const_iterator end() const
	{
		return const_iterator(this, cells_.end());
	}

private:
	const MESH& m_;
	CellMarker<MESH, CELL> marker_;
	std::unordered_map<uint32, CELL> cells_;
	std::string name_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CELLS_SET_H_
