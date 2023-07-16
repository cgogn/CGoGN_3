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

#ifndef CGOGN_CORE_TYPES_MAPS_DART_MARKER_H_
#define CGOGN_CORE_TYPES_MAPS_DART_MARKER_H_

namespace cgogn
{

template <typename MESH>
struct mesh_traits;

/*****************************************************************************/

template <typename MAP>
class DartMarker
{
private:
	const MAP& map_;
	typename mesh_traits<MAP>::MarkAttribute* mark_attribute_;

public:
	DartMarker(const MAP& map) : map_(map)
	{
		mark_attribute_ = get_dart_mark_attribute(map_);
	}

	~DartMarker()
	{
		unmark_all();
		release_dart_mark_attribute(map_, mark_attribute_);
	}

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(DartMarker);

	inline void mark(Dart d)
	{
		(*mark_attribute_)[d.index] = 1u;
	}
	inline void unmark(Dart d)
	{
		(*mark_attribute_)[d.index] = 0u;
	}

	inline bool is_marked(Dart d) const
	{
		return (*mark_attribute_)[d.index] != 0u;
	}

	inline void unmark_all()
	{
		mark_attribute_->fill(0u);
	}
};

template <typename MAP>
class DartMarkerStore
{
private:
	const MAP& map_;
	typename mesh_traits<MAP>::MarkAttribute* mark_attribute_;
	std::vector<Dart> marked_darts_;

public:
	DartMarkerStore(const MAP& map) : map_(map)
	{
		mark_attribute_ = get_dart_mark_attribute(map_);
		marked_darts_.reserve(512u);
	}

	~DartMarkerStore()
	{
		unmark_all();
		release_dart_mark_attribute(map_, mark_attribute_);
	}

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(DartMarkerStore);

	inline void mark(Dart d)
	{
		if (!is_marked(d))
		{
			(*mark_attribute_)[d.index] = 1u;
			marked_darts_.push_back(d);
		}
	}

	inline void unmark(Dart d)
	{
		auto it = std::find(marked_darts_.begin(), marked_darts_.end(), d);
		if (it != marked_darts_.end())
		{
			(*mark_attribute_)[d.index] = 0u;
			std::swap(*it, marked_darts_.back());
			marked_darts_.pop_back();
		}
	}

	inline bool is_marked(Dart d) const
	{
		return (*mark_attribute_)[d.index] != 0u;
	}

	inline void unmark_all()
	{
		for (Dart d : marked_darts_)
			(*mark_attribute_)[d.index] = 0u;
		marked_darts_.clear();
	}

	inline const std::vector<Dart>& marked_darts() const
	{
		return marked_darts_;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MAPS_DART_MARKER_H_
