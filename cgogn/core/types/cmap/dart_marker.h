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

#ifndef CGOGN_CORE_TYPES_CMAP_DART_MARKER_H_
#define CGOGN_CORE_TYPES_CMAP_DART_MARKER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap_base.h>

namespace cgogn
{

class CGOGN_CORE_EXPORT DartMarker
{
	const CMapBase& map_;
	CMapBase::MarkAttributePtr mark_attribute_;

public:

	DartMarker(const CMapBase& map);
	virtual ~DartMarker();

	inline void mark(Dart d) { (*mark_attribute_)[d.index] = 1u; }
	inline void unmark(Dart d) { (*mark_attribute_)[d.index] = 0u; }

	inline bool is_marked(Dart d) const
	{
		return (*mark_attribute_)[d.index] != 0u;
	}

	virtual inline void unmark_all()
	{
		mark_attribute_->fill(0u);
	}
};

class CGOGN_CORE_EXPORT DartMarkerStore : public DartMarker
{
	std::vector<Dart> marked_darts_;

public:

	DartMarkerStore(const CMapBase& map);
	~DartMarkerStore() override;

	inline void mark(Dart d)
	{
		if (!is_marked(d))
		{
			DartMarker::mark(d);
			marked_darts_.push_back(d);
		}
	}

	inline void unmark(Dart d)
	{
		auto it = std::find(marked_darts_.begin(), marked_darts_.end(), d);
		if (it != marked_darts_.end())
		{
			DartMarker::unmark(d);
			std::swap(*it, marked_darts_.back());
			marked_darts_.pop_back();
		}
	}

	inline void unmark_all() override
	{
		for (Dart d : marked_darts_)
			DartMarker::unmark(d);
		marked_darts_.clear();
	}

	inline const std::vector<Dart>& marked_darts() const
	{
		return marked_darts_;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_DART_MARKER_H_
