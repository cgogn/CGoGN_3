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

#ifndef CGOGN_CORE_TYPES_MARKER_H_
#define CGOGN_CORE_TYPES_MARKER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/cmap/cmap_ops.h>

#include <cgogn/core/functions/traversals/global.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// typename mesh_traits<MESH>::template MarkAttributePtr<T> add_mark_attribute(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::MarkAttributePtr
add_mark_attribute(const MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!m.template is_embedded<CELL>())
	{
        MESH& mesh = const_cast<MESH&>(m);
		mesh.template init_embedding<CELL>();
		foreach_cell(mesh, [&] (CELL c) -> bool
		{
			create_embedding(mesh, c);
			return true;
		}, true);
	}
	return m.attribute_containers_[CELL::ORBIT].add_mark_attribute();
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::MarkAttributePtr
add_mark_attribute(const MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return add_mark_attribute<CELL>(m.mesh());
}

/*****************************************************************************/
/*                                   CellMarker                              */
/*****************************************************************************/

template <typename MESH, typename CELL>
class CellMarker
{
protected:

	const MESH& mesh_;
	typename mesh_traits<MESH>::MarkAttributePtr mark_attribute_;

public:

	CellMarker(const MESH& mesh) : mesh_(mesh)
	{
		mark_attribute_ = add_mark_attribute<CELL>(mesh);
	}

	virtual ~CellMarker()
	{}

	inline void mark(CELL c) { (*mark_attribute_)[index_of(mesh_, c)] = 1u; }
	inline void unmark(CELL c) { (*mark_attribute_)[index_of(mesh_, c)] = 0u; }

	inline bool is_marked(CELL c) const
	{
		return (*mark_attribute_)[index_of(mesh_, c)] != 0u;
	}

	inline void unmark_all()
	{
		std::fill(mark_attribute_->begin(), mark_attribute_->end(), 0u);
	}
};

template <typename MESH, typename CELL>
class CellMarkerStore : public CellMarker<MESH, CELL>
{
	std::vector<uint32> marked_cells_;

public:

	inline CellMarkerStore(const MESH& mesh) : CellMarker<MESH, CELL>(mesh)
	{}

	~CellMarkerStore() override
	{}

	inline void mark(CELL c)
	{
		if (!is_marked(c))
		{
			CellMarker<MESH, CELL>::mark(c);
			marked_cells_.push_back(index_of(this->mesh_, c));
		}
	}

	inline void unmark(CELL c)
	{
		auto it = std::find(marked_cells_.begin(), marked_cells_.end(), index_of(this->mesh_, c));
		if (it !=  marked_cells_.end())
		{
			CellMarker<MESH, CELL>::unmark(c);
			std::swap(*it, marked_cells_.back());
			marked_cells_.pop_back();
		}
	}

	inline void unmark_all()
	{
		for (uint32 i : marked_cells_)
            (*this->mark_attribute_)[i] = 0u;
		marked_cells_.clear();
	}

	inline const std::vector<uint32>& marked_cells() const
	{
		return marked_cells_;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MARKER_H_
