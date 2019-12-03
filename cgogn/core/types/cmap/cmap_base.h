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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_H_
#define CGOGN_CORE_CMAP_CMAP_BASE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/vector.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/cmap/cell.h>

#include <cgogn/core/utils/type_traits.h>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/assert.h>

#include <array>
#include <sstream>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMapBase
{
	// using AttributeContainer = AttributeContainerT<Vector>;
	using AttributeContainer = AttributeContainerT<ChunkArray>;

	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	using AttributeGen = AttributeContainer::AttributeGen;
	using MarkAttribute = AttributeContainer::MarkAttribute;

	// Dart container
	mutable AttributeContainer topology_;
	// shortcuts to relations Dart attributes
	std::vector<std::shared_ptr<Attribute<Dart>>> relations_;
	// shortcuts to embedding indices Dart attributes
	std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS> cells_indices_;
	// shortcut to boundary marker Dart attribute
    MarkAttribute* boundary_marker_;

	// Cells attributes containers
	mutable std::array<AttributeContainer, NB_ORBITS> attribute_containers_;

	CMapBase();
	virtual ~CMapBase();

public:

	std::shared_ptr<Attribute<Dart>> add_relation(const std::string& name)
	{
		return relations_.emplace_back(topology_.add_attribute<Dart>(name));
	}

    std::shared_ptr<Attribute<Dart>> add_or_get_relation(const std::string& name)
    {
        auto rel = topology_.add_attribute<Dart>(name);
        if(rel == nullptr)
            return topology_.get_attribute<Dart>(name);
        return relations_.emplace_back(rel);
    }

	// template <typename CELL>
	// inline void unset_index(Dart d)
	// {
	// 	static const Orbit orbit = CELL::ORBIT;
	// 	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	// 	cgogn_message_assert(is_indexed<CELL>(), "Trying to access the cell index of an unindexed cell type");
	// 	const uint32 old = (*cells_indices_[orbit])[d.index];
	// 	if (old != INVALID_INDEX)
	// 		attribute_containers_[orbit].unref_index(old);	// unref the old emb
	// 	(*cells_indices_[orbit])[d.index] = INVALID_INDEX;	// affect the index to the dart
	// }

	inline Dart add_dart()
	{
		uint32 index = topology_.new_index();
		Dart d(index);
		for (auto rel : relations_)
			(*rel)[d.index] = d;
		for (auto emb : cells_indices_)
			if (emb)
				(*emb)[d.index] = INVALID_INDEX;
		return d;
	}

	inline void remove_dart(Dart d)
	{
		for (uint32 orbit = 0; orbit < NB_ORBITS; ++orbit)
		{
			if (cells_indices_[orbit])
			{
				uint32 index = (*cells_indices_[orbit])[d.index];
				if (index != INVALID_INDEX)
					attribute_containers_[orbit].unref_index(index);
			}
		}
		topology_.release_index(d.index);
	}
};


} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
