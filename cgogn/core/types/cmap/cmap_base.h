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

	inline uint32 nb_darts() const
	{
		return topology_.nb_elements();
	}

	inline void set_boundary(Dart d, bool b)
	{
		(*boundary_marker_)[d.index] = b ? 1u : 0u;
	}

	inline bool is_boundary(Dart d) const
	{
		return (*boundary_marker_)[d.index] != 0u;
	}

	template <typename CELL>
	inline bool is_indexed() const
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
		return cells_indices_[orbit] != nullptr;
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

	template <typename CELL>
	inline void init_cells_indexing()
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
		if (!is_indexed<CELL>())
		{
			std::ostringstream oss;
			oss << "__index_" << orbit_name(orbit);
			cells_indices_[orbit] = topology_.add_attribute<uint32>(oss.str());
			cells_indices_[orbit]->fill(INVALID_INDEX);
		}
	}

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

	template <typename FUNC>
	void foreach_dart(const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		for (uint32 i = topology_.first_index(), last_index = topology_.last_index(); i < last_index; i = topology_.next_index(i))
			if (!f(Dart(i)))
				break;
	}

	inline Dart begin() const { return Dart(topology_.first_index()); }
	inline Dart end() const { return Dart(topology_.last_index()); }
	inline Dart next(Dart d) const { return Dart(topology_.next_index(d.index)); }
};


} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
