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

#include <cgogn/core/types/container/chunk_array_container.h>
#include <cgogn/core/types/cmap/cell.h>

#include <cgogn/core/utils/type_traits.h>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/utils/tuples.h>

#include <array>
#include <sstream>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMapBase
{
	using AttributeContainer = ChunkArrayContainer<32>;

	template <typename T>
	using Attribute = AttributeContainer::ChunkArray<T>;
	template <typename T>
	using AttributePtr = AttributeContainer::ChunkArrayPtr<T>;

	using AttributeGen = AttributeContainer::ChunkArrayGen;
	using AttributeGenPtr = AttributeContainer::ChunkArrayGenPtr;

	// Dart container
	mutable AttributeContainer topology_;
	// shortcuts to relations Dart attributes
	std::vector<AttributePtr<Dart>> relations_;
	// shortcuts to embedding indices Dart attributes
	std::array<AttributePtr<uint32>, NB_ORBITS> embeddings_;
	// shortcut to boundary marker Dart attribute
	AttributePtr<uint8> boundary_marker_;

	// Cells attributes containers
	mutable std::array<AttributeContainer, NB_ORBITS> attribute_containers_;

	CMapBase();
	virtual ~CMapBase();

protected:

	AttributePtr<Dart> add_relation(const std::string& name)
	{
		AttributePtr<Dart> rel = topology_.add_chunk_array<Dart>(name);
		relations_.push_back(rel);
		return rel;
	}

public:

	uint32 nb_darts()
	{
		return topology_.size();
	}

	void set_boundary(Dart d, bool b)
	{
		(*boundary_marker_)[d.index] = b ? 1u : 0u;
	}

	bool is_boundary(Dart d) const
	{
		return (*boundary_marker_)[d.index] != 0u;
	}

	template <typename CELL>
	bool is_embedded() const
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
		return embeddings_[orbit] != nullptr;
	}

	template <typename CELL>
	uint32 embedding(CELL c) const
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
		return (*embeddings_[orbit])[c.dart.index];
	}

	template <typename CELL>
	void set_embedding(Dart d, uint32 emb)
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
		(*embeddings_[orbit])[d.index] = emb;
	}

	template <typename CELL>
	void copy_embedding(Dart dest, Dart src)
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
		set_embedding<CELL>(dest, embedding(CELL(src)));
	}

	template <typename CELL>
	void create_embedding()
	{
		static const Orbit orbit = CELL::ORBIT;
		static_assert (orbit < NB_ORBITS, "Unknown orbit parameter");
		std::ostringstream oss;
		oss << "emb_" << orbit_name(orbit);
		AttributePtr<uint32> emb = topology_.add_chunk_array<uint32>(oss.str());
		embeddings_[orbit] = emb;
		for (uint32& i : *emb)
			i = INVALID_INDEX;
	}

	Dart add_dart()
	{
		uint32 index = topology_.get_index();
		Dart d(index);
		for (auto rel : relations_)
			(*rel)[d.index] = d;
		return d;
	}

	template <typename FUNC>
	void foreach_dart(const FUNC& f) const
	{
		static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		for (uint32 i = 0; i < topology_.size(); ++i)
			if (!f(Dart(i)))
				break;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
