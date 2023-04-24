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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_H_
#define CGOGN_CORE_CMAP_CMAP_BASE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/container/vector.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/types/cmap/cell.h>
#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>


#include <any>
#include <array>
#include <unordered_map>

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

	/*************************************************************************/
	// Map-wise attributes container
	/*************************************************************************/
	std::unordered_map<std::string, std::any> attributes_;

	/*************************************************************************/
	// Dart attributes container
	/*************************************************************************/
	mutable AttributeContainer darts_;

	// shortcuts to topological relations attributes
	std::vector<std::shared_ptr<Attribute<Dart>>> relations_;
	// shortcuts to cells indices attributes
	std::array<std::shared_ptr<Attribute<uint32>>, NB_ORBITS> cells_indices_;

	// shortcut to boundary marker attribute
	MarkAttribute* boundary_marker_;

	/*************************************************************************/
	// Cells attributes containers
	/*************************************************************************/
	mutable std::array<AttributeContainer, NB_ORBITS> attribute_containers_;

	enum TraversalPolicy
	{
		AUTO,
		DART_MARKING
	};

	CMapBase();
	~CMapBase();

	// Map-wise attributes
	template <typename T>
	T& get_attribute(const std::string& name)
	{
		auto [it, inserted] = attributes_.try_emplace(name, T());
		return std::any_cast<T&>(it->second);
	}

	inline std::shared_ptr<Attribute<Dart>> add_relation(const std::string& name)
	{
		return relations_.emplace_back(darts_.add_attribute<Dart>(name));
	}

	inline Dart begin() const
	{
		return Dart(darts_.first_index());
	}
	
	inline Dart end() const
	{
		return Dart(darts_.last_index());
	}

	inline Dart next(Dart d) const
	{
		return Dart(darts_.next_index(d.index));
	}

};

void CGOGN_CORE_EXPORT dump_map_darts(const CMapBase& m);

Dart CGOGN_CORE_EXPORT add_dart(CMapBase& m);

void CGOGN_CORE_EXPORT remove_dart(CMapBase& m, Dart d);


inline uint32 nb_darts(const CMapBase& m)
{
	return m.darts_.nb_elements();
}

inline void set_boundary(const CMapBase& m, Dart d, bool b)
{
	(*m.boundary_marker_)[d.index] = b ? 1u : 0u;
}


} // namespace cgogn

#include <cgogn/core/types/cmap/cmap_base_index.hpp>
#include <cgogn/core/types/cmap/cmap_base_attribute.hpp>
#include <cgogn/core/types/cmap/cmap_base_orbit_traversals.hpp>
#include <cgogn/core/types/cmap/cmap_base_local_traversals.hpp>
#include <cgogn/core/types/cmap/cmap_base_global_traversals.hpp>
#include <cgogn/core/types/cmap/cmap_base_functions.hpp>

#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
