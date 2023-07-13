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

#ifndef CGOGN_CORE_TYPES_MAPS_MAP_BASE_H_
#define CGOGN_CORE_TYPES_MAPS_MAP_BASE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/container/vector.h>
#include <cgogn/core/types/maps/cell.h>
#include <cgogn/core/types/maps/dart_marker.h>
#include <cgogn/core/utils/thread.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>

#include <any>
#include <array>
#include <iostream>
#include <sstream>
#include <unordered_map>

namespace cgogn
{

template <typename MESH>
struct mesh_traits;

// The base class for all map types
struct MapBase
{
	// using AttributeContainer = AttributeContainerT<Vector>;
	using AttributeContainer = AttributeContainerT<ChunkArray>;

	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	using AttributeGen = AttributeContainer::AttributeGen;
	using MarkAttribute = AttributeContainer::MarkAttribute;

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

	MapBase();
	~MapBase();

	/*************************************************************************/
	// Map-wise attributes container
	/*************************************************************************/
	std::unordered_map<std::string, std::any> attributes_;

	template <typename T>
	T& get_attribute(const std::string& name)
	{
		auto [it, inserted] = attributes_.try_emplace(name, T());
		return std::any_cast<T&>(it->second);
	}

	/*************************************************************************/
	// Add a topological relation (only used in derived maps constructors)
	/*************************************************************************/

	inline std::shared_ptr<Attribute<Dart>> add_relation(const std::string& name)
	{
		return relations_.emplace_back(darts_.add_attribute<Dart>(name));
	}

	/*************************************************************************/
	// Basic functions for iterating over darts
	/*************************************************************************/

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

/*************************************************************************/
// Darts basic functions
/*************************************************************************/

Dart add_dart(MapBase& m);
void remove_dart(MapBase& m, Dart d);

inline uint32 nb_darts(const MapBase& m)
{
	return m.darts_.nb_elements();
}

// returns the number of darts of a given cell orbit
// do not downgrade the concrete MESH type to MapBase because foreach_dart_of_orbit function is specialized for derived
// MESH types
template <typename CELL, typename MESH>
auto nb_darts_of_orbit(const MESH& m, CELL c) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>, uint32>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	uint32 result = 0;
	foreach_dart_of_orbit(m, c, [&](Dart) -> bool {
		++result;
		return true;
	});
	return result;
}

/*************************************************************************/
// Boundary marker management
/*************************************************************************/

inline void set_boundary(const MapBase& m, Dart d, bool b)
{
	(*m.boundary_marker_)[d.index] = b ? 1u : 0u;
}

inline bool is_boundary(const MapBase& m, Dart d)
{
	return (*m.boundary_marker_)[d.index] != 0u;
}

// indicates if a cell is incident to the boundary of the mesh
// do not downgrade the concrete MESH type to MapBase in case of an overload of index_of for a derived MESH type
template <typename MESH, typename CELL>
auto is_incident_to_boundary(const MESH& m, CELL c) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>, bool>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	bool result = false;
	foreach_dart_of_orbit(m, c, [&m, &result](Dart d) -> bool {
		if (is_boundary(m, d))
		{
			result = true;
			return false;
		}
		return true;
	});
	return result;
}

/*************************************************************************/
// Cells indexing management
/*************************************************************************/

// get a new index for the given CELL type
template <typename CELL>
uint32 new_index(const MapBase& m)
{
	return m.attribute_containers_[CELL::ORBIT].new_index();
}

// set the given index to the given dart (updates the ref count of the index in the container)
template <typename CELL>
void set_index(MapBase& m, Dart d, uint32 index)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");

	const uint32 old = (*m.cells_indices_[orbit])[d.index];
	// ref_index() is done before unref_index() to avoid deleting the index line if old == index
	if (index != INVALID_INDEX)
		m.attribute_containers_[orbit].ref_index(index); // ref the new index
	if (old != INVALID_INDEX)
		m.attribute_containers_[orbit].unref_index(old); // unref the old index
	(*m.cells_indices_[orbit])[d.index] = index;		 // affect the index to the dart
}

// copy the index of the given CELL type from the src dart to the dest dart
// do not downgrade the concrete MESH type to MapBase in case of an overload of index_of for a derived MESH type
template <typename CELL, typename MESH>
auto copy_index(MESH& m, Dart dest, Dart src) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	set_index<CELL>(m, dest, index_of(m, CELL(src)));
}

// get the index of the given cell
template <typename CELL>
uint32 index_of(const MapBase& m, CELL c)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	return (*m.cells_indices_[orbit])[c.dart.index];
}

// get the cell of the given index
template <typename CELL>
CELL of_index(const MapBase& m, uint32 i)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		if (!is_boundary(m, d))
		{
			const CELL c(d);
			if (index_of(m, c) == i)
				return c;
		}
	}
	return CELL();
}

// indicates if the map is currently globally indexing the given CELL type
template <typename CELL>
bool is_indexed(const MapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}

// indicates if the map is currently globally indexing the given orbit (same as previous one but not templated)
inline bool is_indexed(const MapBase& m, Orbit orbit)
{
	cgogn_message_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	return m.cells_indices_[orbit] != nullptr;
}

// creates the Dart attribute needed to store the indices of the given CELL type
template <typename CELL>
void init_cells_indexing(MapBase& m)
{
	static const Orbit orbit = CELL::ORBIT;
	static_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	if (!is_indexed<CELL>(m))
	{
		std::ostringstream oss;
		oss << "__index_" << orbit_name(orbit);
		m.cells_indices_[orbit] = m.darts_.add_attribute<uint32>(oss.str());
		m.cells_indices_[orbit]->fill(INVALID_INDEX);
	}
}

// creates the Dart attribute needed to store the indices of the given orbit (same as previous one but not templated)
inline void init_cells_indexing(MapBase& m, Orbit orbit)
{
	cgogn_message_assert(orbit < NB_ORBITS, "Unknown orbit parameter");
	if (!is_indexed(m, orbit))
	{
		std::ostringstream oss;
		oss << "__index_" << orbit_name(orbit);
		m.cells_indices_[orbit] = m.darts_.add_attribute<uint32>(oss.str());
		m.cells_indices_[orbit]->fill(INVALID_INDEX);
	}
}

// for all the darts of the given cell, sets the index of the given CELL type to the given value
// do not downgrade the concrete MESH type to MapBase because foreach_dart_of_orbit function is specialized for derived
// MESH types
template <typename CELL, typename MESH>
auto set_index(MESH& m, CELL c, uint32 index) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	cgogn_message_assert(is_indexed<CELL>(m), "Trying to access the cell index of an unindexed cell type");
	foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
		set_index<CELL>(m, d, index);
		return true;
	});
}

// gives an index to all the cells of the given CELL type
// do not downgrade the concrete MESH type to MapBase because set_index function needs it
template <typename CELL, typename MESH>
auto index_cells(MESH& m) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		init_cells_indexing<CELL>(m);

	MapBase& base = static_cast<MapBase&>(m);
	DartMarker dm(m);
	for (Dart d = base.begin(), end = base.end(); d != end; d = base.next(d))
	{
		if (!is_boundary(m, d) && !dm.is_marked(d))
		{
			const CELL c(d);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				dm.mark(d);
				return true;
			});

			if (index_of(m, c) == INVALID_INDEX)
				set_index(m, c, new_index<CELL>(m));
		}
	}
}

/*************************************************************************/
// Attributes management
/*************************************************************************/

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>* = nullptr>
std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>> add_attribute(MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		index_cells<CELL>(m);
	MapBase& mb = static_cast<MapBase&>(m);
	return mb.attribute_containers_[CELL::ORBIT].template add_attribute<T>(name);
}

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>* = nullptr>
std::shared_ptr<MapBase::Attribute<T>> get_attribute(const MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::ORBIT].template get_attribute<T>(name);
}

template <typename CELL>
void remove_attribute(MapBase& m, const std::shared_ptr<MapBase::AttributeGen>& attribute)
{
	m.attribute_containers_[CELL::ORBIT].remove_attribute(attribute);
}

template <typename CELL>
void remove_attribute(MapBase& m, MapBase::AttributeGen* attribute)
{
	m.attribute_containers_[CELL::ORBIT].remove_attribute(attribute);
}

template <typename CELL, typename FUNC>
void foreach_attribute(const MapBase& m, const FUNC& f)
{
	using AttributeGen = MapBase::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeGen>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::ORBIT])
		f(a);
}

template <typename T, typename CELL, typename FUNC>
void foreach_attribute(const MapBase& m, const FUNC& f)
{
	using AttributeT = MapBase::Attribute<T>;
	using AttributeGen = MapBase::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeT>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::ORBIT])
	{
		std::shared_ptr<AttributeT> at = std::dynamic_pointer_cast<AttributeT>(a);
		if (at)
			f(at);
	}
}

template <typename T>
T& get_attribute(MapBase& m, const std::string& name)
{
	return m.get_attribute<T>(name);
}

template <typename CELL, typename MESH>
auto get_mark_attribute(const MESH& m)
	-> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>, typename mesh_traits<MESH>::MarkAttribute*>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		index_cells<CELL>(const_cast<MESH&>(m));
	const MapBase& mb = static_cast<const MapBase&>(m);
	return mb.attribute_containers_[CELL::ORBIT].get_mark_attribute();
}

template <typename CELL>
void release_mark_attribute(const MapBase& m, MapBase::MarkAttribute* attribute)
{
	return m.attribute_containers_[CELL::ORBIT].release_mark_attribute(attribute);
}

inline typename MapBase::MarkAttribute* get_dart_mark_attribute(const MapBase& m)
{
	return m.darts_.get_mark_attribute();
}

inline void release_dart_mark_attribute(const MapBase& m, MapBase::MarkAttribute* attribute)
{
	return m.darts_.release_mark_attribute(attribute);
}

/*************************************************************************/
// Global cells traversals
/*************************************************************************/

template <typename MESH, typename FUNC>
auto foreach_cell(const MESH& m, const FUNC& func) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	foreach_cell(m, func, MapBase::TraversalPolicy::AUTO);
}

// call the given function on each cell of the given MESH
// the type of the cell is deduced from the function parameter type
// do not downgrade the concrete MESH type to MapBase because several function are specialized for derived MESH types
template <typename MESH, typename FUNC>
auto foreach_cell(const MESH& m, const FUNC& f, MapBase::TraversalPolicy traversal_policy)
	-> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	using CELL = func_parameter_type<FUNC>;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, CELL>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (traversal_policy == MapBase::TraversalPolicy::AUTO && is_indexed<CELL>(m))
	{
		CellMarker<MESH, CELL> cm(m);
		for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		{
			const CELL c(d);
			if (!is_boundary(m, d) && !cm.is_marked(c))
			{
				cm.mark(c);
				if (!f(c))
					break;
			}
		}
	}
	else
	{
		DartMarker dm(m);
		for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		{
			if (!is_boundary(m, d) && !dm.is_marked(d))
			{
				const CELL c(d);
				foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
					dm.mark(d);
					return true;
				});
				if (!f(c))
					break;
			}
		}
	}
}

template <typename MESH, typename FUNC>
auto parallel_foreach_cell(const MESH& m, const FUNC& f) -> std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>
{
	using CELL = func_parameter_type<FUNC>;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, CELL>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	ThreadPool* pool = thread_pool();
	uint32 nb_workers = pool->nb_workers();
	if (nb_workers == 0)
		return foreach_cell(m, f);

	using VecCell = std::vector<uint32>;
	using Future = std::future<void>;

	std::array<std::vector<VecCell*>, 2> cells_buffers;
	std::array<std::vector<Future>, 2> futures;
	cells_buffers[0].reserve(nb_workers);
	cells_buffers[1].reserve(nb_workers);
	futures[0].reserve(nb_workers);
	futures[1].reserve(nb_workers);

	Buffers<uint32>* buffers = uint32_buffers();

	Dart it = m.begin();
	Dart last = m.end();

	uint32 i = 0u; // buffer id (0/1)
	uint32 j = 0u; // thread id (0..nb_workers)

	if (is_indexed<CELL>(m))
	{
		CellMarker<MESH, CELL> cm(m);
		while (it.index < last.index)
		{
			// fill buffer
			cells_buffers[i].push_back(buffers->buffer());
			VecCell& cells = *cells_buffers[i].back();
			cells.reserve(PARALLEL_BUFFER_SIZE);
			for (uint32 k = 0u; k < PARALLEL_BUFFER_SIZE && it.index < last.index;)
			{
				CELL c(it);
				if (!is_boundary(m, it) && !cm.is_marked(c))
				{
					cm.mark(c);
					cells.push_back(c.dart.index);
					++k;
				}
				it = m.next(it);
			}
			// launch thread
			futures[i].push_back(pool->enqueue([&cells, &f]() {
				for (uint32 index : cells)
					f(CELL(Dart(index)));
			}));
			// next thread
			if (++j == nb_workers)
			{ // again from 0 & change buffer
				j = 0u;
				i = (i + 1u) % 2u;
				for (auto& fu : futures[i])
					fu.wait();
				for (auto& b : cells_buffers[i])
					buffers->release_buffer(b);
				futures[i].clear();
				cells_buffers[i].clear();
			}
		}
	}
	else
	{
		DartMarker dm(m);
		while (it.index < last.index)
		{
			// fill buffer
			cells_buffers[i].push_back(buffers->buffer());
			VecCell& cells = *cells_buffers[i].back();
			cells.reserve(PARALLEL_BUFFER_SIZE);
			for (uint32 k = 0u; k < PARALLEL_BUFFER_SIZE && it.index < last.index;)
			{
				if (!is_boundary(m, it) && !dm.is_marked(it))
				{
					CELL c(it);
					foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
						dm.mark(d);
						return true;
					});
					cells.push_back(c.dart.index);
					++k;
				}
				it = m.next(it);
			}
			// launch thread
			futures[i].push_back(pool->enqueue([&cells, &f]() {
				for (uint32 index : cells)
					f(CELL(Dart(index)));
			}));
			// next thread
			if (++j == nb_workers)
			{ // again from 0 & change buffer
				j = 0u;
				i = (i + 1u) % 2u;
				for (auto& fu : futures[i])
					fu.wait();
				for (auto& b : cells_buffers[i])
					buffers->release_buffer(b);
				futures[i].clear();
				cells_buffers[i].clear();
			}
		}
	}

	// clean all at the end
	for (auto& fu : futures[0u])
		fu.wait();
	for (auto& b : cells_buffers[0u])
		buffers->release_buffer(b);
	for (auto& fu : futures[1u])
		fu.wait();
	for (auto& b : cells_buffers[1u])
		buffers->release_buffer(b);
}

/*************************************************************************/
// Clear map
/*************************************************************************/

void clear(MapBase& m, bool keep_attributes = true);

/*************************************************************************/
// Copy map
/*************************************************************************/

template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>* = nullptr>
void copy(MESH& dst, const MESH& src)
{
	clear(dst, false);
	for (uint32 orbit = 0; orbit < NB_ORBITS; ++orbit)
	{
		if (src.cells_indices_[orbit] != nullptr)
			init_cells_indexing(dst, Orbit(orbit));
	}
	dst.darts_.copy(src.darts_);
	for (uint32 i = 0; i < NB_ORBITS; ++i)
		dst.attribute_containers_[i].copy(src.attribute_containers_[i]);
	dst.boundary_marker_ = dst.darts_.get_mark_attribute();
	dst.boundary_marker_->copy(*src.boundary_marker_);
}

/*************************************************************************/
// Debugging helper functions
/*************************************************************************/

void dump_map_darts(const MapBase& m);

// checks if the indexing of the map for the given CELL type is valid
// do not downgrade the concrete MESH type to MapBase because several functions are specialized for derived MESH types
template <typename CELL, typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, MapBase&>>* = nullptr>
bool check_indexing(MESH& m, bool verbose = true)
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");

	if (!is_indexed<CELL>(m))
		return true;

	bool result = true;

	auto counter = add_attribute<uint32, CELL>(m, "__cell_counter");
	counter->fill(0);

	foreach_cell(
		m,
		[&](CELL c) -> bool {
			const uint32 index = index_of(m, c);

			++(*counter)[index];

			bool valid_index = index != INVALID_INDEX;
			if (verbose && !valid_index)
				std::cerr << "Cell " << c << " (" << cell_name<CELL>(m) << ") has invalid index" << std::endl;

			bool all_darts_same_index = true;
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				const uint32 index_d = index_of(m, CELL(d));
				if (index_d != index)
				{
					if (verbose)
						std::cerr << "Cell " << c << " (" << cell_name<CELL>(m) << ") has darts with different indices"
								  << std::endl;
					all_darts_same_index = false;
				}
				return true;
			});

			result &= valid_index && all_darts_same_index;
			return true;
		},
		MapBase::TraversalPolicy::DART_MARKING);

	// check that all lines of the attribute container are used
	for (uint32 i = m.attribute_containers_[CELL::ORBIT].first_index(),
				end = m.attribute_containers_[CELL::ORBIT].last_index();
		 i != end; i = m.attribute_containers_[CELL::ORBIT].next_index(i))
	{
		if ((*counter)[i] == 0)
		{
			if (verbose)
				std::cerr << "Cell index " << i << " is not used in container " << cell_name<CELL>(m) << std::endl;
			result = false;
		}
		else
		{
			if ((*counter)[i] >= 2ul)
			{
				if (verbose)
					std::cerr << "Multiple cells with same index " << i << " in container " << cell_name<CELL>(m)
							  << std::endl;
				result = false;
			}
		}
	}

	remove_attribute<CELL>(m, counter);

	return result;
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MAPS_MAP_BASE_H_
