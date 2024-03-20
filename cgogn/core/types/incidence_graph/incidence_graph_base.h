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

#ifndef CGOGN_CORE_TYPES_INCIDENCE_GRAPH_BASE_H_
#define CGOGN_CORE_TYPES_INCIDENCE_GRAPH_BASE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/container/vector.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/core/utils/assert.h>
#include <cgogn/core/utils/thread.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/core/utils/type_traits.h>

#include <any>
#include <array>
#include <unordered_map>

namespace cgogn
{

template <typename MESH>
struct mesh_traits;

struct IncidenceGraphBase
{
	// using AttributeContainer = AttributeContainerT<Vector>;
	using AttributeContainer = AttributeContainerT<ChunkArray>;

	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	using AttributeGen = AttributeContainer::AttributeGen;
	using MarkAttribute = AttributeContainer::MarkAttribute;

	/*************************************************************************/
	// Cells
	/*************************************************************************/

	struct Vertex
	{
		static const uint32 CELL_INDEX = 0;
		uint32 index_;
		inline Vertex() : index_(INVALID_INDEX)
		{
		}
		inline Vertex(uint32 id) : index_(id)
		{
		}
		operator uint32() const
		{
			return index_;
		}
		bool operator<(Vertex v) const
		{
			return index_ < v.index_;
		}
		bool operator==(Vertex v) const
		{
			return index_ == v.index_;
		}
		bool operator!=(Vertex v) const
		{
			return index_ != v.index_;
		}
		inline bool is_valid() const
		{
			return index_ != INVALID_INDEX;
		}
	};

	struct Edge
	{
		static const uint32 CELL_INDEX = 1;
		uint32 index_;
		inline Edge() : index_(INVALID_INDEX)
		{
		}
		inline Edge(uint32 id) : index_(id)
		{
		}
		operator uint32() const
		{
			return index_;
		}
		bool operator<(Edge e) const
		{
			return index_ < e.index_;
		}
		bool operator==(Edge e) const
		{
			return index_ == e.index_;
		}
		bool operator!=(Edge e) const
		{
			return index_ != e.index_;
		}
		inline bool is_valid() const
		{
			return index_ != INVALID_INDEX;
		}
	};

	struct Face
	{
		static const uint32 CELL_INDEX = 2;
		uint32 index_;
		inline Face() : index_(INVALID_INDEX)
		{
		}
		inline Face(uint32 id) : index_(id)
		{
		}
		operator uint32() const
		{
			return index_;
		}
		bool operator<(Face f) const
		{
			return index_ < f.index_;
		}
		bool operator==(Face f) const
		{
			return index_ == f.index_;
		}
		bool operator!=(Face f) const
		{
			return index_ != f.index_;
		}
		inline bool is_valid() const
		{
			return index_ != INVALID_INDEX;
		}
	};

	/*************************************************************************/
	// Cells attributes containers
	/*************************************************************************/
	mutable std::array<AttributeContainer, 3> attribute_containers_;

	// shortcuts to topological relations attributes
	std::shared_ptr<Attribute<std::vector<Edge>>> vertex_incident_edges_;
	std::shared_ptr<Attribute<std::pair<Vertex, Vertex>>> edge_incident_vertices_;
	std::shared_ptr<Attribute<std::vector<Face>>> edge_incident_faces_;
	std::shared_ptr<Attribute<std::vector<Edge>>> face_incident_edges_;
	std::shared_ptr<Attribute<std::vector<uint8>>> face_incident_edges_dir_;

	/*************************************************************************/
	// Graph-wise attributes container
	/*************************************************************************/
	std::unordered_map<std::string, std::any> attributes_;

	template <typename T>
	T& get_attribute(const std::string& name)
	{
		auto [it, inserted] = attributes_.try_emplace(name, T());
		return std::any_cast<T&>(it->second);
	}

	IncidenceGraphBase()
	{
		create_incidence_attributes();
	};

	~IncidenceGraphBase()
	{
	}

	void create_incidence_attributes()
	{
		vertex_incident_edges_ =
			attribute_containers_[Vertex::CELL_INDEX].add_attribute<std::vector<Edge>>("incident_edges");
		edge_incident_vertices_ =
			attribute_containers_[Edge::CELL_INDEX].add_attribute<std::pair<Vertex, Vertex>>("incident_vertices");
		edge_incident_faces_ =
			attribute_containers_[Edge::CELL_INDEX].add_attribute<std::vector<Face>>("incident_faces");
		face_incident_edges_ =
			attribute_containers_[Face::CELL_INDEX].add_attribute<std::vector<Edge>>("incident_edges");
		face_incident_edges_dir_ =
			attribute_containers_[Face::CELL_INDEX].add_attribute<std::vector<uint8>>("incident_edges_dir");
	}
};

template <typename MESH, typename CELL>
struct CellMarker;

template <typename MESH, typename CELL>
struct CellMarkerStore;

/*************************************************************************/
// Cells basic functions
/*************************************************************************/

template <typename CELL>
CELL add_cell(IncidenceGraphBase& ig)
{
	return CELL(ig.attribute_containers_[CELL::CELL_INDEX].new_index());
}

template <typename CELL>
void remove_cell(IncidenceGraphBase& ig, CELL c)
{
	ig.attribute_containers_[CELL::CELL_INDEX].release_index(c.index_);
}

/*************************************************************************/
// Cells indexing management
/*************************************************************************/

template <typename CELL>
uint32 new_index(const IncidenceGraphBase& ig)
{
	return ig.attribute_containers_[CELL::CELL_INDEX].new_index();
}

template <typename CELL>
uint32 index_of(const IncidenceGraphBase& /*ig*/, CELL c)
{
	return c.index_;
}

template <typename CELL>
CELL of_index(const IncidenceGraphBase& /*ig*/, uint32 index)
{
	return CELL(index);
}

template <typename CELL>
bool is_indexed(const IncidenceGraphBase& /*ig*/)
{
	return true;
}

/*************************************************************************/
// Attributes management
/*************************************************************************/

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>>* = nullptr>
std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>> add_attribute(MESH& m, const std::string& name)
{
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	IncidenceGraphBase& mb = static_cast<IncidenceGraphBase&>(m);
	return mb.attribute_containers_[CELL::CELL_INDEX].template add_attribute<T>(name);
}

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>>* = nullptr>
std::shared_ptr<IncidenceGraphBase::Attribute<T>> get_attribute(const MESH& m, const std::string& name)
{
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::CELL_INDEX].template get_attribute<T>(name);
}

template <typename CELL>
void remove_attribute(IncidenceGraphBase& m, const std::shared_ptr<IncidenceGraphBase::AttributeGen>& attribute)
{
	m.attribute_containers_[CELL::CELL_INDEX].remove_attribute(attribute);
}

template <typename CELL>
void remove_attribute(IncidenceGraphBase& m, IncidenceGraphBase::AttributeGen* attribute)
{
	m.attribute_containers_[CELL::CELL_INDEX].remove_attribute(attribute);
}

template <typename CELL, typename FUNC>
void foreach_attribute(const IncidenceGraphBase& m, const FUNC& f)
{
	using AttributeGen = IncidenceGraphBase::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeGen>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::CELL_INDEX])
		f(a);
}

template <typename T, typename CELL, typename FUNC>
void foreach_attribute(const IncidenceGraphBase& m, const FUNC& f)
{
	using AttributeT = IncidenceGraphBase::Attribute<T>;
	using AttributeGen = IncidenceGraphBase::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeT>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::CELL_INDEX])
	{
		std::shared_ptr<AttributeT> at = std::dynamic_pointer_cast<AttributeT>(a);
		if (at)
			f(at);
	}
}

template <typename T>
T& get_attribute(IncidenceGraphBase& ig, const std::string& name)
{
	return ig.get_attribute<T>(name);
}

template <typename CELL, typename MESH>
auto get_mark_attribute(const MESH& m)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>, typename mesh_traits<MESH>::MarkAttribute*>
{
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::CELL_INDEX].get_mark_attribute();
}

template <typename CELL, typename MESH>
auto release_mark_attribute(const MESH& m, typename mesh_traits<MESH>::MarkAttribute* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>>
{
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::CELL_INDEX].release_mark_attribute(attribute);
}

/*************************************************************************/
// Global cells traversals
/*************************************************************************/

template <typename MESH, typename FUNC>
auto foreach_cell(const MESH& ig, const FUNC& f) -> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>>
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	for (uint32 i = ig.attribute_containers_[CELL::CELL_INDEX].first_index(),
				end = ig.attribute_containers_[CELL::CELL_INDEX].last_index();
		 i != end; i = ig.attribute_containers_[CELL::CELL_INDEX].next_index(i))
	{
		CELL c(i);
		if (/*c.is_valid() && */ !f(c))
			break;
	}
}

template <typename MESH, typename FUNC>
auto parallel_foreach_cell(const MESH& m, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>>
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
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

	uint32 it = m.attribute_containers_[CELL::CELL_INDEX].first_index();
	uint32 last = m.attribute_containers_[CELL::CELL_INDEX].last_index();

	uint32 i = 0u; // buffer id (0/1)
	uint32 j = 0u; // thread id (0..nb_workers)

	while (it < last)
	{
		// fill buffer
		cells_buffers[i].push_back(buffers->buffer());
		VecCell& cells = *cells_buffers[i].back();
		cells.reserve(PARALLEL_BUFFER_SIZE);
		for (uint32 k = 0u; k < PARALLEL_BUFFER_SIZE && it < last; ++k)
		{
			cells.push_back(it);
			it = m.attribute_containers_[CELL::CELL_INDEX].next_index(it);
		}
		// launch thread
		futures[i].push_back(pool->enqueue([&cells, &f]() {
			for (uint32 index : cells)
				f(CELL(index));
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
// Clear incidence graph
/*************************************************************************/

inline void clear(IncidenceGraphBase& m, bool keep_attributes = true)
{
	if (!keep_attributes)
	{
		m.vertex_incident_edges_.reset();
		m.edge_incident_vertices_.reset();
		m.edge_incident_faces_.reset();
		m.face_incident_edges_.reset();
		m.face_incident_edges_dir_.reset();
	}

	for (IncidenceGraphBase::AttributeContainer& container : m.attribute_containers_)
	{
		container.clear_attributes();
		if (!keep_attributes)
		{
			// if there are still shared_ptr somewhere, some attributes may not be removed
			container.remove_attributes();
		}
	}

	if (!keep_attributes)
		m.create_incidence_attributes();
}

/*************************************************************************/
// Copy incidence graph
/*************************************************************************/

template <typename MESH>
auto copy(MESH& /*dst*/, const MESH& /*src*/) -> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraphBase&>>
{
	// TODO
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_INCIDENCE_GRAPH_BASE_H_
