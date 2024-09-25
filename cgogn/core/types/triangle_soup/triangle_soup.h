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

#ifndef CGOGN_CORE_TYPES_TRIANGLE_SOUP_H_
#define CGOGN_CORE_TYPES_TRIANGLE_SOUP_H_

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/utils/assert.h>
#include <cgogn/core/utils/thread.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/utils/numerics.h>

#include <array>

namespace cgogn
{

struct TriangleSoup
{
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

	struct Face
	{
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

	using Cells = std::tuple<Vertex, Face>;

	/*************************************************************************/
	// Cells attributes containers
	/*************************************************************************/
	// 0 is for vertices
	// 1 is for faces
	mutable std::array<AttributeContainer, 2> attribute_containers_;

	// shortcuts to topological relations attributes
	std::shared_ptr<Attribute<std::array<Vertex, 3>>> face_incident_vertices_;

	template <typename CELL>
	static uint32 constexpr cell_container_index()
	{
		static_assert(is_in_tuple<CELL, Cells>::value, "CELL not supported in this MESH");
		if constexpr (std::is_same_v<CELL, Vertex>)
			return 0;
		else if constexpr (std::is_same_v<CELL, Face>)
			return 1;
	}

	TriangleSoup()
	{
		face_incident_vertices_ = attribute_containers_[1].add_attribute<std::array<Vertex, 3>>("incident_vertices");
	}

	~TriangleSoup()
	{
	}
};

template <>
struct mesh_traits<TriangleSoup>
{
	static constexpr const char* name = "TriangleSoup";
	static constexpr const uint8 dimension = 2;

	using Vertex = TriangleSoup::Vertex;
	using Face = TriangleSoup::Face;

	using Cells = std::tuple<Vertex, Face>;
	static constexpr const char* cell_names[] = {"Vertex", "Face"};

	template <typename T>
	using Attribute = TriangleSoup::Attribute<T>;
	using AttributeGen = TriangleSoup::AttributeGen;
	using MarkAttribute = TriangleSoup::MarkAttribute;
};

/*************************************************************************/
// Cells indexing management
/*************************************************************************/

template <typename CELL>
uint32 new_index(const TriangleSoup& ts)
{
	static_assert(is_in_tuple<CELL, TriangleSoup::Cells>::value, "CELL not supported in this MESH");
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	return ts.attribute_containers_[container_index].new_index();
}

template <typename CELL>
uint32 index_of(const TriangleSoup& /*ts*/, CELL c)
{
	static_assert(is_in_tuple<CELL, TriangleSoup::Cells>::value, "CELL not supported in this MESH");
	return c;
}

template <typename CELL>
bool is_indexed(const TriangleSoup& /*ts*/)
{
	static_assert(is_in_tuple<CELL, TriangleSoup::Cells>::value, "CELL not supported in this MESH");
	return true;
}

/*************************************************************************/
// Attributes management
/*************************************************************************/

template <typename T, typename CELL>
std::shared_ptr<TriangleSoup::Attribute<T>> add_attribute(TriangleSoup& ts, const std::string& name)
{
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	return ts.attribute_containers_[container_index].template add_attribute<T>(name);
}

template <typename T, typename CELL>
std::shared_ptr<TriangleSoup::Attribute<T>> get_attribute(const TriangleSoup& ts, const std::string& name)
{
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	return ts.attribute_containers_[container_index].template get_attribute<T>(name);
}

template <typename CELL>
void remove_attribute(TriangleSoup& ts, const std::shared_ptr<TriangleSoup::AttributeGen>& attribute)
{
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	ts.attribute_containers_[container_index].remove_attribute(attribute);
}

template <typename CELL>
void remove_attribute(TriangleSoup& ts, TriangleSoup::AttributeGen* attribute)
{
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	ts.attribute_containers_[container_index].remove_attribute(attribute);
}

template <typename CELL, typename FUNC>
void foreach_attribute(const TriangleSoup& m, const FUNC& f)
{
	using AttributeGen = TriangleSoup::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeGen>&>::value,
				  "Wrong function attribute parameter type");
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[container_index])
		f(a);
}

template <typename T, typename CELL, typename FUNC>
void foreach_attribute(const TriangleSoup& m, const FUNC& f)
{
	using AttributeT = TriangleSoup::Attribute<T>;
	using AttributeGen = TriangleSoup::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeT>&>::value,
				  "Wrong function attribute parameter type");
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[container_index])
	{
		std::shared_ptr<AttributeT> at = std::dynamic_pointer_cast<AttributeT>(a);
		if (at)
			f(at);
	}
}

template <typename CELL>
TriangleSoup::MarkAttribute* get_mark_attribute(const TriangleSoup& m)
{
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	return m.attribute_containers_[container_index].get_mark_attribute();
}

template <typename CELL>
void release_mark_attribute(const TriangleSoup& m, TriangleSoup::MarkAttribute* attribute)
{
	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	return m.attribute_containers_[container_index].release_mark_attribute(attribute);
}

/*************************************************************************/
// Global cells traversals
/*************************************************************************/

template <typename MESH, typename FUNC>
auto foreach_cell(const MESH& m, const FUNC& f) -> std::enable_if_t<std::is_convertible_v<MESH&, TriangleSoup&>>
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();
	for (uint32 i = m.attribute_containers_[container_index].first_index(),
				end = m.attribute_containers_[container_index].last_index();
		 i != end; i = m.attribute_containers_[container_index].next_index(i))
	{
		CELL c(i);
		if (!f(c))
			break;
	}
}

template <typename MESH, typename FUNC>
auto parallel_foreach_cell(const MESH& m, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, TriangleSoup&>>
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

	constexpr uint32 container_index = TriangleSoup::cell_container_index<CELL>();

	uint32 it = m.attribute_containers_[container_index].first_index();
	uint32 last = m.attribute_containers_[container_index].last_index();

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
			it = m.attribute_containers_[container_index].next_index(it);
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
// Clear mesh
/*************************************************************************/

void clear(TriangleSoup& ts, bool keep_attributes = true);

/*************************************************************************/
// Copy mesh
/*************************************************************************/

void copy(TriangleSoup& dst, const TriangleSoup& src);

/*************************************************************************/
// Specific mesh info overload
/*************************************************************************/

inline uint32 codegree(const TriangleSoup& ts, TriangleSoup::Face f)
{
	return 3;
}

/*************************************************************************/
// Local vertex traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& ts, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, TriangleSoup&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, Face>)
	{
		const std::array<Vertex, 3>& vertices = (*ts.face_incident_vertices_)[c];
		for (Vertex v : vertices)
			if (!func(v))
				break;
	}
}

/*************************************************************************/
// Operators
/*************************************************************************/

TriangleSoup::Face add_face(TriangleSoup& ts, TriangleSoup::Vertex v1, TriangleSoup::Vertex v2,
							TriangleSoup::Vertex v3);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_TRIANGLE_SOUP_H_
