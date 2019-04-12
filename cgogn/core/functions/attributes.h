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

#ifndef CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_
#define CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/cmap/cmap_ops.h>

#include <cgogn/core/utils/tuples.h>

#include <string>

namespace cgogn
{

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// typename mesh_traits<MESH>::template AttributePtr<T> add_attribute(MESH& m, const std::string& name);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL, typename MESH,
		  typename = typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type>
typename mesh_traits<MESH>::template AttributePtr<T>
add_attribute(MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!m.template is_embedded<CELL>())
	{
		m.template create_embedding<CELL>();
		create_embeddings<CELL>(m);
	}
	return m.attribute_containers_[CELL::ORBIT].template add_chunk_array<T>(name);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// typename mesh_traits<MESH>::template AttributePtr<T> get_attribute(const MESH& m, const std::string& name);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL, typename MESH,
		  typename = typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type>
typename mesh_traits<MESH>::template AttributePtr<T>
get_attribute(const MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::ORBIT].template get_chunk_array<T>(name);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// void remove_attribute(MESH& m, typename mesh_traits<MESH>::template AttributePtr<T> attribute)

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH,
		  typename = typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type>
void remove_attribute(MESH& m, typename mesh_traits<MESH>::AttributeGenPtr attribute)
{
	static_assert (is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL note supported in this MESH");
	m.attribute_containers_[CELL::ORBIT].remove_chunk_array(attribute);
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 index_of(MESH& m, CELL c);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
uint32
index_of(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return m.embedding(c);
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
uint32
index_of(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return index_of(m.mesh(), c);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// T& value(MESH& m, typename mesh_traits<MESH>::template AttributePtr<T> attribute, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename T, typename CELL, typename MESH>
const T&
value(const MESH& m, typename mesh_traits<MESH>::template AttributePtr<const T> attribute, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return (*attribute)[index_of(m, c)];
}

template <typename T, typename CELL, typename MESH>
T&
value(const MESH& m, typename mesh_traits<MESH>::template AttributePtr<T> attribute, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return (*attribute)[index_of(m, c)];
}

/*****************************************************************************/

// template <typename T, typename FUNC>
// void foreach_value(typename mesh_traits<MESH>::template AttributePtr<T> attribute, const FUNC& f);

/*****************************************************************************/

////////////////
// ChunkArray //
////////////////

template <typename T, typename FUNC>
void
foreach_value(CMapBase::AttributePtr<T> attribute, const FUNC& f)
{
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	static_assert(is_ith_func_parameter_same<FUNC, 0, T&>::value, "Wrong function parameter type");
	static_assert(is_ith_func_parameter_same<FUNC, 1, uint32>::value, "Wrong function parameter type");

	for (typename CMapBase::Attribute<T>::iterator it = attribute->begin(), end = attribute->end(); it != end; ++it)
	{
		if (!f(*it, it.index()))
			break;
	}
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_
