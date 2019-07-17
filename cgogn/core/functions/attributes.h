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

#include <cgogn/core/functions/traversals/global.h>

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
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::template AttributePtr<T>
add_attribute(MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!m.template is_embedded<CELL>())
	{
		m.template init_embedding<CELL>();
		foreach_cell(m, [&] (CELL c) -> bool
		{
			create_embedding(m, c);
			return true;
		}, true);
	}
	return m.attribute_containers_[CELL::ORBIT].template add_attribute<T>(name);
}

//////////////
// MESHVIEW //
//////////////

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::template AttributePtr<T>
add_attribute(MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return add_attribute<T, CELL>(m.mesh(), name);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// typename mesh_traits<MESH>::template AttributePtr<T> get_attribute(const MESH& m, const std::string& name);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::template AttributePtr<T>
get_attribute(const MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::ORBIT].template get_attribute<T>(name);
}

//////////////
// MESHVIEW //
//////////////

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
typename mesh_traits<MESH>::template AttributePtr<T>
get_attribute(const MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return get_attribute<T, CELL>(m.mesh(), name);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// void remove_attribute(MESH& m, typename mesh_traits<MESH>::template AttributePtr<T> attribute)

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
void
remove_attribute(MESH& m, typename mesh_traits<MESH>::AttributeGenPtr attribute)
{
	static_assert (is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL note supported in this MESH");
	m.attribute_containers_[CELL::ORBIT].remove_attribute(attribute);
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
remove_attribute(MESH& m, typename mesh_traits<MESH>::AttributeGenPtr attribute)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	remove_attribute<CELL>(m.mesh(), attribute);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH, typename FUNC>
// void foreach_attribute(const MESH& m, const FUNC& f);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type* = nullptr>
void
foreach_attribute(const MESH& m, const FUNC& f)
{
	using AttributeGenPtr = typename mesh_traits<MESH>::AttributeGenPtr;
	using AttributeT = typename mesh_traits<MESH>::template Attribute<T>;
	using AttributePtrT = typename mesh_traits<MESH>::template AttributePtr<T>;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, const AttributePtrT&>::value, "Wrong function attribute parameter type");
	for (const AttributeGenPtr& a : m.attribute_containers_[CELL::ORBIT])
	{
		AttributePtrT at = std::dynamic_pointer_cast<AttributeT>(a);
		if (at)
			f(at);
	}
}

//////////////
// MESHVIEW //
//////////////

template <typename T, typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
foreach_attribute(const MESH& m, const FUNC& f)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_attribute<T, CELL>(m.mesh(), f);
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
inline
uint32 index_of(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return m.embedding(c);
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
inline
uint32 index_of(const MESH& m, CELL c)
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
inline
T& value(const MESH& m, const typename mesh_traits<MESH>::template AttributePtr<T>& attribute, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return (*attribute)[index_of(m, c)];
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_
