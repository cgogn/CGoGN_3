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

#ifndef CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_
#define CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/core/utils/tuples.h>

#include <string>

namespace cgogn
{

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>> add_attribute(MESH& m, const std::string& name);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL, typename MESH,
		  typename std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>* = nullptr>
std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>> add_attribute(MESH& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	if (!is_indexed<CELL>(m))
		index_cells<CELL>(m);
	CMapBase& mb = static_cast<CMapBase&>(m);
	return mb.attribute_containers_[CELL::ORBIT].template add_attribute<T>(name);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>> get_attribute(const MESH& m, const std::string&
// name);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL>
std::shared_ptr<CMapBase::Attribute<T>> get_attribute(const CMapBase& m, const std::string& name)
{
	return m.attribute_containers_[CELL::ORBIT].template get_attribute<T>(name);
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void remove_attribute(MESH& m, std::shared_ptr<typename mesh_traits<MESH>::AttributeGen> attribute)

// template <typename CELL, typename MESH>
// void remove_attribute(MESH& m, typename mesh_traits<MESH>::AttributeGen* attribute)

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL>
void remove_attribute(CMapBase& m, std::shared_ptr<CMapBase::AttributeGen> attribute)
{
	m.attribute_containers_[CELL::ORBIT].remove_attribute(attribute);
}

template <typename CELL>
void remove_attribute(CMapBase& m, CMapBase::AttributeGen* attribute)
{
	m.attribute_containers_[CELL::ORBIT].remove_attribute(attribute);
}

/*****************************************************************************/

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_attribute(const MESH& m, const FUNC& f);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename CELL, typename FUNC>
void foreach_attribute(const CMapBase& m, const FUNC& f)
{
	using AttributeGen = CMapBase::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeGen>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::ORBIT])
		f(a);
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH, typename FUNC>
// void foreach_attribute(const MESH& m, const FUNC& f);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T, typename CELL, typename FUNC>
void foreach_attribute(const CMapBase& m, const FUNC& f)
{
	using AttributeT = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeT>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::ORBIT])
	{
		std::shared_ptr<AttributeT> at = std::dynamic_pointer_cast<AttributeT>(a);
		if (at)
			f(at);
	}
}

/*****************************************************************************/

// template <typename T, typename CELL, typename MESH>
// T& value(MESH& m, typename mesh_traits<MESH>::template AttributePtr<T> attribute, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename T, typename CELL, typename MESH>
inline T& value(const MESH& m, const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute,
				CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return (*attribute)[index_of(m, c)];
}

template <typename T, typename CELL, typename MESH>
inline T& value(const MESH& m, typename mesh_traits<MESH>::template Attribute<T>* attribute, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return (*attribute)[index_of(m, c)];
}

template <typename T, typename CELL, typename MESH>
inline const T& value(const MESH& m, const typename mesh_traits<MESH>::template Attribute<T>* attribute, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return (*attribute)[index_of(m, c)];
}

/*****************************************************************************/

// template <typename T, typename MESH>
// T& get_attribute(MESH& m, const std::string& name);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename T>
T& get_attribute(CMapBase& m, const std::string& name)
{
	return m.get_attribute<T>(name);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_ATTRIBUTES_H_
