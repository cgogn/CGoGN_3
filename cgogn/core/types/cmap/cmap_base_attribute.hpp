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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_ATTRIBUTE_HPP_
#define CGOGN_CORE_CMAP_CMAP_BASE_ATTRIBUTE_HPP_

namespace cgogn
{


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


} // namespace cgogn

#endif // CGOGN_CORE_CMAP_CMAP_BASE_ATTRIBUTE_HPP_
