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

#ifndef CGOGN_CORE_INCIDENCE_GRAPH_ATT_H_
#define CGOGN_CORE_INCIDENCE_GRAPH_ATT_H_

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>
#include <cgogn/core/utils/thread.h>
#include <cgogn/core/utils/thread_pool.h>

namespace cgogn
{



template <typename T, typename CELL>
std::shared_ptr<typename mesh_traits<IncidenceGraph>::template Attribute<T>> add_attribute(IncidenceGraph& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<IncidenceGraph>::Cells>::value, "CELL not supported in this MESH");
	IncidenceGraph& mb = static_cast<IncidenceGraph&>(m);
	return mb.attribute_containers_[CELL::CELL_INDEX].template add_attribute<T>(name);
}


template <typename T, typename CELL>
std::shared_ptr<IncidenceGraph::Attribute<T>> get_attribute(const IncidenceGraph& m, const std::string& name)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this MESH");
	return m.attribute_containers_[CELL::CELL_INDEX].template get_attribute<T>(name);
}

template <typename CELL>
auto get_mark_attribute(const IncidenceGraph& ig)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this MESH");
	return ig.attribute_containers_[CELL::CELL_INDEX].get_mark_attribute();
}


template <typename CELL>
void release_mark_attribute(const IncidenceGraph& ig, IncidenceGraph::MarkAttribute* attribute)
{
	return ig.attribute_containers_[CELL::CELL_INDEX].release_mark_attribute(attribute);
}

template <typename CELL>
void remove_attribute(IncidenceGraph& m, const std::shared_ptr<IncidenceGraph::AttributeGen>& attribute)
{
	m.attribute_containers_[CELL::CELL_INDEX].remove_attribute(attribute);
}

template <typename CELL>
void remove_attribute(IncidenceGraph& m, IncidenceGraph::AttributeGen* attribute)
{
	m.attribute_containers_[CELL::CELL_INDEX].remove_attribute(attribute);
}

template <typename CELL, typename FUNC>
void foreach_attribute(const IncidenceGraph& m, const FUNC& f)
{
	using AttributeGen = IncidenceGraph::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeGen>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::CELL_INDEX])
		f(a);
}


template <typename T, typename CELL, typename FUNC>
void foreach_attribute(const IncidenceGraph& m, const FUNC& f)
{
	using AttributeT = IncidenceGraph::Attribute<T>;
	using AttributeGen = IncidenceGraph::AttributeGen;
	static_assert(is_func_parameter_same<FUNC, const std::shared_ptr<AttributeT>&>::value,
				  "Wrong function attribute parameter type");
	for (const std::shared_ptr<AttributeGen>& a : m.attribute_containers_[CELL::CELL_INDEX])
	{
		std::shared_ptr<AttributeT> at = std::dynamic_pointer_cast<AttributeT>(a);
		if (at)
			f(at);
	}
}

}


#endif // CGOGN_CORE_INCIDENCE_GRAPH_OPS_H_
