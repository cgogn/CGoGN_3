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

#ifndef CGOGN_CORE_TYPES_MESH_VIEWS_CELL_CACHE_H_
#define CGOGN_CORE_TYPES_MESH_VIEWS_CELL_CACHE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/tuples.h>

namespace cgogn
{

template <class>
struct tuple_of_vectors_of_T_from_tuple_of_T;
template <template <typename... Args> class tuple, typename... T>
struct tuple_of_vectors_of_T_from_tuple_of_T<tuple<T...>>
{
	using type = std::tuple<std::vector<T>...>;
};

template <typename MESH>
class CellCache
{
	using CellVectors = typename tuple_of_vectors_of_T_from_tuple_of_T<typename mesh_traits<MESH>::Cells>::type;

	const MESH& m_;
	CellVectors cells_;

public:
	CellCache(const MESH& m) : m_(m)
	{
	}

	operator MESH&()
	{
		return m_;
	}
	operator const MESH&() const
	{
		return m_;
	}

	template <typename CELL>
	const std::vector<CELL>& cell_vector() const
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		return std::get<tuple_type_index<std::vector<CELL>, CellVectors>::value>(cells_);
	}

	template <typename CELL>
	std::vector<CELL>& cell_vector()
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		return std::get<tuple_type_index<std::vector<CELL>, CellVectors>::value>(cells_);
	}

	template <typename CELL>
	uint32 size() const
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		return std::get<tuple_type_index<std::vector<CELL>, CellVectors>::value>(cells_).size();
	}

	template <typename CELL>
	typename std::vector<CELL>::const_iterator begin() const
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		return cell_vector<CELL>().begin();
	}

	template <typename CELL>
	typename std::vector<CELL>::const_iterator end() const
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		return cell_vector<CELL>().end();
	}

	template <typename CELL>
	void build()
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		std::vector<CELL>& cells = cell_vector<CELL>();
		cells.clear();
		foreach_cell(m_, [&](CELL c) -> bool {
			cells.push_back(c);
			return true;
		});
	}

	template <typename CELL, typename FUNC>
	void build(const FUNC& filter)
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		static_assert(is_func_parameter_same<FUNC, CELL>::value, "Wrong function cell parameter type");
		static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
		std::vector<CELL>& cells = cell_vector<CELL>();
		cells.clear();
		foreach_cell(m_, [&](CELL c) -> bool {
			if (filter(c))
				cells.push_back(c);
			return true;
		});
	}

	template <typename CELL>
	void add(CELL c)
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		cell_vector<CELL>().push_back(c);
	}

	template <typename CELL>
	void clear()
	{
		static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
		cell_vector<CELL>().clear();
	}
};

template <typename MESH>
struct mesh_traits<CellCache<MESH>> : public mesh_traits<MESH>
{
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MESH_VIEWS_CELL_CACHE_H_
