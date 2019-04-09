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

#ifndef CGOGN_CORE_TYPES_MESH_VIEWS_CELL_FILTER_H_
#define CGOGN_CORE_TYPES_MESH_VIEWS_CELL_FILTER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/tuples.h>

#include <functional>

namespace cgogn
{

template <class> struct FunctionsFromTuple;
template <template <typename ...Args> class tuple, typename ...T>
struct FunctionsFromTuple<tuple<T...>>
{
	using type = std::tuple<std::function<bool(T)>...>;
};

template <typename MESH>
class CellFilter
{
	using CellFilters = typename FunctionsFromTuple<typename mesh_traits<MESH>::Cells>::type;

	MESH& m_;
	CellFilters filters_;

	template <typename CELL>
	const std::function<bool(CELL)>& cell_filter() const
	{
		return std::get<tuple_type_index<std::function<bool(CELL)>, CellFilters>::value>(filters_);
	}

	template <typename CELL>
	std::function<bool(CELL)>& cell_filter()
	{
		return std::get<tuple_type_index<std::function<bool(CELL)>, CellFilters>::value>(filters_);
	}

public:

	static const bool is_mesh_view = true;

	CellFilter(MESH& m) : m_(m) {}

	MESH& mesh() { return m_; }
	const MESH& mesh() const { return m_; }

	template <typename CELL, typename FilterFunction>
	void set_filter(const FilterFunction&& filter)
	{
		cell_filter<CELL>() = [filter] (CELL c) -> bool { return filter(c); };
	}

	template <typename CELL>
	bool filter(CELL c) const
	{
		return cell_filter<CELL>()(c);
	}
};

template <typename MESH>
struct mesh_traits<CellFilter<MESH>> : public mesh_traits<MESH>
{};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MESH_VIEWS_CELL_FILTER_H_
