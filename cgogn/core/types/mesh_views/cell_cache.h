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

#ifndef CGOGN_CORE_TYPES_MESH_VIEWS_CELL_CACHE_H_
#define CGOGN_CORE_TYPES_MESH_VIEWS_CELL_CACHE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/functions/traversals/global.h>

namespace cgogn
{

template <class> struct VectorsFromTuple;
template <template <typename ...Args> class tuple, typename ...T>
struct VectorsFromTuple<tuple<T...>>
{
	using type = std::tuple<std::vector<T>...>;
};

template <typename MESH>
class CellCache
{
	using CellVectors = typename VectorsFromTuple<typename mesh_traits<MESH>::Cells>::type;

	MESH& m_;
	CellVectors cells_;

	template <typename CELL>
	const std::vector<CELL>& cell_vector() const
	{
		return std::get<tuple_type_index<std::vector<CELL>, CellVectors>::value>(cells_);
	}

	template <typename CELL>
	std::vector<CELL>& cell_vector()
	{
		return std::get<tuple_type_index<std::vector<CELL>, CellVectors>::value>(cells_);
	}

public:

	static const bool is_mesh_view = true;
	using MeshType = MESH;

	template <typename CELL>
	typename std::vector<CELL>::const_iterator begin() const
	{
		return cell_vector<CELL>().begin();
	}

	template <typename CELL>
	typename std::vector<CELL>::const_iterator end() const
	{
		return cell_vector<CELL>().end();
	}

	CellCache(MESH& m) : m_(m) {}

	MESH& mesh() { return m_; }
	const MESH& mesh() const { return m_; }

	template <typename CELL>
	void build()
	{
		std::vector<CELL>& cells = cell_vector<CELL>();
		cells.clear();
		foreach_cell(m_, [&] (CELL c) -> bool { cells.push_back(c); return true; });
	}
};

template <typename MESH>
struct mesh_traits<CellCache<MESH>> : public mesh_traits<MESH>
{};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MESH_VIEWS_CELL_CACHE_H_
