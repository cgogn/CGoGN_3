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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_GLOBAL_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_GLOBAL_H_

#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/cell_marker.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename FUNC>
// void foreach_cell(MESH& m, const FUNC& f);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

template <typename MESH, typename FUNC,
		  typename = typename std::enable_if<std::is_base_of<CMapBase, MESH>::value>::type>
void
foreach_cell(const MESH& m, const FUNC& f, bool force_dart_marking = false)
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (!force_dart_marking && m.template is_embedded<CELL>())
	{
		CellMarker<CELL> cm(m);
		m.foreach_dart([&] (Dart d) -> bool
		{
			const CELL c(d);
			if (!m.is_boundary(d) && !cm.is_marked(c))
			{
				cm.mark(c);
				return f(c);
			}
			return true;
		});
	}
	else
	{
		DartMarker dm(m);
		m.foreach_dart([&] (Dart d) -> bool
		{
			if (!m.is_boundary(d) && !dm.is_marked(d))
			{
				const CELL c(d);
				m.foreach_dart_of_orbit(c, [&] (Dart d) -> bool { dm.mark(d); return true; });
				return f(c);
			}
			return true;
		});
	}
}

///////////////
// CellCache //
///////////////

template <typename MESH>
class CellCache;

template <typename MESH, typename FUNC>
void
foreach_cell(const CellCache<MESH>& cc, const FUNC& f)
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	for (auto it = cc.template begin<CELL>(), end = cc.template end<CELL>(); it != end; it++)
		if (!f(*it))
			break;
}

////////////////
// CellFilter //
////////////////

template <typename MESH>
class CellFilter;

template <typename MESH, typename FUNC>
void
foreach_cell(const CellFilter<MESH>& cf, const FUNC& f)
{
	using CELL = func_parameter_type<FUNC>;
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	foreach_cell(cf.mesh(), [&] (CELL c) -> bool
	{
		if (cf.filter(c))
			return f(c);
		return true;
	});
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_GLOBAL_H_
