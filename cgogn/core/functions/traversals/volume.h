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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/tuples.h>
#include <cgogn/core/utils/type_traits.h>

#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_volume(const MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////////////////////////
// CMapBase (or convertible) //
///////////////////////////////

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_volume(const MESH& m, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Volume = typename mesh_traits<MESH>::Volume;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_convertible_v<MESH&, CMap2&> && mesh_traits<MESH>::dimension == 2)
	{
		func(Volume(c.dart));
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Edge>)
	{
		Dart d = c.dart;
		do
		{
			if (!is_boundary(m, d))
			{
				if (!func(Volume(d)))
					break;
			}
			d = phi3(m, phi2(m, d));
		} while (d != c.dart);
	}
	else if constexpr (std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3 &&
					   std::is_same_v<CELL, typename mesh_traits<MESH>::Face>)
	{
		Dart d = c.dart;
		if (!is_boundary(m, d))
			if (!func(Volume(d)))
				return;
		d = phi3(m, d);
		if (!is_boundary(m, d))
			func(Volume(d));
	}
	else
	{
		if (is_indexed<Volume>(m))
		{
			CellMarkerStore<MESH, Volume> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Volume v(d);
				if constexpr (mesh_traits<MESH>::dimension == 3) // volumes can be boundary cells
				{
					if (!marker.is_marked(v) && !is_boundary(m, d))
					{
						marker.mark(v);
						return func(v);
					}
				}
				else
				{
					if (!marker.is_marked(v))
					{
						marker.mark(v);
						return func(v);
					}
				}
				return true;
			});
		}
		else
		{
			DartMarkerStore<MESH> marker(m);
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				Volume v(d);
				if constexpr (mesh_traits<MESH>::dimension == 3) // volumes can be boundary cells
				{
					if (!marker.is_marked(d) && !is_boundary(m, d))
					{
						foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(v);
					}
				}
				else
				{
					if (!marker.is_marked(d))
					{
						foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
							marker.mark(d);
							return true;
						});
						return func(v);
					}
				}
				return true;
			});
		}
	}
}

////////////////////////////
// CMap3 (or convertible) //
////////////////////////////

template <typename MESH, typename FUNC>
auto foreach_incident_volume(const MESH& m, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3>
{
	using Volume = typename mesh_traits<MESH>::Volume;
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	if (is_indexed<Volume>(m))
	{
		CellMarkerStore<MESH, Volume> marker(m);
		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
			Volume vol(d);
			if (!marker.is_marked(vol) && !is_boundary(m, d))
			{
				marker.mark(vol);
				return func(vol);
			}
			return true;
		});
	}
	else
	{
		DartMarkerStore<MESH> marker(m);
		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
			Volume vol(d);
			if (!marker.is_marked(d) && !is_boundary(m, d))
			{
				foreach_dart_of_orbit(m, typename mesh_traits<MESH>::Vertex2(d), [&](Dart d) -> bool {
					marker.mark(d);
					return true;
				});
				return func(vol);
			}
			return true;
		});
	}
}

template <typename MESH, typename FUNC>
auto foreach_incident_volume(const CMap3& m, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMap3&> && mesh_traits<MESH>::dimension == 3>
{
	using Volume = typename mesh_traits<MESH>::Volume;
	static_assert(is_func_parameter_same<FUNC, Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = f.dart;
	if (!is_boundary(m, it))
		if (!func(Volume(it)))
			return;
	it = phi3(m, it);
	if (!is_boundary(m, it))
		func(Volume(it));
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Volume> incident_volumes(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::Volume> incident_volumes(const MESH& m, CELL c)
{
	using Volume = typename mesh_traits<MESH>::Volume;
	if constexpr (mesh_traits<MESH>::dimension == 2)
		return {Volume(c.dart)};
	else
	{
		std::vector<Volume> volumes;
		volumes.reserve(32u);
		foreach_incident_volume(m, c, [&](Volume v) -> bool {
			volumes.push_back(v);
			return true;
		});
		return volumes;
	}
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_
