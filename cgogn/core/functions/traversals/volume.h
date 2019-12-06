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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>
#include <cgogn/core/functions/traversals/dart.h>
#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::vector<typename mesh_traits<MESH>::Volume> incident_volumes(const MESH& m, CELL c);

/*****************************************************************************/

///////////
// CMap2 //
///////////

std::vector<CMap2::Volume>
CGOGN_CORE_EXPORT incident_volumes(const CMap2& m, CMap2::Vertex v);

std::vector<CMap2::Volume>
CGOGN_CORE_EXPORT incident_volumes(const CMap2& m, CMap2::Edge e);

std::vector<CMap2::Volume>
CGOGN_CORE_EXPORT incident_volumes(const CMap2& m, CMap2::Face f);

///////////
// CMap3 //
///////////

std::vector<CMap3::Volume>
CGOGN_CORE_EXPORT incident_volumes(const CMap3& m, CMap3::Vertex v);

std::vector<CMap3::Volume>
CGOGN_CORE_EXPORT incident_volumes(const CMap3& m, CMap3::Face f);

std::vector<CMap3::Volume>
CGOGN_CORE_EXPORT incident_volumes(const CMap3& m, CMap3::Volume v);

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
std::vector<typename mesh_traits<MESH>::Volume>
incident_volumes(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return incident_volumes(m.mesh(), c);
}

/*****************************************************************************/

// template <typename CELL, typename MESH, typename FUNC>
// void foreach_incident_volume(const MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////
// CMap2 //
///////////

template <typename FUNC>
void foreach_incident_volume(const CMap2& m, CMap2::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	func(CMap2::Volume(v.dart));
}

template <typename FUNC>
void foreach_incident_volume(const CMap2& m, CMap2::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	func(CMap2::Volume(e.dart));
}

template <typename FUNC>
void foreach_incident_volume(const CMap2& m, CMap2::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	func(CMap2::Volume(f.dart));
}

///////////
// CMap3 //
///////////

template <typename FUNC>
void foreach_incident_volume(const CMap3& m, CMap3::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	DartMarkerStore marker(m);
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d) && !is_boundary(m,d))
		{
			foreach_dart_of_orbit(m,CMap3::Vertex2(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap3::Volume(d));
		}
		return true;
	});
}

template <typename FUNC>
void foreach_incident_volume(const CMap3& m, CMap3::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = e.dart;
	do
	{
		if (!is_boundary(m,it))
		{
			if (!func(CMap3::Volume(it)))
				break;
		}
		it = phi3(m,phi2(m,it));
	} while (it != e.dart);
}

template <typename FUNC>
void foreach_incident_volume(const CMap3& m, CMap3::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Volume>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = f.dart;
	if (!is_boundary(m,it))
		if (!func(CMap3::Volume(it)))
			return;
	it = phi3(m,it);
	if (!is_boundary(m,it))
		func(CMap3::Volume(it));
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
foreach_incident_volume(const MESH& m, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	foreach_incident_volume(*m.mesh(), c, func);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_VOLUME_H_
