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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>
#include <cgogn/core/functions/traversals/dart.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Face> incident_faces(MESH& m, CELL c);

/*****************************************************************************/

///////////
// CMap2 //
///////////

std::vector<CMap2::Face>
CGOGN_CORE_EXPORT incident_faces(const CMap2& m, CMap2::Vertex v);

CMap2::Face
CGOGN_CORE_EXPORT incident_face(const CMap2& m, CMap2::HalfEdge h);

std::vector<CMap2::Face>
CGOGN_CORE_EXPORT incident_faces(const CMap2& m, CMap2::Edge e);

std::vector<CMap2::Face>
CGOGN_CORE_EXPORT incident_faces(const CMap2& m, CMap2::Volume v);

///////////
// CMap3 //
///////////

std::vector<CMap3::Face>
CGOGN_CORE_EXPORT incident_faces(const CMap3& m, CMap3::Vertex v);

std::vector<CMap3::Face>
CGOGN_CORE_EXPORT incident_faces(const CMap3& m, CMap3::Edge e);

std::vector<CMap3::Face>
CGOGN_CORE_EXPORT incident_faces(const CMap3& m, CMap3::Volume v);

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
std::vector<typename mesh_traits<MESH>::Face>
incident_faces(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return incident_faces(m.mesh(), c);
}

/*****************************************************************************/

// template <typename MESH, typename CELL, typename FUNC>
// void foreach_incident_face(MESH& m, CELL c, const FUNC& f);

/*****************************************************************************/

///////////
// CMap2 //
///////////

template <typename FUNC>
void foreach_incident_face(const CMap2& m, CMap2::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
        if (!m.mesh()->is_boundary(d))
			return func(CMap2::Face(d));
		return true;
	});
}

template <typename FUNC>
void foreach_incident_face(const CMap2& m, CMap2::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_orbit(m,e, [&] (Dart d) -> bool
	{
        if (!m.mesh()->is_boundary(d))
			return func(CMap2::Face(d));
		return true;
	});
}

template <typename FUNC>
void foreach_incident_face(const CMap2& m, CMap2::Volume v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap2::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
    DartMarkerStore marker(*m.mesh());
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
        if (!marker.is_marked(d) && !m.mesh()->is_boundary(d))
		{
			foreach_dart_of_orbit(m,CMap2::Face(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap2::Face(d));
		}
		return true;
	});
}

///////////
// CMap3 //
///////////

template <typename FUNC>
void foreach_incident_face(const CMap3& m, CMap3::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
    DartMarkerStore marker(*m.mesh());
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			foreach_dart_of_orbit(m,CMap3::Face(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap3::Face(d));
		}
		return true;
	});
}

template <typename FUNC>
void foreach_incident_face(const CMap3& m, CMap3::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = e.dart;
	do
	{
		if (!func(CMap3::Face(it))) break;
		it = phi3(m,phi2(m,it));
	} while (it != e.dart);
}

template <typename FUNC>
void foreach_incident_face(const CMap3& m, CMap3::Volume v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, CMap3::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
    DartMarkerStore marker(*m.mesh());
	foreach_dart_of_orbit(m,v, [&] (Dart d) -> bool
	{
		if (!marker.is_marked(d))
		{
			// TODO: could mark only the darts of CMap2::Face(d)
			foreach_dart_of_orbit(m,CMap3::Face(d), [&] (Dart d) -> bool { marker.mark(d); return true; });
			return func(CMap3::Face(d));
		}
		return true;
	});
}

//////////////
// MESHVIEW //
//////////////

template <typename CELL, typename MESH, typename FUNC,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
foreach_incident_face(const MESH& m, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
    foreach_incident_face(*m.mesh(), c, func);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_
