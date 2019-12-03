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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_INFO_H_
#define CGOGN_CORE_FUNCTIONS_MESH_INFO_H_

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/core/functions/attributes.h>

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/volume.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::string cell_name(const MESH& m);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename MESH>
std::string cell_name(const MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	return mesh_traits<MESH>::cell_names[tuple_type_index<CELL, typename mesh_traits<MESH>::Cells>::value];
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 nb_cells(const MESH& m);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename MESH>
uint32 nb_cells(const MESH& m)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	uint32 result = 0;
	foreach_cell(m, [&] (CELL) -> bool { ++result; return true; });
	return result;
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// uint32 degree(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH>
uint32 degree(const MESH& m, typename mesh_traits<MESH>::Vertex v)
{
	uint32 result = 0;
	foreach_incident_edge(m, v, [&] (typename mesh_traits<MESH>::Edge) -> bool { ++result; return true; });
	return result;
}

template <typename MESH>
uint32 degree(const MESH& m, typename mesh_traits<MESH>::Edge e)
{
	uint32 result = 0;
	foreach_incident_face(m, e, [&] (typename mesh_traits<MESH>::Face) -> bool { ++result; return true; });
	return result;
}

template <typename MESH>
uint32 degree(const MESH& m, typename mesh_traits<MESH>::Face f)
{
	uint32 result = 0;
	foreach_incident_volume(m, f, [&] (typename mesh_traits<MESH>::Volume) -> bool { ++result; return true; });
	return result;
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// uint32 codegree(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH>
uint32 codegree(const MESH& m, typename mesh_traits<MESH>::Edge e)
{
	uint32 result = 0;
	foreach_incident_vertex(m, e, [&] (typename mesh_traits<MESH>::Vertex) -> bool { ++result; return true; });
	return result;
}

template <typename MESH>
uint32 codegree(const MESH& m, typename mesh_traits<MESH>::Face f)
{
	uint32 result = 0;
	foreach_incident_edge(m, f, [&] (typename mesh_traits<MESH>::Edge) -> bool { ++result; return true; });
	return result;
}

template <typename MESH>
uint32 codegree(const MESH& m, typename mesh_traits<MESH>::Volume v)
{
	uint32 result = 0;
	foreach_incident_face(m, v, [&] (typename mesh_traits<MESH>::Face) -> bool { ++result; return true; });
	return result;
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// bool is_incident_to_boundary(const MESH& m, CELL c);

/*****************************************************************************/

template <typename MESH, typename CELL,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
bool
is_incident_to_boundary(const MESH& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	bool result = false;
	foreach_dart_of_orbit(m,c, [&m, &result] (Dart d) -> bool
	{
        if (is_boundary(m,d)) { result = true; return false; }
		return true;
	});
	return result;
}

/*****************************************************************************/

// template <typename MESH>
// bool edge_can_collapse(const MESH& m, typename mesh_traits<MESH>::Edge e);

/*****************************************************************************/

///////////
// CMap2 //
///////////

inline bool edge_can_collapse(const CMap2& m, CMap2::Edge e)
{
	using Vertex = CMap2::Vertex;
	using Face = CMap2::Face;

	auto vertices = incident_vertices(m, e);

	if (is_incident_to_boundary(m, vertices[0]) || is_incident_to_boundary(m, vertices[1]))
		return false;

	uint32 val_v1 = degree(m, vertices[0]);
	uint32 val_v2 = degree(m, vertices[1]);

	if (val_v1 + val_v2 < 8 || val_v1 + val_v2 > 14)
		return false;

	Dart e1 = e.dart;
	Dart e2 = phi2(m,e.dart);

	if (codegree(m, Face(e1)) == 3)
	{
		if (degree(m, Vertex(phi_1(m,e1))) < 4)
			return false;
	}

	if (codegree(m, Face(e2)) == 3)
	{
		if (degree(m, Vertex(phi_1(m,e2))) < 4)
			return false;
	}

	auto next_edge = [&m] (Dart d) { return phi2(m,phi_1(m,d)); };

	// Check vertex sharing condition
	std::vector<uint32> vn1;
	Dart it = next_edge(next_edge(e1));
	Dart end = phi1(m,e2);
	do
	{
		vn1.push_back(index_of(m, Vertex(phi1(m,it))));
		it = next_edge(it);
	} while (it != end);
	it = next_edge(next_edge(e2));
	end = phi1(m,e1);
	do
	{
		auto vn1it = std::find(vn1.begin(), vn1.end(), index_of(m, Vertex(phi1(m,it))));
		if (vn1it != vn1.end())
			return false;
		it = next_edge(it);
	} while (it != end);

	return true;
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
bool
edge_can_collapse(const MESH& m, typename mesh_traits<MESH>::Edge e)
{
	return edge_can_collapse(m.mesh(), e);
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_INFO_H_
