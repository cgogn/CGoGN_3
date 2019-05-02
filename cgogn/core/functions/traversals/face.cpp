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

#include <cgogn/core/functions/traversals/face.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH, typename CELL>
// std::vector<typename mesh_traits<MESH>::Face> incident_faces(MESH& m, CELL c);

/*****************************************************************************/

///////////
// CMap2 //
///////////

std::vector<CMap2::Face> incident_faces(const CMap2& m, CMap2::Vertex v)
{
	std::vector<CMap2::Face> faces;
	m.foreach_dart_of_orbit(v, [&] (Dart d) -> bool
	{
		if (!m.is_boundary(d))
			faces.push_back(CMap2::Face(d));
		return true;
	});
	return faces;
}

std::vector<CMap2::Face> incident_faces(const CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Face> faces;
	m.foreach_dart_of_orbit(e, [&] (Dart d) -> bool
	{
		if (!m.is_boundary(d))
			faces.push_back(CMap2::Face(d));
		return true;
	});
	return faces;
}

} // namespace cgogn
