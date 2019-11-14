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
	faces.reserve(8u);
	foreach_incident_face(m, v, [&] (CMap2::Face f) -> bool { faces.push_back(f); return true; });
	return faces;
}

CMap2::Face incident_face(const CMap2& m, CMap2::Edge1 e1)
{
	return CMap2::Face(e1.dart);
}

std::vector<CMap2::Face> incident_faces(const CMap2& m, CMap2::Edge e)
{
	std::vector<CMap2::Face> faces;
	faces.reserve(2u);
	foreach_incident_face(m, e, [&] (CMap2::Face f) -> bool { faces.push_back(f); return true; });
	return faces;
}

std::vector<CMap2::Face> incident_faces(const CMap2& m, CMap2::Volume v)
{
	std::vector<CMap2::Face> faces;
	faces.reserve(32u);
	foreach_incident_face(m, v, [&] (CMap2::Face f) -> bool { faces.push_back(f); return true; });
	return faces;
}

///////////
// CMap3 //
///////////

std::vector<CMap3::Face> incident_faces(const CMap3& m, CMap3::Vertex v)
{
	std::vector<CMap3::Face> faces;
	faces.reserve(16u);
	foreach_incident_face(m, v, [&] (CMap3::Face f) -> bool { faces.push_back(f); return true; });
	return faces;
}

std::vector<CMap3::Face> incident_faces(const CMap3& m, CMap3::Edge e)
{
	std::vector<CMap3::Face> faces;
	faces.reserve(16u);
	foreach_incident_face(m, e, [&] (CMap3::Face f) -> bool { faces.push_back(f); return true; });
	
	return faces;
}

std::vector<CMap3::Face> incident_faces(const CMap3& m, CMap3::Volume v)
{
	std::vector<CMap3::Face> faces;
	faces.reserve(32u);
	foreach_incident_face(m, v, [&] (CMap3::Face f) -> bool { faces.push_back(f); return true; });
	return faces;
}

} // namespace cgogn
