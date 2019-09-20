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

#include <cgogn/core/functions/traversals/volume.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::vector<typename mesh_traits<MESH>::Vertex> incident_volumes(const MESH& m, CELL c);

/*****************************************************************************/

///////////
// CMap2 //
///////////

std::vector<CMap2::Volume> incident_volumes(const CMap2& m, CMap2::Vertex v)
{
	return { CMap2::Volume(v.dart) };
}

std::vector<CMap2::Volume> incident_volumes(const CMap2& m, CMap2::Edge e)
{
	return { CMap2::Volume(e.dart) };
}

std::vector<CMap2::Volume> incident_volumes(const CMap2& m, CMap2::Face f)
{
	return { CMap2::Volume(f.dart) };
}

///////////
// CMap3 //
///////////

std::vector<CMap3::Volume> incident_volumes(const CMap3& m, CMap3::Vertex v)
{
	std::vector<CMap3::Volume> volumes;
	volumes.reserve(32u);
	foreach_incident_volume(m, v, [&] (CMap3::Volume vol) -> bool { volumes.push_back(vol); return true; });
	return volumes;
}

std::vector<CMap3::Volume> incident_volumes(const CMap3& m, CMap3::Edge e)
{
	std::vector<CMap3::Volume> volumes;
	volumes.reserve(16u);
	foreach_incident_volume(m, e, [&] (CMap3::Volume v) -> bool { volumes.push_back(v); return true; });
	return volumes;
}

std::vector<CMap3::Volume> incident_volumes(const CMap3& m, CMap3::Face f)
{
	std::vector<CMap3::Volume> volumes;
	volumes.reserve(2u);
	foreach_incident_volume(m, f, [&] (CMap3::Volume v) -> bool { volumes.push_back(v); return true; });
	return volumes;
}

} // namespace cgogn
