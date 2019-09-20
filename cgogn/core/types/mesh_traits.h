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

#ifndef CGOGN_CORE_TYPES_MESH_TRAITS_H_
#define CGOGN_CORE_TYPES_MESH_TRAITS_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/cmap/graph.h>

namespace cgogn
{

template <typename MESH>
struct mesh_traits;

template <>
struct mesh_traits<CMap0>
{
	static constexpr const char* name = "CMap0";
	static constexpr const uint8 dimension = 0;

	using Vertex = typename CMap0::Vertex;

	using Cells = std::tuple<Vertex>;
	static constexpr const char* cell_names[] = { "Vertex" };

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

template <>
struct mesh_traits<CMap1>
{
	static constexpr const char* name = "CMap1";
	static constexpr const uint8 dimension = 1;

	using Vertex = CMap1::Vertex;
	using Edge = CMap1::Edge;
	using Face = CMap1::Face;

	using Cells = std::tuple<Vertex, Edge, Face>;
	static constexpr const char* cell_names[] = { "Vertex", "Edge", "Face" };

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

template <>
struct mesh_traits<CMap2>
{
	static constexpr const char* name = "CMap2";
	static constexpr const uint8 dimension = 2;

	using Vertex = CMap2::Vertex;
	using Edge = CMap2::Edge;
	using Face = CMap2::Face;
	using Volume = CMap2::Volume;

	using Cells = std::tuple<Vertex, Edge, Face, Volume>;
	static constexpr const char* cell_names[] = { "Vertex", "Edge", "Face", "Volume" };

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

template <>
struct mesh_traits<CMap3>
{
	static constexpr const char* name = "CMap3";
	static constexpr const uint8 dimension = 3;

	using Vertex = CMap3::Vertex;
	using Edge = CMap3::Edge;
	using Face = CMap3::Face;
	using Volume = CMap3::Volume;

	using Cells = std::tuple<Vertex, Edge, Face, Volume>;
	static constexpr const char* cell_names[] = { "Vertex", "Edge", "Face", "Volume" };

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

template <>
struct mesh_traits<Graph>
{
	static constexpr const char* name = "Graph";
	static constexpr const uint8 dimension = 1;

	using Vertex = Graph::Vertex;
	using Edge = Graph::Edge;

	using Cells = std::tuple<Vertex, Edge>;
	static constexpr const char* cell_names[] = { "Vertex", "Edge" };

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MESH_TRAITS_H_
