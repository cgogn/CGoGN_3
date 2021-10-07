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

#ifndef CGOGN_CORE_TYPES_MESH_TRAITS_H_
#define CGOGN_CORE_TYPES_MESH_TRAITS_H_

#include <cgogn/core/cgogn_core_export.h>
#include <cgogn/core/utils/numerics.h>

namespace cgogn
{

template <typename MESH>
struct mesh_traits
{
	// static constexpr const char* name;
	// static constexpr const uint8 dimension;

	// Cell Types
	// using Vertex = Graph::Vertex;
	// using HalfEdge = Graph::HalfEdge;
	// using Edge = Graph::Edge;
	// ...

	// Cells tuple definition
	// using Cells = std::tuple<Vertex, HalfEdge, Edge>;

	// Cells names
	// static constexpr const char* cell_names[] = {"Vertex", "HalfEdge", "Edge"};

	// Attribute types
	// template <typename T>
	// using Attribute;
	// using AttributeGen;
	// using MarkAttribute;
};

template <typename MESH>
constexpr uint8 dimension_of(const MESH&)
{
	return mesh_traits<MESH>::dimension;
}

template <typename MESH>
constexpr bool is_dimension_of(const MESH&, uint8 dim)
{
	return mesh_traits<MESH>::dimension == dim;
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_MESH_TRAITS_H_
