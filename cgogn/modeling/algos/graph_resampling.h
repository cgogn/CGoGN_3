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

#ifndef CGOGN_MODELING_ALGOS_GRAPH_RESAMPLING_H_
#define CGOGN_MODELING_ALGOS_GRAPH_RESAMPLING_H_

#include <cgogn/core/types/maps/cmap/graph.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

using geometry::Vec3;
using geometry::Scalar;

///////////
// Graph //
///////////

void resample_graph(Graph& g, Graph::Attribute<Vec3>* g_vertex_position, Graph::Attribute<Scalar>* g_vertex_radius,
					Graph& new_g, Graph::Attribute<Vec3>* new_g_vertex_position,
					Graph::Attribute<Scalar>* new_g_vertex_radius, Scalar density);

void resample_branch(Graph& g, std::pair<Dart, Dart> g_branch, std::pair<Scalar, Scalar> g_branch_arclength,
					 Graph& new_g, Graph::Edge new_g_edge, Graph::Attribute<Vec3>* g_vertex_position,
					 Graph::Attribute<Scalar>* g_vertex_radius, Graph::Attribute<Vec3>* new_g_vertex_position,
					 Graph::Attribute<Scalar>* new_g_vertex_radius, Graph::Attribute<Scalar>* new_g_vertex_arclength,
					 Scalar density);

////////////////////
// IncidenceGraph //
////////////////////

void resample_graph(IncidenceGraph& g, IncidenceGraph::Attribute<Vec3>* g_vertex_position,
					IncidenceGraph::Attribute<Scalar>* g_vertex_radius, IncidenceGraph& new_g,
					IncidenceGraph::Attribute<Vec3>* new_g_vertex_position,
					IncidenceGraph::Attribute<Scalar>* new_g_vertex_radius, Scalar density);

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_GRAPH_RESAMPLING_H_
