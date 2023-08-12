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

#include <cgogn/modeling/algos/graph_resampling.h>

#include <cgogn/core/functions/traversals/halfedge.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/distance.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/picking.h>

#include <cgogn/modeling/algos/graph_utils.h>

namespace cgogn
{

namespace modeling
{

namespace internal
{

Scalar branch_length(const Graph& g, std::pair<Dart, Dart> g_branch, const Graph::Attribute<Vec3>* g_vertex_position)
{
	Scalar branch_length = 0;
	Dart cur = g_branch.first;
	Dart next = alpha0(g, cur);
	do
	{
		branch_length += (value<Vec3>(g, g_vertex_position, Graph::Vertex(next)) -
						  value<Vec3>(g, g_vertex_position, Graph::Vertex(cur)))
							 .norm();
		cur = alpha1(g, next);
		next = alpha0(g, cur);
	} while (alpha_1(g, cur) != g_branch.second);
	return branch_length;
}

std::pair<Dart, Scalar> arclength_value(const Graph& g, std::pair<Dart, Dart> g_branch,
										const Graph::Attribute<Vec3>* g_vertex_position, Scalar branch_length, Scalar l)
{
	// find the edge in which the given value lies
	Scalar cur_arclength = 0;
	Scalar next_arclength = 0;
	Dart cur = g_branch.first;
	Dart next = alpha0(g, cur);
	do
	{
		cur_arclength = next_arclength;
		next_arclength += (value<Vec3>(g, g_vertex_position, Graph::Vertex(next)) -
						   value<Vec3>(g, g_vertex_position, Graph::Vertex(cur)))
							  .norm() /
						  branch_length;
		if (next_arclength < l)
		{
			cur = alpha1(g, next);
			next = alpha0(g, cur);
		}
	} while (next_arclength < l);
	// find where in this edge the l value lies
	Scalar alpha = (l - cur_arclength) / (next_arclength - cur_arclength);

	return {cur, 1.0 - alpha};
}

} // namespace internal

///////////
// Graph //
///////////

void resample_graph(Graph& g, Graph::Attribute<Vec3>* g_vertex_position, Graph::Attribute<Scalar>* g_vertex_radius,
					Graph& new_g, Graph::Attribute<Vec3>* new_g_vertex_position,
					Graph::Attribute<Scalar>* new_g_vertex_radius, Scalar density)
{
	GraphData g_data;
	get_graph_data(g, g_data);

	auto g_vertex_new_g_vertices = add_attribute<std::vector<Dart>, Graph::Vertex>(g, "new_g_vertices");
	auto new_g_vertex_arclength = add_attribute<Scalar, Graph::Vertex>(new_g, "new_g_vertex_arclength");

	for (auto& branch : g_data.branches)
	{
		// create an edge in the new graph corresponding to the graph branch
		Dart ed = add_dart(new_g);
		Dart edd = add_dart(new_g);
		alpha0_sew(new_g, ed, edd);

		Graph::Vertex v1(ed);
		Graph::Vertex v2(edd);

		set_index(new_g, v1, new_index<Graph::Vertex>(new_g));
		set_index(new_g, v2, new_index<Graph::Vertex>(new_g));

		// store the new edge extremity darts in the graph branch end vertices
		value<std::vector<Dart>>(g, g_vertex_new_g_vertices, Graph::Vertex(branch.first.dart_)).push_back(ed);
		value<std::vector<Dart>>(g, g_vertex_new_g_vertices, Graph::Vertex(branch.second.dart_)).push_back(edd);

		// set the position of the new graph vertices
		value<Vec3>(new_g, new_g_vertex_position, v1) =
			value<Vec3>(g, g_vertex_position, Graph::Vertex(branch.first.dart_));
		value<Vec3>(new_g, new_g_vertex_position, v2) =
			value<Vec3>(g, g_vertex_position, Graph::Vertex(branch.second.dart_));
		// set the radius of the new graph vertices
		value<Scalar>(new_g, new_g_vertex_radius, v1) =
			value<Scalar>(g, g_vertex_radius, Graph::Vertex(branch.first.dart_));
		value<Scalar>(new_g, new_g_vertex_radius, v2) =
			value<Scalar>(g, g_vertex_radius, Graph::Vertex(branch.second.dart_));
		// set the arclength of the new graph vertices
		value<Scalar>(new_g, new_g_vertex_arclength, v1) = 0.0;
		value<Scalar>(new_g, new_g_vertex_arclength, v2) = 1.0;

		resample_branch(g, {branch.first.dart_, branch.second.dart_}, {0.0, 1.0}, new_g, Graph::Edge(ed),
						g_vertex_position, g_vertex_radius, new_g_vertex_position, new_g_vertex_radius,
						new_g_vertex_arclength.get(), density);

		// average new branch vertices data along the old branch
		Scalar branch_length = internal::branch_length(g, {branch.first.dart_, branch.second.dart_}, g_vertex_position);
		std::vector<Graph::Vertex> new_vertices;
		std::vector<Scalar> filtered_arclength;
		Dart cur = ed;
		Dart next = alpha0(new_g, cur);
		if (next != edd)
		{
			do
			{
				new_vertices.push_back(Graph::Vertex(next));
				cur = alpha1(new_g, next);
				next = alpha0(new_g, cur);
			} while (next != edd);
			for (Graph::Vertex v : new_vertices)
			{
				std::vector<Graph::Vertex> av = adjacent_vertices_through_edge(new_g, v);
				filtered_arclength.push_back(0.5 * (value<Scalar>(new_g, new_g_vertex_arclength, av[0]) +
													value<Scalar>(new_g, new_g_vertex_arclength, av[1])));
			}
			for (uint32 i = 0; i < new_vertices.size(); ++i)
			{
				auto location = internal::arclength_value(g, {branch.first.dart_, branch.second.dart_},
														  g_vertex_position, branch_length, filtered_arclength[i]);
				Dart e1 = location.first;
				Dart e2 = alpha0(g, e1);
				Scalar alpha = location.second;
				Vec3 pos = alpha * value<Vec3>(g, g_vertex_position, Graph::Vertex(e1)) +
						   (1.0 - alpha) * value<Vec3>(g, g_vertex_position, Graph::Vertex(e2));
				Scalar radius = alpha * value<Scalar>(g, g_vertex_radius, Graph::Vertex(e1)) +
								(1.0 - alpha) * value<Scalar>(g, g_vertex_radius, Graph::Vertex(e2));

				value<Vec3>(new_g, new_g_vertex_position, new_vertices[i]) = pos;
				value<Scalar>(new_g, new_g_vertex_radius, new_vertices[i]) = radius;
			}
		}
	}

	for (auto& intersection : g_data.intersections)
	{
		const std::vector<Dart>& new_g_vertices =
			value<std::vector<Dart>>(g, g_vertex_new_g_vertices, Graph::Vertex(intersection));
		Dart d = new_g_vertices[0];
		for (uint32 i = 1; i < new_g_vertices.size(); ++i)
			merge_vertices(new_g, Graph::Vertex(d), Graph::Vertex(new_g_vertices[i]));
	}

	remove_attribute<Graph::Vertex>(g, g_vertex_new_g_vertices);
	remove_attribute<Graph::Vertex>(new_g, new_g_vertex_arclength);

	// remove short length-1 branches
	foreach_cell(new_g, [&](Graph::Edge e) -> bool {
		std::vector<Graph::Vertex> vertices = incident_vertices(new_g, e);
		Scalar length = geometry::length(new_g, e, new_g_vertex_position);
		uint32 deg1 = degree(new_g, vertices[0]);
		uint32 deg2 = degree(new_g, vertices[1]);
		if ((deg1 > 2 && deg2 > 2) || (deg1 > 2 && deg2 == 1) || (deg1 == 1 && deg2 > 2)) // if branch with only 1 edge
		{
			if (value<Scalar>(new_g, new_g_vertex_radius, vertices[0]) + // and radiuses sum is greater than edge length
					value<Scalar>(new_g, new_g_vertex_radius, vertices[1]) >
				length)
			{
				Vec3 mid_pos = geometry::centroid<Vec3>(new_g, e, new_g_vertex_position);
				Scalar mid_radius = geometry::centroid<Scalar>(new_g, e, new_g_vertex_radius);
				Graph::Vertex v = collapse_edge(new_g, e);
				value<Vec3>(new_g, new_g_vertex_position, v) = mid_pos;
				value<Scalar>(new_g, new_g_vertex_radius, v) = mid_radius;
			}
		}
		return true;
	});
}

void resample_branch(Graph& g, std::pair<Dart, Dart> g_branch, std::pair<Scalar, Scalar> g_branch_arclength,
					 Graph& new_g, Graph::Edge new_g_edge, Graph::Attribute<Vec3>* g_vertex_position,
					 Graph::Attribute<Scalar>* g_vertex_radius, Graph::Attribute<Vec3>* new_g_vertex_position,
					 Graph::Attribute<Scalar>* new_g_vertex_radius, Graph::Attribute<Scalar>* new_g_vertex_arclength,
					 Scalar density)
{
	Scalar branch_length = internal::branch_length(g, g_branch, g_vertex_position);

	// if the radius of the current extremities are disjoint, bisect
	if (value<Scalar>(g, g_vertex_radius, Graph::Vertex(g_branch.first)) +
			value<Scalar>(g, g_vertex_radius, Graph::Vertex(g_branch.second)) <
		density * branch_length)
	{
		auto mid_location = internal::arclength_value(g, g_branch, g_vertex_position, branch_length, 0.5);
		Dart e1 = mid_location.first;
		Dart e2 = alpha0(g, e1);
		Scalar alpha = mid_location.second;
		Vec3 mid_pos = alpha * value<Vec3>(g, g_vertex_position, Graph::Vertex(e1)) +
					   (1.0 - alpha) * value<Vec3>(g, g_vertex_position, Graph::Vertex(e2));
		Scalar mid_radius = alpha * value<Scalar>(g, g_vertex_radius, Graph::Vertex(e1)) +
							(1.0 - alpha) * value<Scalar>(g, g_vertex_radius, Graph::Vertex(e2));

		// insert a vertex in the graph at this location
		Graph::Vertex g_v = cut_edge(g, Graph::Edge(e1));
		value<Vec3>(g, g_vertex_position, g_v) = mid_pos;
		value<Scalar>(g, g_vertex_radius, g_v) = mid_radius;

		// insert a vertex in the new graph at this location
		Graph::Vertex new_g_v = cut_edge(new_g, new_g_edge);
		value<Vec3>(new_g, new_g_vertex_position, new_g_v) = mid_pos;
		value<Scalar>(new_g, new_g_vertex_radius, new_g_v) = mid_radius;

		// recursive calls on the right and left branches
		Scalar g_v_arclength = g_branch_arclength.first + (g_branch_arclength.second - g_branch_arclength.first) * 0.5;
		value<Scalar>(new_g, new_g_vertex_arclength, new_g_v) = g_v_arclength;
		resample_branch(g, {g_branch.first, g_v.dart_}, {g_branch_arclength.first, g_v_arclength}, new_g, new_g_edge,
						g_vertex_position, g_vertex_radius, new_g_vertex_position, new_g_vertex_radius,
						new_g_vertex_arclength, density);
		resample_branch(g, {alpha1(g, g_v.dart_), g_branch.second}, {g_v_arclength, g_branch_arclength.second}, new_g,
						Graph::Edge(alpha1(new_g, new_g_v.dart_)), g_vertex_position, g_vertex_radius,
						new_g_vertex_position, new_g_vertex_radius, new_g_vertex_arclength, density);
	}
}

////////////////////
// IncidenceGraph //
////////////////////

void resample_graph(IncidenceGraph& g, IncidenceGraph::Attribute<Vec3>* g_vertex_position,
					IncidenceGraph::Attribute<Scalar>* g_vertex_radius, IncidenceGraph& new_g,
					IncidenceGraph::Attribute<Vec3>* new_g_vertex_position,
					IncidenceGraph::Attribute<Scalar>* new_g_vertex_radius, Scalar density)
{
	IncidenceGraphData g_data;
	get_incidenceGraph_data(g, g_data);
}

} // namespace modeling

} // namespace cgogn
