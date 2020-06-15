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

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>
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

void compute_graph_radius_from_surface(Graph& g, Graph::Attribute<Vec3>* g_vertex_position,
									   Graph::Attribute<Scalar>* g_vertex_radius, const CMap2& s,
									   const CMap2::Attribute<Vec3>* s_vertex_position)
{
	using SelectedFace = std::tuple<CMap2::Face, Vec3, Scalar>;
	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		const Vec3& p = value<Vec3>(g, g_vertex_position, v);
		Vec3 cp = geometry::closest_point_on_surface(s, s_vertex_position, p);
		value<Scalar>(g, g_vertex_radius, v) = (cp - p).norm();
		return true;
	});
}

void resample_graph(Graph& g, Graph::Attribute<Vec3>* g_vertex_position, Graph::Attribute<Scalar>* g_vertex_radius,
					Graph& new_g, Graph::Attribute<Vec3>* new_g_vertex_position,
					Graph::Attribute<Scalar>* new_g_vertex_radius)
{
	GraphData g_data;
	get_graph_data(g, g_data);

	auto g_vertex_new_g_vertices = add_attribute<std::vector<Dart>, Graph::Vertex>(g, "new_g_vertices");

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
		value<std::vector<Dart>>(g, g_vertex_new_g_vertices, Graph::Vertex(branch.first.dart)).push_back(ed);
		value<std::vector<Dart>>(g, g_vertex_new_g_vertices, Graph::Vertex(branch.second.dart)).push_back(edd);

		// set the position of the new graph vertices
		value<Vec3>(new_g, new_g_vertex_position, v1) =
			value<Vec3>(g, g_vertex_position, Graph::Vertex(branch.first.dart));
		value<Vec3>(new_g, new_g_vertex_position, v2) =
			value<Vec3>(g, g_vertex_position, Graph::Vertex(branch.second.dart));
		// set the radius of the new graph vertices
		value<Scalar>(new_g, new_g_vertex_radius, v1) =
			value<Scalar>(g, g_vertex_radius, Graph::Vertex(branch.first.dart));
		value<Scalar>(new_g, new_g_vertex_radius, v2) =
			value<Scalar>(g, g_vertex_radius, Graph::Vertex(branch.second.dart));

		resample_branch(g, {branch.first.dart, branch.second.dart}, new_g, Graph::Edge(ed), g_vertex_position,
						g_vertex_radius, new_g_vertex_position, new_g_vertex_radius);
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

void resample_branch(Graph& g, std::pair<Dart, Dart> g_branch, Graph& new_g, Graph::Edge new_g_edge,
					 Graph::Attribute<Vec3>* g_vertex_position, Graph::Attribute<Scalar>* g_vertex_radius,
					 Graph::Attribute<Vec3>* new_g_vertex_position, Graph::Attribute<Scalar>* new_g_vertex_radius)
{
	// compute the length of the branch
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

	// if the radius of the current extremities are disjoint, bisect
	if (value<Scalar>(g, g_vertex_radius, Graph::Vertex(g_branch.first)) +
			value<Scalar>(g, g_vertex_radius, Graph::Vertex(g_branch.second)) <
		0.9 * branch_length)
	{
		// find the edge in which the 0.5 value lies
		Scalar cur_length = 0;
		Scalar next_length = 0;
		cur = g_branch.first;
		next = alpha0(g, cur);
		do
		{
			cur_length = next_length;
			next_length += (value<Vec3>(g, g_vertex_position, Graph::Vertex(next)) -
							value<Vec3>(g, g_vertex_position, Graph::Vertex(cur)))
							   .norm() /
						   branch_length;
			if (next_length < 0.5)
			{
				cur = alpha1(g, next);
				next = alpha0(g, cur);
			}
		} while (next_length < 0.5);
		// find where in this edge the 0.5 value lies
		Scalar alpha = (0.5 - cur_length) / (next_length - cur_length);
		Vec3 mid_pos = alpha * value<Vec3>(g, g_vertex_position, Graph::Vertex(next)) +
					   (1.0 - alpha) * value<Vec3>(g, g_vertex_position, Graph::Vertex(cur));
		Scalar mid_radius = alpha * value<Scalar>(g, g_vertex_radius, Graph::Vertex(next)) +
							(1.0 - alpha) * value<Scalar>(g, g_vertex_radius, Graph::Vertex(cur));
		// insert a vertex in the graph at this location
		Graph::Vertex g_v = cut_edge(g, Graph::Edge(cur));
		value<Vec3>(g, g_vertex_position, g_v) = mid_pos;
		value<Scalar>(g, g_vertex_radius, g_v) = mid_radius;
		// insert a vertex in the new graph at this location
		Graph::Vertex new_g_v = cut_edge(new_g, new_g_edge);
		value<Vec3>(new_g, new_g_vertex_position, new_g_v) = mid_pos;
		value<Scalar>(new_g, new_g_vertex_radius, new_g_v) = mid_radius;
		// recursive calls on the right and left branches
		resample_branch(g, {g_v.dart, g_branch.first}, new_g, Graph::Edge(new_g_v.dart), g_vertex_position,
						g_vertex_radius, new_g_vertex_position, new_g_vertex_radius);
		resample_branch(g, {alpha1(g, g_v.dart), g_branch.second}, new_g, Graph::Edge(alpha1(new_g, new_g_v.dart)),
						g_vertex_position, g_vertex_radius, new_g_vertex_position, new_g_vertex_radius);
	}
}

} // namespace modeling

} // namespace cgogn