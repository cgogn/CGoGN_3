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

#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/core/types/cmap/cmap_info.h>
#include <cgogn/core/types/cmap/cmap_ops.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// add_vertex(MESH& m, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Vertex add_vertex(Graph& g, bool set_indices)
{
	Dart d = add_dart(g);
	Dart dd = add_dart(g);
	alpha0_sew(g, d, dd);
	alpha1_sew(g, d, dd);

	Graph::Vertex v(d);

	if (set_indices)
	{
		if (is_indexed<Graph::Vertex>(g))
			set_index(g, v, new_index<Graph::Vertex>(g));
		if (is_indexed<Graph::HalfEdge>(g))
		{
			set_index(g, Graph::HalfEdge(d), new_index<Graph::HalfEdge>(g));
			set_index(g, Graph::HalfEdge(dd), new_index<Graph::HalfEdge>(g));
		}
		if (is_indexed<Graph::Edge>(g))
			set_index(g, Graph::Edge(v.dart), new_index<Graph::Edge>(g));
	}

	return v;
}

/*****************************************************************************/

// template <typename MESH>
// void
// remove_vertex(MESH& m, typename mesh_traits<MESH>::Vertex v, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

void remove_vertex(Graph& g, Graph::Vertex v, bool set_indices)
{
	Dart dd = alpha0(g, v.dart);
	cgogn_message_assert(is_boundary(g, dd), "Vertex is still connected to another vertex");
	remove_dart(g, v.dart);
	remove_dart(g, dd);

	if (set_indices)
	{
	}
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// connect_vertices(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool
// set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Edge connect_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices)
{
	auto is_isolated = [](Graph& g, Graph::Vertex v) -> bool { return alpha0(g, v.dart) == alpha1(g, v.dart); };

	Dart d = v1.dart;
	Dart e = v2.dart;
	Dart dd = alpha0(g, d);
	Dart ee = alpha0(g, e);
	if (is_isolated(g, v1))
	{
		if (is_isolated(g, v2))
		{
			alpha1_unsew(g, d);
			alpha1_unsew(g, e);
			remove_dart(g, dd);
			remove_dart(g, ee);
			alpha0_sew(g, d, e);
			if (set_indices)
			{
				if (is_indexed<Graph::Edge>(g))
					copy_index<Graph::Edge>(g, e, d);
			}
			return Graph::Edge(d);
		}
		else
		{
			alpha1_unsew(g, d);
			alpha1_sew(g, e, dd);
			if (set_indices)
			{
				if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g, dd, e);
			}
			return Graph::Edge(d);
		}
	}
	else
	{
		if (is_isolated(g, v2))
		{
			alpha1_unsew(g, e);
			alpha1_sew(g, d, ee);
			if (set_indices)
			{
				if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g, ee, d);
			}
			return Graph::Edge(ee);
		}
		else
		{
			Dart ddd = add_dart(g);
			Dart eee = add_dart(g);
			alpha0_sew(g, ddd, eee);
			alpha1_sew(g, d, ddd);
			alpha1_sew(g, e, eee);
			if (set_indices)
			{
				if (is_indexed<Graph::Vertex>(g))
				{
					copy_index<Graph::Vertex>(g, ddd, d);
					copy_index<Graph::Vertex>(g, eee, e);
				}
				if (is_indexed<Graph::HalfEdge>(g))
				{
					set_index(g, Graph::HalfEdge(ddd), new_index<Graph::HalfEdge>(g));
					set_index(g, Graph::HalfEdge(eee), new_index<Graph::HalfEdge>(g));
				}
				if (is_indexed<Graph::Edge>(g))
					set_index(g, Graph::Edge(ddd), new_index<Graph::Edge>(g));
			}
			return Graph::Edge(ddd);
		}
	}
}

/*****************************************************************************/

// template <typename MESH>
// void
// disconnect_vertices(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

void disconnect_vertices(Graph& g, Graph::Edge e, bool set_indices)
{
	auto is_isolated = [](Graph& g, Graph::Vertex v) -> bool { return alpha0(g, v.dart) == alpha1(g, v.dart); };

	Dart x = e.dart;
	Dart y = alpha0(g, x);
	cgogn_message_assert(!(is_isolated(g, Graph::Vertex(x)) || is_isolated(g, Graph::Vertex(y))),
						 "Given edge does not connect 2 vertices");
	if (alpha1(g, x) == x)
	{
		if (alpha1(g, y) == y)
		{
			alpha0_unsew(g, x);
			Dart xx = add_dart(g);
			Dart yy = add_dart(g);
			alpha0_sew(g, x, xx);
			alpha1_sew(g, x, xx);
			alpha0_sew(g, y, yy);
			alpha1_sew(g, y, yy);
			if (set_indices)
			{
				if (is_indexed<Graph::Vertex>(g))
				{
					copy_index<Graph::Edge>(g, xx, x);
					copy_index<Graph::Edge>(g, yy, y);
				}
				if (is_indexed<Graph::HalfEdge>(g))
				{
					set_index(g, Graph::HalfEdge(xx), new_index<Graph::HalfEdge>(g));
					set_index(g, Graph::HalfEdge(yy), new_index<Graph::HalfEdge>(g));
				}
				if (is_indexed<Graph::Edge>(g))
				{
					copy_index<Graph::Edge>(g, xx, x);
					set_index(g, Graph::Edge(y), new_index<Graph::Edge>(g));
				}
			}
		}
		else
		{
			alpha1_unsew(g, y);
			alpha1_sew(g, x, y);
			if (set_indices)
			{
				if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g, y, x);
			}
		}
	}
	else
	{
		if (alpha1(g, y) == y)
		{
			alpha1_unsew(g, x);
			alpha1_sew(g, y, x);
			if (set_indices)
			{
				if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g, x, y);
			}
		}
		else
		{
			alpha0_unsew(g, x);
			alpha1_unsew(g, x);
			alpha1_unsew(g, y);
			remove_dart(g, x);
			remove_dart(g, y);
			if (set_indices)
			{
			}
		}
	}
}

} // namespace cgogn
