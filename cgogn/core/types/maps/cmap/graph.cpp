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

#include <cgogn/core/types/maps/cmap/graph.h>

namespace cgogn
{

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

void merge_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices)
{
	alpha1_sew(g, v1.dart, v2.dart);

	if (set_indices)
	{
		if (is_indexed<Graph::Vertex>(g))
			set_index(g, v1, index_of(g, v1));
	}
}

Graph::Edge connect_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices)
{
	static auto is_isolated = [](Graph& g, Graph::Vertex v) -> bool { return alpha0(g, v.dart) == alpha1(g, v.dart); };

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

void disconnect_vertices(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart x = e.dart;
	Dart y = alpha0(g, x);
	cgogn_message_assert(!(is_vertex_isolated(g, Graph::Vertex(x)) || is_vertex_isolated(g, Graph::Vertex(y))),
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

Graph::Vertex cut_edge(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart e0 = e.dart;
	Dart e1 = alpha0(g, e0);

	Dart v0 = add_dart(g);
	Dart v1 = add_dart(g);

	alpha1_sew(g, v0, v1);
	alpha0_unsew(g, e0);
	alpha0_sew(g, e0, v0);
	alpha0_sew(g, e1, v1);

	if (set_indices)
	{
		if (is_indexed<Graph::Vertex>(g))
			set_index(g, Graph::Vertex(v0), new_index<Graph::Vertex>(g));
		if (is_indexed<Graph::HalfEdge>(g))
		{
			set_index(g, Graph::HalfEdge(v0), new_index<Graph::HalfEdge>(g));
			set_index(g, Graph::HalfEdge(v1), new_index<Graph::HalfEdge>(g));
		}
		if (is_indexed<Graph::Edge>(g))
		{
			copy_index<Graph::Edge>(g, v0, e0);
			set_index(g, Graph::Edge(e1), new_index<Graph::Edge>(g));
		}
	}

	return Graph::Vertex(v0);
}

Graph::Vertex collapse_edge(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart d = e.dart;
	Dart d1 = alpha_1(g, d);
	Dart dd = alpha0(g, d);
	Dart dd1 = alpha_1(g, dd);

	Graph::Vertex v;

	if (dd1 != dd)
	{
		v.dart = dd1;
		alpha1_unsew(g, dd);
	}
	if (d1 != d)
	{
		v.dart = d1;
		alpha1_unsew(g, d);
	}
	if (d1 != d && dd1 != dd)
		alpha1_sew(g, d1, dd1);
	remove_dart(g, d);
	remove_dart(g, dd);

	if (set_indices && v.is_valid())
	{
		if (is_indexed<Graph::Vertex>(g))
			set_index(g, v, index_of(g, v));
	}

	return v;
}

void alpha0_sew(Graph& m, Dart d, Dart e)
{
	(*m.alpha0_)[d.index] = e;
	(*m.alpha0_)[e.index] = d;
}

void alpha0_unsew(Graph& m, Dart d)
{
	Dart e = alpha0(m, d);
	(*m.alpha0_)[d.index] = d;
	(*m.alpha0_)[e.index] = e;
}

void alpha1_sew(Graph& m, Dart d, Dart e)
{
	Dart f = alpha1(m, d);
	Dart g = alpha1(m, e);
	(*m.alpha1_)[d.index] = g;
	(*m.alpha1_)[e.index] = f;
	(*m.alpha_1_)[g.index] = d;
	(*m.alpha_1_)[f.index] = e;
}

void alpha1_unsew(Graph& m, Dart d)
{
	Dart e = alpha1(m, d);
	Dart f = alpha_1(m, d);
	(*m.alpha1_)[f.index] = e;
	(*m.alpha1_)[d.index] = d;
	(*m.alpha_1_)[e.index] = f;
	(*m.alpha_1_)[d.index] = d;
}

} // namespace cgogn
