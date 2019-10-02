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

#include <cgogn/core/functions/mesh_ops/vertex.h>
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

Graph::Vertex
add_vertex(Graph& g, bool set_indices)
{
	Dart d = g.add_dart();
	Dart dd = g.add_dart();
	g.alpha0_sew(d, dd);
	g.set_boundary(dd, true);

	Graph::Vertex v(d);

	if (set_indices)
	{
		if (g.is_embedded<Graph::Vertex>())
			create_embedding(g, v);
		if (g.is_embedded<Graph::Edge>())
			create_embedding(g, Graph::Edge(v.dart));
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

void
remove_vertex(Graph& g, Graph::Vertex v, bool set_indices)
{
	Dart dd = g.alpha0(v.dart);
	cgogn_message_assert(g.is_boundary(dd), "Vertex is still connected to another vertex");
	g.remove_dart(v.dart);
	g.remove_dart(dd);

	if (set_indices)
	{}
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// connect_vertices(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Edge
connect_vertices(Graph& g, Graph::Vertex v1, Graph::Vertex v2, bool set_indices)
{
	Dart d = v1.dart;
	Dart e = v2.dart;
	Dart dd = g.alpha0(d);
	Dart ee = g.alpha0(e);
	if (g.is_boundary(dd))
	{
		if (g.is_boundary(ee))
		{
			g.remove_dart(dd);
			g.remove_dart(ee);
			g.alpha0_sew(d, e);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Edge>())
					g.copy_embedding<Graph::Edge>(e, d);
			}
			return Graph::Edge(d);
		}
		else
		{
			g.set_boundary(dd, false);
			g.alpha1_sew(e, dd);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Vertex>())
					g.copy_embedding<Graph::Vertex>(dd, e);
			}
			return Graph::Edge(d);
		}
	}
	else
	{
		if (g.is_boundary(ee))
		{
			g.set_boundary(ee, false);
			g.alpha1_sew(d, ee);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Vertex>())
					g.copy_embedding<Graph::Vertex>(ee, d);	
			}
			return Graph::Edge(ee);
		}
		else
		{
			Dart dd = g.add_dart();
			Dart ee = g.add_dart();
			g.alpha1_sew(d, dd);
			g.alpha1_sew(e, ee);
			g.alpha0_sew(dd, ee);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Vertex>())
                {
                    g.copy_embedding<Graph::Vertex>(dd, d);
                    g.copy_embedding<Graph::Vertex>(ee, e);
				}
				if (g.is_embedded<Graph::Edge>())
                    create_embedding(g, Graph::Edge(dd));
			}
			return Graph::Edge(dd);
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

void
disconnect_vertices(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart x = e.dart;
	Dart y = g.alpha0(x);
	cgogn_message_assert(!(g.is_boundary(x) || g.is_boundary(y)), "Given edge does not connect 2 vertices");
	if (g.alpha1(x) == x)
	{
		if (g.alpha1(y) == y)
		{
			g.alpha0_unsew(x);
			Dart xx = g.add_dart();
			Dart yy = g.add_dart();
			g.alpha0_sew(x, xx);
			g.alpha0_sew(y, yy);
			g.set_boundary(xx, true);
			g.set_boundary(yy, true);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Edge>())
				{
					g.copy_embedding<Graph::Edge>(g.alpha0(x), x);
					create_embedding(g, Graph::Edge(y));
				}
			}
		}
		else
		{
			g.alpha1_unsew(y);
			g.set_boundary(y, true);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Vertex>())
					g.unset_embedding<Graph::Vertex>(g.alpha0(x));
			}
		}
	}
	else
	{
		if (g.alpha1(y) == y)
		{
			g.alpha1_unsew(x);
			g.set_boundary(x, true);
			if (set_indices)
			{
				if (g.is_embedded<Graph::Vertex>())
					g.unset_embedding<Graph::Vertex>(g.alpha0(y));
			}
		}
		else
		{
			g.alpha0_unsew(x);
			g.alpha1_unsew(x);
			g.alpha1_unsew(y);
			g.remove_dart(x);
			g.remove_dart(y);
			if (set_indices)
			{}
		}
	}
}

} // namespace cgogn
