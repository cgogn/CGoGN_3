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
#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/cmapbase_infos.h>

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
    Dart d = g.mesh().add_dart();
    Dart dd = g.mesh().add_dart();
	g.alpha0_sew(d, dd);
	g.alpha1_sew(d, dd);

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

void
remove_vertex(Graph& g, Graph::Vertex v, bool set_indices)
{
	Dart dd = g.alpha0(v.dart);
	cgogn_message_assert(is_boundary(g,dd), "Vertex is still connected to another vertex");
    g.mesh().remove_dart(v.dart);
    g.mesh().remove_dart(dd);

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
	if (g.is_isolated(v1))
	{
		if (g.is_isolated(v2))
		{
			g.alpha1_unsew(d);
			g.alpha1_unsew(e);
            g.mesh().remove_dart(dd);
            g.mesh().remove_dart(ee);
			g.alpha0_sew(d, e);
			if (set_indices)
			{
                if (is_indexed<Graph::Edge>(g))
					copy_index<Graph::Edge>(g,e, d);
			}
			return Graph::Edge(d);
		}
		else
		{
			g.alpha1_unsew(d);
			g.alpha1_sew(e, dd);
			if (set_indices)
			{
                if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g,dd, e);
			}
			return Graph::Edge(d);
		}
	}
	else
	{
		if (g.is_isolated(v2))
		{
			g.alpha1_unsew(e);
			g.alpha1_sew(d, ee);
			if (set_indices)
			{
                if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g,ee, d);	
			}
			return Graph::Edge(ee);
		}
		else
		{
            Dart dd = g.mesh().add_dart();
            Dart ee = g.mesh().add_dart();
			g.alpha0_sew(dd, ee);
			g.alpha1_sew(d, dd);
			g.alpha1_sew(e, ee);
			if (set_indices)
			{
                if (is_indexed<Graph::Vertex>(g))
                {
                    copy_index<Graph::Vertex>(g,dd, d);
                    copy_index<Graph::Vertex>(g,ee, e);
				}
                if (is_indexed<Graph::HalfEdge>(g))
                {
					set_index(g, Graph::HalfEdge(dd), new_index<Graph::HalfEdge>(g));
					set_index(g, Graph::HalfEdge(ee), new_index<Graph::HalfEdge>(g));
				}
                if (is_indexed<Graph::Edge>(g))
					set_index(g, Graph::Edge(dd), new_index<Graph::Edge>(g));
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
	cgogn_message_assert(!(g.is_isolated(Graph::Vertex(x)) || g.is_isolated(Graph::Vertex(y))), "Given edge does not connect 2 vertices");
	if (g.alpha1(x) == x)
	{
		if (g.alpha1(y) == y)
		{
			g.alpha0_unsew(x);
            Dart xx = g.mesh().add_dart();
            Dart yy = g.mesh().add_dart();
			g.alpha0_sew(x, xx);
			g.alpha1_sew(x, xx);
			g.alpha0_sew(y, yy);
			g.alpha1_sew(y, yy);
			if (set_indices)
			{
                if (is_indexed<Graph::Vertex>(g))
				{
					copy_index<Graph::Edge>(g,xx, x);
					copy_index<Graph::Edge>(g,yy, y);
				}
                if (is_indexed<Graph::HalfEdge>(g))
                {
					set_index(g, Graph::HalfEdge(xx), new_index<Graph::HalfEdge>(g));
					set_index(g, Graph::HalfEdge(yy), new_index<Graph::HalfEdge>(g));
				}
                if (is_indexed<Graph::Edge>(g))
				{
					copy_index<Graph::Edge>(g,xx, x);
					set_index(g, Graph::Edge(y), new_index<Graph::Edge>(g));
				}
			}
		}
		else
		{
			g.alpha1_unsew(y);
			g.alpha1_sew(x, y);
			if (set_indices)
			{
                if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g,y, x);
			}
		}
	}
	else
	{
		if (g.alpha1(y) == y)
		{
			g.alpha1_unsew(x);
			g.alpha1_sew(y, x);
			if (set_indices)
			{
                if (is_indexed<Graph::Vertex>(g))
					copy_index<Graph::Vertex>(g,x, y);
			}
		}
		else
		{
			g.alpha0_unsew(x);
			g.alpha1_unsew(x);
			g.alpha1_unsew(y);
            g.mesh().remove_dart(x);
            g.mesh().remove_dart(y);
			if (set_indices)
			{}
		}
	}
}

} // namespace cgogn
