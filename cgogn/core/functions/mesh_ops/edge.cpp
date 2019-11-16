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

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// cut_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// Graph //
///////////

Graph::Vertex
cut_edge(Graph& g, Graph::Edge e, bool set_indices)
{
	Dart e0 = e.dart;
	Dart e1 = g.alpha0(e0);

	Dart v0 = g.add_dart();
	Dart v1 = g.add_dart();

	g.alpha1_sew(v0, v1);
	g.alpha0_unsew(e0);
	g.alpha0_sew(e0, v0);
	g.alpha0_sew(e1, v1);

	if (set_indices)
	{
		if (g.is_indexed<Graph::Vertex>())
			set_index(g, Graph::Vertex(v0), new_index<Graph::Vertex>(g));
		if (g.is_indexed<Graph::HalfEdge>())
		{
			set_index(g, Graph::HalfEdge(v0), new_index<Graph::HalfEdge>(g));
			set_index(g, Graph::HalfEdge(v1), new_index<Graph::HalfEdge>(g));
		}
		if (g.is_indexed<Graph::Edge>())
		{
			g.copy_index<Graph::Edge>(v0, e0);
			set_index(g, Graph::Edge(e1), new_index<Graph::Edge>(g));
		}
	}

	return Graph::Vertex(v0);
}

///////////
// CMap1 //
///////////

CMap1::Vertex
cut_edge(CMap1& m, CMap1::Edge e, bool set_indices)
{
	Dart d = m.add_dart();
	m.phi1_sew(e.dart, d);
	CMap1::Vertex v(d);

	if (set_indices)
	{
		if (m.is_indexed<CMap1::Vertex>())
			set_index(m, v, new_index<CMap1::Vertex>(m));
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (m.is_indexed<CMap1::Face>())
			m.copy_index<CMap1::Face>(d, e.dart);
	}

	return v;
}

///////////
// CMap2 //
///////////

CMap2::Vertex
cut_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart d1 = e.dart;
	Dart d2 = m.phi2(d1);
	m.phi2_unsew(d1);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d1), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d2), false);
	m.phi2_sew(d1, nv2.dart);
	m.phi2_sew(d2, nv1.dart);
	m.set_boundary(nv1.dart, m.is_boundary(d1));
	m.set_boundary(nv2.dart, m.is_boundary(d2));
	CMap2::Vertex v(nv1.dart);

	if (set_indices)
	{
		if (m.is_indexed<CMap2::Vertex>())
			set_index(m, v, new_index<CMap2::Vertex>(m));
		if (m.is_indexed<CMap2::HalfEdge>())
		{
			set_index(m, CMap2::HalfEdge(nv1.dart), new_index<CMap2::HalfEdge>(m));
			set_index(m, CMap2::HalfEdge(nv2.dart), new_index<CMap2::HalfEdge>(m));
		}
		if (m.is_indexed<CMap2::Edge>())
		{
			m.copy_index<CMap2::Edge>(nv2.dart, d1);
			set_index(m, CMap2::Edge(nv1.dart), new_index<CMap2::Edge>(m));
		}
		if (m.is_indexed<CMap2::Face>())
		{
			m.copy_index<CMap2::Face>(nv1.dart, d1);
			m.copy_index<CMap2::Face>(nv2.dart, d2);
		}
		if (m.is_indexed<CMap2::Volume>())
		{
			m.copy_index<CMap2::Volume>(nv1.dart, d1);
			m.copy_index<CMap2::Volume>(nv2.dart, d2);
		}
	}

	return v;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Vertex
// collapse_edge(MESH& m, typename mesh_traits<MESH>::Edge e, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

CMap1::Vertex
collapse_edge(CMap1& m, CMap1::Edge e, bool set_indices)
{
	Dart d = m.phi_1(e.dart);
	m.phi1_unsew(d);
	m.remove_dart(e.dart);
	CMap1::Vertex v(m.phi1(d));

	if (set_indices)
	{}

	return v;
}

///////////
// CMap2 //
///////////

CMap2::Vertex
collapse_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart dd = e.dart;
	Dart dd_1 = m.phi_1(dd);
	Dart dd_12 = m.phi2(dd_1);
	Dart ee = m.phi2(dd);
	Dart ee_1 = m.phi_1(ee);
	Dart ee_12 = m.phi2(ee_1);

	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);

	if (codegree(m, CMap2::Face(dd_1)) == 2u)
	{
		Dart dd1 = m.phi1(dd_1);
		Dart dd12 = m.phi2(dd1);
		m.phi2_unsew(dd1);
		m.phi2_unsew(dd_1);
		m.phi2_sew(dd12, dd_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(dd1), false);
	}

	if (codegree(m, CMap2::Face(ee_1)) == 2u)
	{
		Dart ee1 = m.phi1(ee_1);
		Dart ee12 = m.phi2(ee1);
		m.phi2_unsew(ee1);
		m.phi2_unsew(ee_1);
		m.phi2_sew(ee12, ee_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(ee1), false);
	}

	CMap2::Vertex v(dd_12);

	if (set_indices)
	{
		if (m.is_indexed<CMap2::Vertex>())
			set_index(m, v, index_of(m, v));
		if (m.is_indexed<CMap2::Edge>())
		{
			m.copy_index<CMap2::Edge>(dd_12, m.phi2(dd_12));
			m.copy_index<CMap2::Edge>(ee_12, m.phi2(ee_12));
		}
	}

	return v;
}

} // namespace cgogn
