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
#include <cgogn/core/functions/mr_cmap_infos.h>

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

	Dart v0 = add_dart(g);
	Dart v1 = add_dart(g);

	g.alpha1_sew(v0, v1);
	g.alpha0_unsew(e0);
	g.alpha0_sew(e0, v0);
	g.alpha0_sew(e1, v1);

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
			copy_index<Graph::Edge>(g,v0, e0);
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
	Dart d = add_dart(m);
	phi1_sew(m,e.dart, d);
	CMap1::Vertex v(d);

	if (set_indices)
	{
		if (is_indexed<CMap1::Vertex>(m))
			set_index(m, v, new_index<CMap1::Vertex>(m));
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (is_indexed<CMap1::Face>(m))
			copy_index<CMap1::Face>(m,d, e.dart);
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
	Dart d2 = phi2(m,d1);
	phi2_unsew(m,d1);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d1), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d2), false);
	phi2_sew(m,d1, nv2.dart);
	phi2_sew(m,d2, nv1.dart);
	set_boundary(m,nv1.dart, is_boundary(m,d1));
	set_boundary(m,nv2.dart, is_boundary(m,d2));
	CMap2::Vertex v(nv1.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
			set_index(m, v, new_index<CMap2::Vertex>(m));
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			set_index(m, CMap2::HalfEdge(nv1.dart), new_index<CMap2::HalfEdge>(m));
			set_index(m, CMap2::HalfEdge(nv2.dart), new_index<CMap2::HalfEdge>(m));
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			copy_index<CMap2::Edge>(m,nv2.dart, d1);
			set_index(m, CMap2::Edge(nv1.dart), new_index<CMap2::Edge>(m));
		}
		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m,nv1.dart, d1);
			copy_index<CMap2::Face>(m,nv2.dart, d2);
		}
		if (is_indexed<CMap2::Volume>(m))
		{
			copy_index<CMap2::Volume>(m,nv1.dart, d1);
			copy_index<CMap2::Volume>(m,nv2.dart, d2);
		}
	}

	return v;
}

///////////
// CMap3 //
///////////

CMap3::Vertex
cut_edge(CMap3& m, CMap3::Edge e, bool set_indices)
{
	Dart d0 = e.dart;
	Dart d23 = phi<23>(m,d0);
	CMap3::Vertex v(cut_edge(static_cast<CMap2&>(m), CMap2::Edge(d0), false).dart);

	while(d23 != e.dart)
	{
		d0 = d23;
		d23 = phi<23>(m,d23);

		cut_edge(static_cast<CMap2&>(m), CMap2::Edge(d0), false);
		
		const Dart d3 = phi3(m,d0);
		phi3_unsew(m,d0);

		phi3_sew(m,d0, phi1(m,d3));
		phi3_sew(m,d3, phi1(m,d0));
	}

	const Dart d3 = phi3(m,e.dart);
	phi3_unsew(m,e.dart);

	phi3_sew(m,e.dart, phi1(m,d3));
	phi3_sew(m,d3, phi1(m,e.dart));

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
			set_index(m, v, new_index<CMap3::Vertex>(m));

		if (is_indexed<CMap3::Edge>(m))
		{
			set_index(m, CMap3::Edge(v.dart), new_index<CMap3::Edge>(m));
			set_index(m, e, index_of(m, e));
		}
		if (is_indexed<CMap3::Face>(m))
		{
			foreach_dart_of_PHI23(m,e.dart, [&](Dart d) -> bool
			{
				copy_index<CMap3::Face>(m,phi1(m,d), d);
				copy_index<CMap3::Face>(m,phi3(m,d), d);
				copy_index<CMap3::Face>(m,phi2(m,d), phi<12>(m,d));
				return true;
			});
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			foreach_dart_of_orbit(m,e, [&](Dart d) -> bool
			{
				if(is_boundary(m,d))
					return true;
				copy_index<CMap3::Volume>(m,phi1(m,d), d);
				copy_index<CMap3::Volume>(m,phi2(m,d), d);
				return true;
			});
		}
	}

	return v;
}

/////////////
// MRCmap3 //
/////////////

MRCmap3::Vertex
cut_edge(MRCmap3& m, MRCmap3::Edge e, bool set_indices)
{
	MRCmap3::Vertex v = cut_edge(m.mesh(),e,false);
	foreach_dart_of_PHI23(m,e.dart, [&](Dart d) -> bool
	{
		m.cph().edge_id(phi1(m,d),m.cph().edge_id(d));
		m.cph().edge_id(phi3(m,d),m.cph().edge_id(d));
		m.cph().edge_id(phi2(m,d),m.cph().edge_id(phi<12>(m,d)));
		m.cph().face_id(phi1(m,d),m.cph().face_id(d));
		m.cph().face_id(phi3(m,d),m.cph().face_id(d));
		m.cph().face_id(phi2(m,d),m.cph().face_id(phi<12>(m,d)));
		m.cph().change_dart_level(phi1(m,d),m.current_level());
		m.cph().change_dart_level(phi3(m,d),m.current_level());
		m.cph().change_dart_level(phi2(m,d),m.current_level());
		return true;
	});
	
	if (set_indices)
	{
		if (is_indexed<MRCmap3::Vertex>(m))
			set_index(m, v, new_index<MRCmap3::Vertex>(m));

		if (is_indexed<MRCmap3::Edge>(m))
		{
			uint32 ne = new_index<MRCmap3::Edge>(m);
			foreach_dart_of_orbit(m,e,[&](Dart d)->bool{
				if(m.cph().dart_level(d) == m.current_level())
					set_index<MRCmap3::Edge>(m,d,ne);
				return true;
			});
			ne = new_index<MRCmap3::Edge>(m);
			foreach_dart_of_orbit(m,MRCmap3::Edge(phi1(m,e.dart)),[&](Dart d)->bool{
				if(m.cph().dart_level(d) == m.current_level())
					set_index<MRCmap3::Edge>(m,d,ne);
				return true;
			});
		}
		if (is_indexed<MRCmap3::Face>(m))
		{
			foreach_dart_of_PHI23(m,e.dart, [&](Dart d) -> bool
			{
				Dart it = phi1(m,d);
                do{
                    it = phi1(m,it);
                }while(m.cph().dart_level(it) < m.current_level()-1 && it != d);
				copy_index<MRCmap3::Face>(m,phi1(m,d), it);
                it = phi2(m,d);
                do{
                    it = phi1(m,it);
                }while(m.cph().dart_level(it) < m.current_level()-1 && it != phi2(m,phi1(m,d)));
				copy_index<MRCmap3::Face>(m,phi2(m,d), it);
				return true;
			});
		}
		if (is_indexed<MRCmap3::Volume>(m))
		{
			foreach_dart_of_PHI23(m,e.dart, [&](Dart d) -> bool
			{
				if (!is_boundary(m,d))
                {
                    Dart it = phi1(m,d);
                    do{
                        it = phi1(m,it);
                    }while(m.cph().dart_level(it) < m.current_level()-1 && it != d);
					copy_index<MRCmap3::Volume>(m,phi1(m,d), it);
                    it = phi2(m,d);
                    do{
                        it = phi1(m,it);
                    }while(m.cph().dart_level(it) < m.current_level()-1 && it != phi2(m,phi1(m,d)));
					copy_index<MRCmap3::Volume>(m,phi2(m,d), it);
                }
				return true;
			});
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
	Dart d = phi_1(m,e.dart);
	phi1_unsew(m,d);
	m.mesh().remove_dart(e.dart);
	CMap1::Vertex v(phi1(m,d));

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
	Dart dd_1 = phi_1(m,dd);
	Dart dd_12 = phi2(m,dd_1);
	Dart ee = phi2(m,dd);
	Dart ee_1 = phi_1(m,ee);
	Dart ee_12 = phi2(m,ee_1);

	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);

	if (codegree(m, CMap2::Face(dd_1)) == 2u)
	{
		Dart dd1 = phi1(m,dd_1);
		Dart dd12 = phi2(m,dd1);
		phi2_unsew(m,dd1);
		phi2_unsew(m,dd_1);
		phi2_sew(m,dd12, dd_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(dd1), false);
	}

	if (codegree(m, CMap2::Face(ee_1)) == 2u)
	{
		Dart ee1 = phi1(m,ee_1);
		Dart ee12 = phi2(m,ee1);
		phi2_unsew(m,ee1);
		phi2_unsew(m,ee_1);
		phi2_sew(m,ee12, ee_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(ee1), false);
	}

	CMap2::Vertex v(dd_12);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
			set_index(m, v, index_of(m, v));
		if (is_indexed<CMap2::Edge>(m))
		{
			copy_index<CMap2::Edge>(m,dd_12, phi2(m,dd_12));
			copy_index<CMap2::Edge>(m,ee_12, phi2(m,ee_12));
		}
	}

	return v;
}

} // namespace cgogn
