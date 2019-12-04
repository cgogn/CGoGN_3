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

#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/phi.h>
#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// add_face(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

CMap1::Face
add_face(CMap1& m, uint32 size, bool set_indices)
{
	Dart d = m.mesh().add_dart();
	for (uint32 i = 1u; i < size; ++i)
	{
		Dart e = m.mesh().add_dart();
		phi1_sew(m,d, e);
	}
	CMap1::Face f(d);

	if (set_indices)
	{
		if (is_indexed<CMap1::Vertex>(m))
		{
			foreach_incident_vertex(m, f, [&] (CMap1::Vertex v) -> bool
			{
				set_index(m, v, new_index<CMap1::Vertex>(m));
				return true;
			});
		}
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (is_indexed<CMap1::Face>(m))
			set_index(m, f, new_index<CMap1::Face>(m));
	}

	return f;
}

///////////
// CMap2 //
///////////

CMap2::Face
add_face(CMap2& m, uint32 size, bool set_indices)
{
	CMap2::Face f = add_face(static_cast<CMap1&>(m), size, false);
	CMap2::Face b = add_face(static_cast<CMap1&>(m), size, false);
	Dart it = b.dart;
	foreach_dart_of_orbit(m,f, [&] (Dart d) -> bool
	{
		set_boundary(m,it, true);
		phi2_sew(m,d, it);
		it = phi_1(m,it);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(m, f, [&] (CMap2::Vertex v) -> bool
			{
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(m, f, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(m, f, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Face>(m))
			set_index(m, f, new_index<CMap2::Face>(m));
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, CMap2::Volume(f.dart), new_index<CMap2::Volume>(m));
	}

	return f;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// remove_face(MESH& m, typename mesh_traits<MESH>::Face f, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap1 //
///////////

void
remove_face(CMap1& m, CMap1::Face f, bool set_indices)
{
	Dart it = phi1(m,f.dart);
	while(it != f.dart)
	{
		Dart next = phi1(m,it);
		m.mesh().remove_dart(it);
		it = next;
	}
	m.mesh().remove_dart(f.dart);

	if (set_indices)
	{}
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Edge
// cut_face(MESH& m, typename mesh_traits<MESH>::Vertex v1, typename mesh_traits<MESH>::Vertex v2, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

CMap2::Edge
cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices)
{
	Dart dd = phi_1(m,v1.dart);
	Dart ee = phi_1(m,v2.dart);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);
	phi1_sew(m,nv1.dart, nv2.dart);
	phi2_sew(m,nv1.dart, nv2.dart);
	set_boundary(m,nv1.dart, is_boundary(m,dd));
	set_boundary(m,nv2.dart, is_boundary(m,ee));
	CMap2::Edge e(nv1.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			copy_index<CMap2::Vertex>(m,nv1.dart, v1.dart);
			copy_index<CMap2::Vertex>(m,nv2.dart, v2.dart);
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			set_index(m, CMap2::HalfEdge(nv1.dart), new_index<CMap2::HalfEdge>(m));
			set_index(m, CMap2::HalfEdge(nv2.dart), new_index<CMap2::HalfEdge>(m));
		}
		if (is_indexed<CMap2::Edge>(m))
			set_index(m, CMap2::Edge(nv1.dart), new_index<CMap2::Edge>(m));
		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m,nv2.dart, v1.dart);
			set_index(m, CMap2::Face(v2.dart), new_index<CMap2::Face>(m));
		}
		if (is_indexed<CMap2::Volume>(m))
		{
			copy_index<CMap2::Volume>(m,nv1.dart, v2.dart);
			copy_index<CMap2::Volume>(m,nv2.dart, v1.dart);
		}
	}

	return e;
}

///////////
// CMap3 //
///////////

CMap3::Edge
cut_face(CMap3& m, CMap3::Vertex v1, CMap3::Vertex v2, bool set_indices)
{
	Dart d = v1.dart;
	Dart e = v2.dart;

	Dart dd = phi<31>(m,v1.dart);
	Dart ee = phi<31>(m,e);

	CMap2::Edge e0 = cut_face(static_cast<CMap2&>(m), CMap2::Vertex(v1.dart), CMap2::Vertex(e), false);
	CMap2::Edge e1 = cut_face(static_cast<CMap2&>(m), CMap2::Vertex(dd), CMap2::Vertex(ee), false);

	phi3_sew(m,phi_1(m,v1.dart), phi_1(m,ee));
	phi3_sew(m,phi_1(m,dd), phi_1(m,e));

	CMap3::Edge edge(phi_1(m,e));

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
		{
			copy_index<CMap3::Vertex>(m,phi_1(m,e), v1.dart);
			copy_index<CMap3::Vertex>(m,phi_1(m,ee), v1.dart);
			copy_index<CMap3::Vertex>(m,phi_1(m,d), e);
			copy_index<CMap3::Vertex>(m,phi_1(m,dd), e);
		}
		if (is_indexed<CMap3::Edge>(m))
			set_index(m, CMap3::Edge(phi_1(m,v1.dart)), new_index<CMap3::Edge>(m));
		if (is_indexed<CMap3::Face>(m))
		{
			copy_index<CMap3::Face>(m,phi_1(m,ee), d);
			copy_index<CMap3::Face>(m,phi_1(m,d), d);
			set_index(m, CMap3::Face(e), new_index<CMap3::Face>(m));
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			copy_index<CMap3::Volume>(m,phi_1(m,d), d);
			copy_index<CMap3::Volume>(m,phi_1(m,e), d);
			copy_index<CMap3::Volume>(m,phi_1(m,dd), dd);
			copy_index<CMap3::Volume>(m,phi_1(m,ee), dd);
		}
	}

	return edge;
}

CMap2::Face close_hole(CMap2& m,Dart d, bool set_indices)
{
	cgogn_message_assert(phi2(m,d) == d, "CMap2: close hole called on a dart that is not a phi2 fix point");

	Dart first = m.mesh().add_dart();	// First edge of the face that will fill the hole
	phi2_sew(m,d, first);			// 2-sew the new edge to the hole

	Dart d_next = d;			// Turn around the hole
	Dart d_phi1;				// to complete the face
	do
	{
		do
		{
			d_phi1 = phi1(m,d_next); // Search and put in d_next
			d_next = phi2(m,d_phi1); // the next dart of the hole
		} while (d_next != d_phi1 && d_phi1 != d);

		if (d_phi1 != d)
		{
			Dart next = m.mesh().add_dart();	// Add a vertex into the built face
			phi1_sew(m,first, next);
			phi2_sew(m,d_next, next);	// and 2-sew the face to the hole
		}
	} while (d_phi1 != d);
	
	CMap2::Face hole(first);

	if (set_indices)
	{
		foreach_dart_of_orbit(m,hole, [&] (Dart hd) -> bool
		{
			Dart hd2 = phi2(m,hd);
			if (is_indexed<CMap2::Vertex>(m))
				copy_index<CMap2::Vertex>(m,hd, phi1(m,hd2));
			if (is_indexed<CMap2::Edge>(m))
				copy_index<CMap2::Edge>(m,hd, hd2);
			if (is_indexed<CMap2::Volume>(m))
				copy_index<CMap2::Volume>(m,hd, hd2);
			return true;
		});
	}

	return hole;
}

uint32 close(CMap2& m,bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	foreach_dart(m,[&] (Dart d) -> bool
	{
		if (phi2(m,d) == d)
			fix_point_darts.push_back(d);
		return true;
	});

	for (Dart d : fix_point_darts)
	{
		if (phi2(m,d) == d)
		{
			CMap2::Face h = close_hole(m,d, set_indices);
			foreach_dart_of_orbit(m,h, [&] (Dart hd) -> bool { set_boundary(m,hd, true); return true; });
			++nb_holes;
		}
	}

	return nb_holes;
}

} // namespace cgogn
