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

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

CMap2::Vertex cut_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart d1 = e.dart;
	Dart d2 = phi2(m, d1);
	phi2_unsew(m, d1);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d1), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(d2), false);
	phi2_sew(m, d1, nv2.dart);
	phi2_sew(m, d2, nv1.dart);
	set_boundary(m, nv1.dart, is_boundary(m, d1));
	set_boundary(m, nv2.dart, is_boundary(m, d2));
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
			copy_index<CMap2::Edge>(m, nv2.dart, d1);
			set_index(m, CMap2::Edge(nv1.dart), new_index<CMap2::Edge>(m));
		}
		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m, nv1.dart, d1);
			copy_index<CMap2::Face>(m, nv2.dart, d2);
		}
		if (is_indexed<CMap2::Volume>(m))
		{
			copy_index<CMap2::Volume>(m, nv1.dart, d1);
			copy_index<CMap2::Volume>(m, nv2.dart, d2);
		}
	}

	return v;
}

CMap2::Vertex collapse_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart dd = e.dart;
	Dart dd_1 = phi_1(m, dd);
	Dart dd_12 = phi2(m, dd_1);
	Dart ee = phi2(m, dd);
	Dart ee_1 = phi_1(m, ee);
	Dart ee_12 = phi2(m, ee_1);

	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	collapse_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);

	if (codegree(m, CMap2::Face(dd_1)) == 2u)
	{
		Dart dd1 = phi1(m, dd_1);
		Dart dd12 = phi2(m, dd1);
		phi2_unsew(m, dd1);
		phi2_unsew(m, dd_1);
		phi2_sew(m, dd12, dd_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(dd1));
	}

	if (codegree(m, CMap2::Face(ee_1)) == 2u)
	{
		Dart ee1 = phi1(m, ee_1);
		Dart ee12 = phi2(m, ee1);
		phi2_unsew(m, ee1);
		phi2_unsew(m, ee_1);
		phi2_sew(m, ee12, ee_12);
		remove_face(static_cast<CMap1&>(m), CMap1::Face(ee1));
	}

	CMap2::Vertex v(dd_12);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
			set_index(m, v, index_of(m, v));
		if (is_indexed<CMap2::Edge>(m))
		{
			copy_index<CMap2::Edge>(m, dd_12, phi2(m, dd_12));
			copy_index<CMap2::Edge>(m, ee_12, phi2(m, ee_12));
		}
	}

	return v;
}

bool flip_edge(CMap2& m, CMap2::Edge e, bool set_indices)
{
	Dart d = e.dart;
	Dart dd = phi2(m, d);
	Dart d1 = phi1(m, d);
	Dart d_1 = phi_1(m, d);
	Dart dd1 = phi1(m, dd);
	Dart dd_1 = phi_1(m, dd);

	// // Cannot flip edge whose incident faces have co-degree 1
	// if (d == d1 || dd == dd1)
	// 	return false;

	// // Both vertices have degree 1 and thus nothing is done // TODO may return true ?
	// if (d == dd_1 && dd == d_1)
	// 	return false;

	// if (d != dd_1)
	// 	phi1_sew(m, d, dd_1); // Detach the edge from its
	// if (dd != d_1)
	// 	phi1_sew(m, dd, d_1); // two incident vertices

	// if (d != dd_1)
	// 	phi1_sew(m, d, d1); // Insert the first end in its new vertices
	// if (dd != d_1)
	// 	phi1_sew(m, dd, dd1); // Insert the second end in its new vertices

	phi1_sew(m, d, dd_1);
	phi1_sew(m, dd, d_1);
	phi1_sew(m, d, d1);
	phi1_sew(m, dd, dd1);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			copy_index<CMap2::Vertex>(m, d, phi1(m, dd));
			copy_index<CMap2::Vertex>(m, dd, phi1(m, d));
		}

		if (is_indexed<CMap2::Face>(m))
		{
			copy_index<CMap2::Face>(m, phi_1(m, d), d);
			copy_index<CMap2::Face>(m, phi_1(m, dd), dd);
		}
	}

	return true;
}

bool edge_can_collapse(const CMap2& m, CMap2::Edge e)
{
	using Vertex = CMap2::Vertex;
	using Face = CMap2::Face;

	auto vertices = incident_vertices(m, e);

	if (is_incident_to_boundary(m, vertices[0]) || is_incident_to_boundary(m, vertices[1]))
		return false;

	uint32 val_v1 = degree(m, vertices[0]);
	uint32 val_v2 = degree(m, vertices[1]);

	if (val_v1 + val_v2 < 8 || val_v1 + val_v2 > 14)
		return false;

	Dart e1 = e.dart;
	Dart e2 = phi2(m, e.dart);
	if (codegree(m, Face(e1)) == 3)
	{
		if (degree(m, Vertex(phi_1(m, e1))) < 4)
			return false;
	}
	if (codegree(m, Face(e2)) == 3)
	{
		if (degree(m, Vertex(phi_1(m, e2))) < 4)
			return false;
	}

	auto next_edge = [&m](Dart d) { return phi<-1, 2>(m, d); };

	// Check vertex sharing condition
	std::vector<uint32> vn1;
	Dart it = next_edge(next_edge(e1));
	Dart end = phi1(m, e2);
	do
	{
		vn1.push_back(index_of(m, Vertex(phi1(m, it))));
		it = next_edge(it);
	} while (it != end);
	it = next_edge(next_edge(e2));
	end = phi1(m, e1);
	do
	{
		auto vn1it = std::find(vn1.begin(), vn1.end(), index_of(m, Vertex(phi1(m, it))));
		if (vn1it != vn1.end())
			return false;
		it = next_edge(it);
	} while (it != end);

	return true;
}

bool edge_can_flip(const CMap2& m, CMap2::Edge e)
{
	if (is_incident_to_boundary(m, e))
		return false;

	Dart e1 = e.dart;
	Dart e2 = phi2(m, e1);

	auto next_edge = [&m](Dart d) { return phi<-1, 2>(m, d); };

	if (codegree(m, CMap2::Face(e1)) == 3 && codegree(m, CMap2::Face(e2)) == 3)
	{
		uint32 idxv2 = index_of(m, CMap2::Vertex(phi_1(m, e2)));
		Dart d = phi_1(m, e1);
		Dart it = d;
		do
		{
			if (index_of(m, CMap2::Vertex(phi1(m, it))) == idxv2)
				return false;
			it = next_edge(it);
		} while (it != d);
	}

	return true;
}

CMap2::Face add_face(CMap2& m, uint32 size, bool set_indices)
{
	CMap2::Face f = add_face(static_cast<CMap1&>(m), size, false);
	CMap2::Face b = add_face(static_cast<CMap1&>(m), size, false);
	Dart it = b.dart;
	foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
		set_boundary(m, it, true);
		phi2_sew(m, d, it);
		it = phi_1(m, it);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(
				m, f,
				[&](CMap2::Vertex v) -> bool {
					set_index(m, v, new_index<CMap2::Vertex>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(
				m, f,
				[&](CMap2::Edge e) -> bool {
					set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(
				m, f,
				[&](CMap2::Edge e) -> bool {
					set_index(m, e, new_index<CMap2::Edge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Face>(m))
			set_index(m, f, new_index<CMap2::Face>(m));
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, CMap2::Volume(f.dart), new_index<CMap2::Volume>(m));
	}

	return f;
}

void merge_incident_faces(CMap2& m, CMap2::Edge e, bool set_indices)
{
	if (is_incident_to_boundary(m, e))
		return;

	Dart d = e.dart;
	Dart d_1 = phi_1(m, d);
	Dart d2 = phi2(m, d);
	Dart d2_1 = phi_1(m, d2);

	phi1_sew(m, d_1, d2);
	phi1_sew(m, d2_1, d);

	remove_face(static_cast<CMap1&>(m), CMap1::Face(d));

	if (set_indices)
	{
		if (is_indexed<CMap2::Face>(m))
			set_index(m, CMap2::Face(d_1), index_of(m, CMap2::Face(d_1)));
	}
}

CMap2::Edge cut_face(CMap2& m, CMap2::Vertex v1, CMap2::Vertex v2, bool set_indices)
{
	Dart dd = phi_1(m, v1.dart);
	Dart ee = phi_1(m, v2.dart);
	CMap1::Vertex nv1 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(dd), false);
	CMap1::Vertex nv2 = cut_edge(static_cast<CMap1&>(m), CMap1::Edge(ee), false);
	phi1_sew(m, nv1.dart, nv2.dart);
	phi2_sew(m, nv1.dart, nv2.dart);
	set_boundary(m, nv1.dart, is_boundary(m, dd));
	set_boundary(m, nv2.dart, is_boundary(m, ee));
	CMap2::Edge e(nv1.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			copy_index<CMap2::Vertex>(m, nv1.dart, v1.dart);
			copy_index<CMap2::Vertex>(m, nv2.dart, v2.dart);
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
			copy_index<CMap2::Face>(m, nv2.dart, v1.dart);
			set_index(m, CMap2::Face(v2.dart), new_index<CMap2::Face>(m));
		}
		if (is_indexed<CMap2::Volume>(m))
		{
			copy_index<CMap2::Volume>(m, nv1.dart, v2.dart);
			copy_index<CMap2::Volume>(m, nv2.dart, v1.dart);
		}
	}

	return e;
}

CMap2::Volume add_pyramid(CMap2& m, uint32 size, bool set_indices)
{
	CMap1::Face first = add_face(static_cast<CMap1&>(m), 3u, false); // First triangle
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next triangles
	{
		CMap1::Face next = add_face(static_cast<CMap1&>(m), 3u, false);
		phi2_sew(m, phi_1(m, current), phi1(m, next.dart));
		current = next.dart;
	}
	phi2_sew(m, phi_1(m, current), phi1(m, first.dart)); // Finish the umbrella
	CMap2::Face base = close_hole(m, first.dart, false); // Add the base face

	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(
				m, vol,
				[&](CMap2::Vertex v) -> bool {
					set_index(m, v, new_index<CMap2::Vertex>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](CMap2::Edge e) -> bool {
					set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
					set_index(m, CMap2::HalfEdge(phi2(m, e.dart)), new_index<CMap2::HalfEdge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](CMap2::Edge e) -> bool {
					set_index(m, e, new_index<CMap2::Edge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Face>(m))
		{
			foreach_incident_face(
				m, vol,
				[&](CMap2::Face f) -> bool {
					set_index(m, f, new_index<CMap2::Face>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

CMap2::Volume add_prism(CMap2& m, uint32 size, bool set_indices)
{
	CMap1::Face first = add_face(static_cast<CMap1&>(m), 4u, false); // first quad
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next quads
	{
		CMap1::Face next = add_face(static_cast<CMap1&>(m), 4u, false);
		phi2_sew(m, phi_1(m, current), phi1(m, next.dart));
		current = next.dart;
	}
	phi2_sew(m, phi_1(m, current), phi1(m, first.dart)); // Finish the sides
	CMap2::Face base = close_hole(m, first.dart, false); // Add the base face
	close_hole(m, phi<1, 1>(m, first.dart), false);		 // Add the top face

	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(
				m, vol,
				[&](CMap2::Vertex v) -> bool {
					set_index(m, v, new_index<CMap2::Vertex>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](CMap2::Edge e) -> bool {
					set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
					set_index(m, CMap2::HalfEdge(phi2(m, e.dart)), new_index<CMap2::HalfEdge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](CMap2::Edge e) -> bool {
					set_index(m, e, new_index<CMap2::Edge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Face>(m))
		{
			foreach_incident_face(
				m, vol,
				[&](CMap2::Face f) -> bool {
					set_index(m, f, new_index<CMap2::Face>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

void remove_volume(CMap2& m, CMap2::Volume v)
{
	std::vector<Dart> darts;
	foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
		darts.push_back(d);
		return true;
	});
	for (Dart d : darts)
		remove_dart(m, d);
}

CMap2::Face close_hole(CMap2& m, Dart d, bool set_indices)
{
	cgogn_message_assert(phi2(m, d) == d, "CMap2: close hole called on a dart that is not a phi2 fix point");

	Dart first = add_dart(m); // First edge of the face that will fill the hole
	phi2_sew(m, d, first);	  // 2-sew the new edge to the hole

	Dart d_next = d; // Turn around the hole
	Dart d_phi1;	 // to complete the face
	do
	{
		do
		{
			d_phi1 = phi1(m, d_next); // Search and put in d_next
			d_next = phi2(m, d_phi1); // the next dart of the hole
		} while (d_next != d_phi1 && d_phi1 != d);

		if (d_phi1 != d)
		{
			Dart next = add_dart(m); // Add a vertex into the built face
			phi1_sew(m, first, next);
			phi2_sew(m, d_next, next); // and 2-sew the face to the hole
		}
	} while (d_phi1 != d);

	CMap2::Face hole(first);

	if (set_indices)
	{
		foreach_dart_of_orbit(m, hole, [&](Dart hd) -> bool {
			Dart hd2 = phi2(m, hd);
			if (is_indexed<CMap2::Vertex>(m))
				copy_index<CMap2::Vertex>(m, hd, phi1(m, hd2));
			if (is_indexed<CMap2::Edge>(m))
				copy_index<CMap2::Edge>(m, hd, hd2);
			if (is_indexed<CMap2::Volume>(m))
				copy_index<CMap2::Volume>(m, hd, hd2);
			return true;
		});
	}

	return hole;
}

uint32 close(CMap2& m, bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		if (phi2(m, d) == d)
			fix_point_darts.push_back(d);

	for (Dart d : fix_point_darts)
	{
		if (phi2(m, d) == d)
		{
			CMap2::Face h = close_hole(m, d, set_indices);
			foreach_dart_of_orbit(m, h, [&](Dart hd) -> bool {
				set_boundary(m, hd, true);
				return true;
			});
			++nb_holes;
		}
	}

	return nb_holes;
}

void reverse_orientation(CMap2& m)
{
	if (is_indexed<CMap2::Vertex>(m))
	{
		auto new_vertex_indices = m.darts_.add_attribute<uint32>("__new_vertex_indices");
		new_vertex_indices->fill(INVALID_INDEX);
		for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
			(*new_vertex_indices)[d.index] = index_of(m, CMap2::Vertex(phi1(m, d)));
		m.cells_indices_[CMap2::Vertex::ORBIT]->swap(new_vertex_indices.get());
		m.darts_.remove_attribute(new_vertex_indices);
	}

	m.phi1_->swap(m.phi_1_.get());
}

/*************************************************************************/
// Debugging helper functions
/*************************************************************************/

bool check_integrity(CMap2& m, bool verbose)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = true;
		relations &= phi2(m, d) != d && phi<2, 2>(m, d) == d;
		relations &= phi<-1, 1>(m, d) == d && phi<1, -1>(m, d) == d;
		if (verbose && !relations)
		{
			std::cerr << "Dart " << d << " has bad relations" << std::endl;
			if (phi2(m, d) == d)
				std::cerr << "  phi2 fixed point" << std::endl;
			if (phi<2, 2>(m, d) != d)
				std::cerr << "  phi2 not involution" << std::endl;
		}

		bool boundary =
			is_boundary(m, d) == is_boundary(m, phi1(m, d)) && (!is_boundary(m, d) || !is_boundary(m, phi2(m, d)));
		if (verbose && !boundary)
			std::cerr << "Dart " << d << " has bad boundary" << std::endl;

		result &= relations && boundary;
	}
	result &= check_indexing<CMap2::Vertex>(m);
	result &= check_indexing<CMap2::HalfEdge>(m);
	result &= check_indexing<CMap2::Edge>(m);
	result &= check_indexing<CMap2::Face>(m);
	result &= check_indexing<CMap2::Volume>(m);
	return result;
}

} // namespace cgogn
