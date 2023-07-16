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
#include <cgogn/core/types/maps/gmap/gmap2.h>

#include <cgogn/core/types/cell_marker.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

GMap2::Vertex cut_edge(GMap2& m, GMap2::Edge e, bool set_indices)
{
	Dart e1 = e.dart;
	Dart e2 = beta2(m, e1);
	Dart ne1 = cut_edge(static_cast<GMap1&>(m), GMap1::Edge{e1}, false).dart;
	Dart ne2 = cut_edge(static_cast<GMap1&>(m), GMap1::Edge{e2}, false).dart;
	beta2_sew(m, ne1, ne2);
	Dart ne1b = beta1(m, ne1);
	Dart ne2b = beta1(m, ne2);
	beta2_sew(m, ne1b, ne2b);

	bool b1 = is_boundary(m, e1);
	bool b2 = is_boundary(m, e2);
	set_boundary(m, ne1, b1);
	set_boundary(m, ne1b, b1);
	set_boundary(m, ne2, b2);
	set_boundary(m, ne2b, b2);

	GMap2::Vertex vert{ne1};
	if (set_indices)
	{
		if (is_indexed<GMap2::Vertex>(m))
			set_index(m, vert, new_index<GMap2::Vertex>(m));

		if (is_indexed<GMap2::Edge>(m))
		{
			uint32 ind1 = index_of(m, GMap2::Edge(e.dart));
			set_index<GMap2::Edge>(m, ne1, ind1);
			set_index<GMap2::Edge>(m, ne2, ind1);
			set_index(m, GMap2::Edge{ne1b}, new_index<GMap2::Edge>(m));
		}

		if (is_indexed<GMap2::Face>(m))
		{
			uint32 ind1 = index_of(m, GMap2::Face(e.dart));
			set_index<GMap2::Edge>(m, ne1, ind1);
			set_index<GMap2::Edge>(m, ne1b, ind1);
			uint32 ind2 = index_of(m, GMap2::Face(e2));
			set_index<GMap2::Edge>(m, ne2, ind2);
			set_index<GMap2::Edge>(m, ne2b, ind2);
		}
	}

	return vert;
}

GMap2::Vertex collapse_edge(GMap2& m, GMap2::Edge e, bool set_indices)
{
	Dart dd = e.dart;
	Dart dd_1 = beta1(m, dd);
	Dart dd_12 = beta2(m, dd_1);
	Dart ee = beta2(m, dd);
	Dart ee_1 = beta1(m, ee);
	Dart ee_12 = beta2(m, ee_1);

	collapse_edge(static_cast<GMap1&>(m), GMap1::Edge(dd), false);
	collapse_edge(static_cast<GMap1&>(m), GMap1::Edge(ee), false);

	uint32 codeg1 = codegree(m, GMap2::Face(dd_1));
	if (codeg1 == 2u)
	{
		Dart dd1 = beta1(m, dd_1);
		Dart dd12 = beta2(m, dd1);
		foreach_dart_of_orbit(m, GMap2::Face(dd_1), [&](Dart di) {
			beta2_unsew(m, di);
			return true;
		});
		beta2_sew(m, dd_12, dd12);
		beta2_sew(m, beta0(m, dd_12), beta0(m, dd12));
		remove_face(static_cast<GMap1&>(m), GMap1::Face(dd_1));
	}

	uint32 codeg2 = codegree(m, GMap2::Face(ee_1));
	if (codeg2 == 2u)
	{
		Dart ee1 = beta1(m, ee_1);
		Dart ee12 = beta2(m, ee1);
		foreach_dart_of_orbit(m, GMap2::Face(ee_1), [&](Dart di) {
			beta2_unsew(m, di);
			return true;
		});
		beta2_sew(m, ee_12, ee12);
		beta2_sew(m, beta0(m, ee_12), beta0(m, ee12));
		remove_face(static_cast<GMap1&>(m), GMap1::Face(ee_1));
	}

	GMap2::Vertex v(dd_12);

	if (set_indices)
	{
		if (is_indexed<GMap2::Vertex>(m))
			set_index(m, v, index_of(m, v));

		if (is_indexed<GMap2::Edge>(m))
		{
			if (codeg1 == 2u)
			{
				GMap2::Edge edg1(dd_12);
				set_index(m, edg1, index_of(m, edg1));
			}
			if (codeg2 == 2u)
			{
				GMap2::Edge edg2(ee_12);
				set_index(m, edg2, index_of(m, edg2));
			}
		}
	}

	return v;
}

bool flip_edge(GMap2& m, GMap2::Edge e, bool set_indices)
{
	auto externals = [&](Dart x) { return std::make_pair(beta1(m, x), beta1(m, beta2(m, x))); };

	auto turn = [&](const std::pair<Dart, Dart>& x) {
		return std::make_pair(beta<0, 1>(m, x.first), beta<1, 0>(m, x.second));
	};

	auto p1_darts = externals(e.dart);
	auto p2_darts = externals(beta<0, 2>(m, e.dart));
	auto pt1 = turn(p1_darts);
	auto pt2 = turn(p2_darts);

	beta1_unsew(m, p1_darts.first);
	beta1_unsew(m, p1_darts.second);
	beta1_unsew(m, p2_darts.first);
	beta1_unsew(m, p2_darts.second);
	// TODO remove unsew in release mode ?
	beta1_sew(m, p1_darts.first, p1_darts.second);
	beta1_sew(m, p2_darts.first, p2_darts.second);

	beta1_unsew(m, pt1.first);
	beta1_unsew(m, pt1.second);
	beta1_unsew(m, pt2.first);
	beta1_unsew(m, pt2.second);
	// TODO remove unsew in release mode ?

	Dart ed = e.dart;
	// is it more efficient like that ?
	// beta1_sew(m, pt1.first, e.dart);
	// beta1_sew(m, pt1.second, beta2(m,e.dart));
	// beta1_sew(m, pt2.first, beta<0,2>(m,e.dart));
	// beta1_sew(m, pt2.second, beta0(m,e.dart));
	beta1_sew(m, pt1.first, ed);
	ed = beta2(m, ed);
	beta1_sew(m, pt1.second, beta2(m, ed));
	ed = beta0(m, ed);
	beta1_sew(m, pt2.first, beta<0, 2>(m, ed));
	ed = beta0(m, ed);
	beta1_sew(m, pt2.second, beta0(m, ed));

	if (set_indices)
	{
		if (is_indexed<GMap2::Vertex>(m))
		{
			copy_index<GMap2::Vertex>(m, e.dart, pt1.first);
			copy_index<GMap2::Vertex>(m, ed, pt2.first);
		}

		if (is_indexed<GMap2::Face>(m))
		{
			copy_index<GMap2::Face>(m, p1_darts.first, pt1.first);
			copy_index<GMap2::Face>(m, beta0(m, p1_darts.first), beta0(m, pt1.first));
			copy_index<GMap2::Face>(m, p2_darts.first, pt2.first);
			copy_index<GMap2::Face>(m, beta0(m, p2_darts.first), beta0(m, pt2.first));
		}
	}

	return true;
}

bool edge_can_collapse(const GMap2& m, GMap2::Edge e)
{
	using Vertex = GMap2::Vertex;
	using Face = GMap2::Face;

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

bool edge_can_flip(const GMap2& m, GMap2::Edge e)
{
	if (is_incident_to_boundary(m, e))
		return false;

	Dart e1 = e.dart;
	Dart e2 = phi2(m, e1);

	auto next_edge = [&m](Dart d) { return phi<-1, 2>(m, d); };

	if (codegree(m, GMap2::Face(e1)) == 3 && codegree(m, GMap2::Face(e2)) == 3)
	{
		uint32 idxv2 = index_of(m, GMap2::Vertex(phi_1(m, e2)));
		Dart d = phi_1(m, e1);
		Dart it = d;
		do
		{
			if (index_of(m, GMap2::Vertex(phi1(m, it))) == idxv2)
				return false;
			it = next_edge(it);
		} while (it != d);
	}
	return true;
}

GMap2::Face add_face(GMap2& m, uint32 size, bool set_indices)
{
	GMap2::Face f = add_face(static_cast<GMap1&>(m), size, false);
	if (set_indices)
	{
		GMap2::Edge ep{f.dart};

		for (uint32 i = 0u; i < size; ++i)
		{
			GMap2::Edge e{phi1(m, ep.dart)};
			if (is_indexed<GMap2::Vertex>(m))
				set_index(m, GMap2::Vertex{e.dart}, new_index<GMap2::Vertex>(m));
			if (is_indexed<GMap2::Edge>(m))
				set_index(m, e, new_index<GMap2::Edge>(m));
			ep = e;
		}
		if (is_indexed<GMap2::Face>(m))
			set_index(m, GMap2::Face{ep.dart}, new_index<GMap2::Face>(m));
	}
	return f;
}

void merge_incident_faces(GMap2& m, GMap2::Edge e, bool set_indices)
{
	using Vertex = GMap2::Vertex;
	using Edge = GMap2::Edge;
	using Face = GMap2::Face;

	// DartMarker df = DartMarker(m);
	//
	// foreach_dart_of_BETA0_BETA1(m, e.dart, [&](Dart x)
	//{
	//	df.mark(x);
	//});

	// Dart d0 = beta2(m, e.dart);
	// do
	//{
	//	d0 = beta1(m, beta0(m, d0));
	//} while (df.is_marked(beta2(m,d0)));
	// d0 = beta0(m, beta1(m, d0));

	// Dart d1 = beta2(m, e.dart);
	// while (df.is_marked(beta2(m, d1)));
	//{
	//	d1 = beta0(m, beta1(m, d0));
	//}
	// d1 = beta1(m, beta0(m, d0));

	Dart d0 = e.dart;
	Dart d1 = beta<0, 2>(m, d0);

	Dart e1 = beta1(m, d0);
	beta1_unsew(m, e1);
	Dart e2 = beta<2, 1>(m, d0);
	beta1_unsew(m, e2);
	beta1_sew(m, e1, e2);

	Dart ee1 = beta1(m, d1);
	beta1_unsew(m, ee1);
	Dart ee2 = beta<2, 1>(m, d1);
	beta1_unsew(m, ee2);
	beta1_sew(m, ee1, ee2);

	remove_edge(m, GMap1::Edge{d0}, false);
	remove_edge(m, GMap1::Edge{d1}, false);

	if (set_indices)
	{
		if (is_indexed<Face>(m))
			set_index(m, Face{e1}, new_index<Face>(m));
	}
}

// WARNING: v1.dart & v2.dart must belong to the same orientation
GMap2::Edge cut_face(GMap2& m, GMap2::Vertex v1, GMap2::Vertex v2, bool set_indices)
{
	using Vertex = GMap2::Vertex;
	using Edge = GMap2::Edge;
	using Face = GMap2::Face;

	Dart d1 = v1.dart;
	Dart dd1 = beta1(m, d1);
	Dart d2 = v2.dart;
	Dart dd2 = beta1(m, d2);

	beta1_unsew(m, d1);
	beta1_unsew(m, d2);

	Dart e1 = add_edge(m, false).dart;
	Dart e2 = add_edge(m, false).dart;

	bool b1 = is_boundary(m, d1);
	set_boundary(m, e1, b1);
	set_boundary(m, beta0(m, e1), b1);
	set_boundary(m, e2, b1);
	set_boundary(m, beta0(m, e2), b1);

	beta1_sew(m, d1, e1);
	beta1_sew(m, d2, e2);
	beta1_sew(m, dd1, beta0(m, e2));
	beta1_sew(m, dd2, beta0(m, e1));

	beta2_sew(m, e1, beta0(m, e2));
	beta2_sew(m, e2, beta0(m, e1));

	if (set_indices)
	{
		if (is_indexed<Vertex>(m))
		{
			uint32 ind1 = index_of(m, v1);
			set_index<Vertex>(m, e1, ind1);
			set_index<Vertex>(m, beta0(m, e2), ind1);
			uint32 ind2 = index_of(m, v2);
			set_index<Vertex>(m, e2, ind2);
			set_index<Vertex>(m, beta0(m, e1), ind2);
		}
		if (is_indexed<Edge>(m))
			set_index(m, Edge{e1}, new_index<Edge>(m));

		if (is_indexed<Face>(m))
		{
			uint32 ind = index_of(m, Face(v1.dart));
			set_index<Face>(m, e1, ind);
			set_index<Face>(m, beta0(m, e2), ind);
			set_index<Face>(m, e2, ind);
			set_index<Face>(m, beta0(m, e1), ind);
		}
	}
	return Edge{e1};
}

GMap2::Volume add_pyramid(GMap2& m, uint32 size, bool set_indices)
{
	GMap1::Face first = add_face(static_cast<GMap1&>(m), 3u, false); // First triangle
																	 //	dump(first.dart);

	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next triangles
	{
		GMap1::Face next = add_face(static_cast<GMap1&>(m), 3u, false);
		phi2_sew(m, phi_1(m, current), phi1(m, next.dart));
		current = next.dart;
	}

	phi2_sew(m, phi_1(m, current), phi1(m, first.dart)); // Finish the umbrella

	GMap2::Face base = close_hole(m, first.dart, false); // Add the base face

	GMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<GMap2::Vertex>(m))
		{
			foreach_incident_vertex(
				m, vol,
				[&](GMap2::Vertex v) -> bool {
					set_index(m, v, new_index<GMap2::Vertex>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::HalfEdge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](GMap2::Edge e) -> bool {
					set_index(m, GMap2::HalfEdge(e.dart), new_index<GMap2::HalfEdge>(m));
					set_index(m, GMap2::HalfEdge(phi2(m, e.dart)), new_index<GMap2::HalfEdge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::Edge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](GMap2::Edge e) -> bool {
					set_index(m, e, new_index<GMap2::Edge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::Face>(m))
		{
			foreach_incident_face(
				m, vol,
				[&](GMap2::Face f) -> bool {
					set_index(m, f, new_index<GMap2::Face>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::Volume>(m))
			set_index(m, vol, new_index<GMap2::Volume>(m));
	}

	return vol;
}

GMap2::Volume add_prism(GMap2& m, uint32 size, bool set_indices)
{
	GMap1::Face first = add_face(static_cast<GMap1&>(m), 4u, false); // first quad
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next quads
	{
		GMap1::Face next = add_face(static_cast<GMap1&>(m), 4u, false);
		phi2_sew(m, phi_1(m, current), phi1(m, next.dart));
		current = next.dart;
	}
	phi2_sew(m, phi_1(m, current), phi1(m, first.dart)); // Finish the sides
	GMap2::Face base = close_hole(m, first.dart, false); // Add the base face
	close_hole(m, phi<1, 1>(m, first.dart), false);		 // Add the top face

	GMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<GMap2::Vertex>(m))
		{
			foreach_incident_vertex(
				m, vol,
				[&](GMap2::Vertex v) -> bool {
					set_index(m, v, new_index<GMap2::Vertex>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::HalfEdge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](GMap2::Edge e) -> bool {
					set_index(m, GMap2::HalfEdge(e.dart), new_index<GMap2::HalfEdge>(m));
					set_index(m, GMap2::HalfEdge(phi2(m, e.dart)), new_index<GMap2::HalfEdge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::Edge>(m))
		{
			foreach_incident_edge(
				m, vol,
				[&](GMap2::Edge e) -> bool {
					set_index(m, e, new_index<GMap2::Edge>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::Face>(m))
		{
			foreach_incident_face(
				m, vol,
				[&](GMap2::Face f) -> bool {
					set_index(m, f, new_index<GMap2::Face>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		if (is_indexed<GMap2::Volume>(m))
			set_index(m, vol, new_index<GMap2::Volume>(m));
	}

	return vol;
}

void remove_volume(GMap2& m, GMap2::Volume v)
{
	std::vector<Dart> darts;
	darts.reserve(96);
	foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
		darts.push_back(d);
		return true;
	});
	for (Dart d : darts)
		remove_dart(m, d);
}

GMap2::Face close_hole(GMap2& m, Dart d, bool set_indices)
{
	cgogn_message_assert(beta2(m, d) == d, "GMap2: close hole called on a dart that is not a phi2 fix point");

	std::vector<Dart> vd;
	vd.reserve(128);
	Dart dd1 = d;
	do
	{
		Dart d_next;
		dd1 = beta0(m, dd1);
		do
		{
			d_next = beta1(m, dd1);
			dd1 = beta2(m, d_next);
		} while (d_next != dd1);
		vd.push_back(dd1);
	} while (dd1 != d);

	GMap1::Face f = add_face(static_cast<GMap1&>(m), vd.size(), false);
	Dart df = f.dart;
	for (Dart dh : vd)
	{
		phi2_sew(m, dh, df);
		df = phi_1(m, df);
	}

	if (set_indices)
	{
		foreach_dart_of_orbit(m, f, [&](Dart fd) -> bool {
			Dart hd = phi2(m, fd);
			if (is_indexed<GMap2::Vertex>(m))
				copy_index<GMap2::Vertex>(m, fd, hd);
			if (is_indexed<GMap2::Edge>(m))
				copy_index<GMap2::Edge>(m, fd, hd);
			if (is_indexed<GMap2::Volume>(m))
				copy_index<GMap2::Volume>(m, fd, hd);
			return true;
		});
	}

	return f;
}

int32 close(GMap2& m, bool set_indices)
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
			GMap2::Face h = close_hole(m, d, set_indices);
			foreach_dart_of_orbit(m, h, [&](Dart hd) -> bool {
				set_boundary(m, hd, true);
				return true;
			});
			++nb_holes;
		}
	}

	return nb_holes;
}

/*************************************************************************/
// Debugging helper functions
/*************************************************************************/

bool check_integrity(GMap2& m, bool verbose)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = true;
		relations &= beta<0, 2, 0, 2>(m, d) == d;
		relations &= beta1(m, d) != d;
		if (verbose && !relations)
		{
			std::cerr << "Dart " << d << " has bad relations" << std::endl;
			if (beta1(m, d) == d)
				std::cerr << "beta1 fixed point" << std::endl;
			if (beta<0, 2, 0, 2>(m, d) != d)
				std::cerr << "edge constraint not respected" << std::endl;
		}
	}
	result &= check_indexing<GMap2::Vertex>(m);
	result &= check_indexing<GMap2::HalfEdge>(m);
	result &= check_indexing<GMap2::Edge>(m);
	result &= check_indexing<GMap2::Face>(m);
	result &= check_indexing<GMap2::Volume>(m);
	return result;
}

} // namespace cgogn
