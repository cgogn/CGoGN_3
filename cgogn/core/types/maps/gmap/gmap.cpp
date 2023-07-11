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
#include <iomanip>
#include <cgogn/core/types/map/gmap/gmap3.h>
#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{

void dump_map_darts(const GMapBase& m)
{
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		std::cout << "index: " << std::setw(5) << d.index << " / ";
		for (auto& r : m.relations_)
			std::cout << r->name() << ": " << std::setw(5) << (*r)[d.index] << " / ";
		//		for (auto& ind : m.cells_indices_)
		//			if (ind)
		//				std::cout << ind->name() << ": " << std::setw(5) << (*ind)[d.index] << " / ";
		for (uint32 orb : m.cells_embedded_orbit_)
		{
			auto& ind = m.cells_indices_[orb];
			std::cout << ind->name() << ": " << std::setw(5) << (*ind)[d.index] << " / ";
		}
		std::cout << std::endl;
	}
}

bool check_integrity(GMap1& m, bool verbose)
{
	bool result = true;
	//for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	//{
	//	//		if (verbose && !relations)
	//	//			std::cerr << "Dart " << d << " has bad relations" << std::endl;

	//	//		result &= relations;
	//}
	result &= check_indexing<GMap1::Vertex>(m);
	result &= check_indexing<GMap1::Edge>(m);
	result &= check_indexing<GMap1::Face>(m);
	return result;
}


GMap0::Edge add_edge(GMap0& m, bool set_indices)
{
	Dart d = add_dart(m);
	Dart e = add_dart(m);
	beta0_sew(m, d, e);

	GMap0::Edge edge{d};

	if (set_indices)
	{
		if (is_indexed<GMap0::Vertex>(m))
		{
			set_index(m, GMap0::Vertex{d}, new_index<GMap0::Vertex>(m));
			set_index(m, GMap0::Vertex{e}, new_index<GMap0::Vertex>(m));
		}

		if (is_indexed<GMap0::Edge>(m))
		{
			set_index(m, edge, new_index<GMap0::Edge>(m));
		}
	}
	return edge;
}

void remove_edge(GMap0& m, GMap0::Edge e, bool set_indice)
{
	remove_dart(m, e.dart);
	remove_dart(m, beta0(m, e.dart));
}

GMap1::Face add_face(GMap1& m, uint32 size, bool set_indices)
{
	using Vertex = GMap1::Vertex;
	using Edge = GMap1::Edge;
	using Face = GMap1::Face;

	Edge e0 = add_edge(m, false);
	Edge ep = e0;
	for (uint32 i = 1u; i < size; ++i)
	{
		Edge e = add_edge(m, false);
		beta1_sew(m, e.dart, beta0(m, ep.dart));
		ep = e;
	}
	beta1_sew(m, e0.dart, beta0(m, ep.dart));


	if (set_indices)
	{
		for (uint32 i = 0u; i < size; ++i)
		{
			Edge e{phi1(m, ep.dart)};
			if (is_indexed<Vertex>(m))
				set_index(m, Vertex{e.dart}, new_index<Vertex>(m));
			if (is_indexed<Edge>(m))
				set_index(m, e, new_index<Edge>(m));

			ep = e;
		}
		if (is_indexed<Face>(m))
			set_index(m, Face{ep.dart}, new_index<Face>(m));
	}
	return Face{e0.dart};
}

void remove_face(GMap1& m, GMap1::Face f)
{
	Dart it = phi1(m, f.dart);
	while (it != f.dart)
	{
		Dart next = phi1(m, it);
		remove_edge(m, GMap1::Edge{it});
		it = next;
	}
	remove_dart(m, f.dart);
}


GMap1::Vertex cut_edge(GMap1& m, GMap1::Edge e, bool set_indices)
{
	using Vertex = GMap1::Vertex;
	using Edge = GMap1::Edge;
	using Face = GMap1::Face;

	Dart e0 = e.dart;
	Dart e1 = beta0(m, e0);

	Dart vd0 = add_dart(m);
	Dart vd1 = add_dart(m);

	beta0_sew(m, e0, vd0);
	beta0_sew(m, e1, vd1);
	beta1_sew(m, vd0, vd1);

	Vertex v(vd0);

	if (set_indices)
	{
		if (is_indexed<Vertex>(m))
			set_index(m, v, new_index<Vertex>(m));
		if (is_indexed<Edge>(m))
		{
			copy_index<Edge>(m, vd0, e.dart);
			copy_index<Edge>(m, vd1, e.dart);
		}

		if (is_indexed<Face>(m))
		{
			copy_index<Face>(m, vd0, e.dart);
			copy_index<Face>(m, vd1, e.dart);
		}
	}

	return v;
}


GMap1::Vertex collapse_edge(GMap1& m, GMap1::Edge e, bool set_indices)
{
	using Vertex = GMap1::Vertex;
	Dart d1 = beta1(m, e.dart);
	Dart d2 = beta1(m, beta0(m, e.dart));

	beta1_unsew(m, d1);
	beta1_unsew(m, d2);

	beta1_sew(m, d1, d2);

	remove_edge(m, e);

	Vertex v(d1);

	if (set_indices)
		copy_index<Vertex>(m, d2, d1);

	return v;
}





bool check_integrity(GMap2& m, bool verbose)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = true;
		relations &= beta<0,2,0,2>(m, d) == d;
		relations &= beta1(m,d) != d;
		if (verbose && !relations)
		{
			std::cerr << "Dart " << d << " has bad relations" << std::endl;
			if (beta1(m, d)==d)
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


GMap2::Vertex cut_edge(GMap2& m, GMap2::Edge e, bool set_indices)
{
	Dart e2 = beta2(m, e.dart);
	Dart d1 = (cut_edge(static_cast<GMap1&>(m), GMap1::Edge{e.dart}, false)).dart;
	Dart d2 = (cut_edge(static_cast<GMap1&>(m), GMap1::Edge{e2}, false)).dart;
	beta2_sew(m, d1, d2);
	beta2_sew(m, beta1(m, d1), beta1(m, d2));

	GMap2::Vertex vert{d1};
	if (set_indices)
	{
		if (is_indexed<GMap2::Vertex>(m))
			set_index(m, vert, new_index<GMap2::Vertex>(m));

		if (is_indexed<GMap2::Edge>(m))
		{
			uint32 ind = index_of(m, GMap2::Edge(e.dart));
			set_index<GMap2::Edge>(m, d1, ind);
			set_index<GMap2::Edge>(m, d2, ind);
			set_index<GMap2::Edge>(m, beta1(m, d1), ind);
			set_index<GMap2::Edge>(m, beta1(m, d2), ind);
		}

		if (is_indexed<GMap2::Face>(m))
		{
			uint32 ind1 = index_of(m, GMap2::Face(e.dart));
			set_index<GMap2::Edge>(m, d1, ind1);
			set_index<GMap2::Edge>(m, beta1(m, d1), ind1);
			uint32 ind2 = index_of(m, GMap2::Face(e2));
			set_index<GMap2::Edge>(m, d2, ind2);
			set_index<GMap2::Edge>(m, beta1(m, d2), ind2);
		}
	}
	return vert;
}

// WARNING v1.dart & v2.dart must belong to the same  orientation
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




void CGOGN_CORE_EXPORT merge_incident_faces(GMap2& m, GMap2::Edge e, bool set_indices)
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

	remove_edge(m, GMap1::Edge{d0},false);
	remove_edge(m, GMap1::Edge{d1}, false);

	if (set_indices)
	{
		if (is_indexed<Face>(m))
			set_index(m, Face{e1}, new_index<Face>(m));
	}

}


GMap2::Face close_hole(GMap2& m, Dart d, bool set_indices)
{
	cgogn_message_assert(phi2(m, d) == d, "CMap2: close hole called on a dart that is not a phi2 fix point");

	std::vector<Dart> vd;
	vd.reserve(128);
	Dart d_next = d; // Turn around the hole
	Dart dd1; // to complete the face
	do
	{
		d_next = beta0(m, d_next);
		do
		{
			dd1 = beta1(m, d_next);
			d_next = beta2(m, dd1);
		} while (d_next != dd1 && dd1 != d);

		if (dd1 != d)
			vd.push_back(d);
	} while (dd1 != d);

	GMap2::Face f = add_face(m, vd.size(), false);
	Dart df = f.dart;
	for (Dart dh : vd)
	{
		phi2_sew(m, dh, df);
		df = phi1(m, df);
	}

	if (set_indices)
	{
		foreach_dart_of_orbit(m,f , [&](Dart fd) -> bool {
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
//				set_boundary(m, hd, true);
				return true;
			});
			++nb_holes;
		}
	}

	return nb_holes;
}


// SAME CODE AS IN CMAP2
bool CGOGN_CORE_EXPORT edge_can_flip(const GMap2& m, GMap2::Edge e)
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


bool CGOGN_CORE_EXPORT flip_edge(GMap2& m, GMap2::Edge e, bool set_indices)
{
	auto externals = [&](Dart x) { return std::make_pair(beta1(m, x), beta1(m, beta2(m, x))); };

	auto turn = [&](const std::pair<Dart,Dart> & x) { return std::make_pair(beta<0,1>(m, x.first), beta<1,0>(m, x.second)); };

	auto p1_darts = externals(e.dart);
	auto p2_darts = externals(beta<0,2>(m, e.dart));
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
	//beta1_sew(m, pt1.first, e.dart);
	//beta1_sew(m, pt1.second, beta2(m,e.dart));
	//beta1_sew(m, pt2.first, beta<0,2>(m,e.dart));
	//beta1_sew(m, pt2.second, beta0(m,e.dart));
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


GMap2::Face CGOGN_CORE_EXPORT add_face(GMap2& m, uint32 size, bool set_indices)
{
	// NOTHING TO DO, THANKS TO FIXED POINTS BOUNDARY
	return add_face(static_cast<GMap1&>(m), size, set_indices);
}

// SAME COE AS IN CMAP2
void CGOGN_CORE_EXPORT remove_volume(GMap2& m, GMap2::Volume v)
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


//void CGOGN_CORE_EXPORT reverse_orientation(GMap2& m);

// SAME CODE  AS IN CMAP2
bool CGOGN_CORE_EXPORT edge_can_collapse(const GMap2& m, GMap2::Edge e)
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


//WIP
//GMap2::Vertex CGOGN_CORE_EXPORT collapse_edge(GMap2& m, GMap2::Edge e, bool set_indices)
//{
//	Dart dd = e.dart;
//	Dart dd_1 = phi_1(m, dd);
//	Dart dd_12 = phi2(m, dd_1);

	//collapse_edge(static_cast<GMap1&>(m), GMap1::Edge(dd), false);


	//if (codegree(m, GMap2::Face(dd_1)) == 2u)
	//{
	//	Dart dd1 = phi1(m, dd_1);
	//	Dart dd12 = phi2(m, dd1);
	//	phi2_unsew(m, dd1);
	//	phi2_unsew(m, dd_1);
	//	phi2_sew(m, dd12, dd_12);
	//	remove_face(static_cast<CMap1&>(m), CMap1::Face(dd1), false);
	//}

	// Dart ee = phi2(m, dd);
	// Dart ee_1 = phi_1(m, ee);
	// Dart ee_12 = phi2(m, ee_1);
	//collapse_edge(static_cast<GMap1&>(m), GMap1::Edge(ee), false);
	//if (codegree(m, GMap2::Face(ee_1)) == 2u)
	//{
	//	Dart ee1 = phi1(m, ee_1);
	//	Dart ee12 = phi2(m, ee1);
	//	phi2_unsew(m, ee1);
	//	phi2_unsew(m, ee_1);
	//	phi2_sew(m, ee12, ee_12);
	//	remove_face(static_cast<CMap1&>(m), CMap1::Face(ee1), false);
	//}

	//CMap2::Vertex v(dd_12);

	//if (set_indices)
	//{
	//	if (is_indexed<CMap2::Vertex>(m))
	//		set_index(m, v, index_of(m, v));
	//	if (is_indexed<CMap2::Edge>(m))
	//	{
	//		copy_index<CMap2::Edge>(m, dd_12, phi2(m, dd_12));
	//		copy_index<CMap2::Edge>(m, ee_12, phi2(m, ee_12));
	//	}
	//}

	//return v;
//}


} // namespace cgogn
                    