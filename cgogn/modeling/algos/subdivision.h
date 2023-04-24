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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_H_

#include <cgogn/modeling/algos/subdivision/basic.h>
#include <cgogn/modeling/algos/subdivision_utils.h>
#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/cmap/cph3.h>

#include <cgogn/geometry/algos/angle.h>

namespace cgogn
{

namespace modeling
{

using geometry::Vec3;

///////////
// CMap3 //
///////////

/// Pre cut hex, input: original dart from edge along the direction of the cut
inline void quadrisect_hex(CMap3& m, CMap3::Volume w)
{
	Dart d0 = phi<1, 2, 1, 1, 1>(m, w.dart);
	std::vector<Dart> path0;
	Dart d1 = d0;
	do
	{
		path0.push_back(d1);
		d1 = phi<1, 2, 1>(m, d1);
	} while (d1 != d0);

	cut_volume(m, path0);
	path0.clear();

	d1 = phi<2, 1>(m, d0);
	Dart d2 = phi<1, 1, 1>(m, d1);
	CMap3::Edge e = cut_face(m, CMap3::Vertex(d1), CMap3::Vertex(d2));

	std::vector<Dart> path1;
	d1 = e.dart;
	d2 = phi3(m, d1);
	do
	{
		path0.push_back(d1);
		path1.push_back(d2);
		d1 = phi<1, 2, 1>(m, d1);
		d2 = phi<1, 2, 1>(m, d2);
	} while (d1 != e.dart);
	cut_volume(m, path0);
	cut_volume(m, path1);

	// Dart d0 = phi<1, 2, 1, 1>(m, w.dart);
	// Dart d1 = d0;
	// Dart d2, d3, d4;
	// do {
	// 	d1 = phi1(m, d1);
	// 	d2 = phi2(m, d1);

	// 	Dart f0 = add_face(static_cast<CMap1&>(m), 4, false).dart;
	// 	Dart f1 = add_face(static_cast<CMap1&>(m), 4, false).dart;
	// 	Dart ee = f0;
	// 	Dart ff = f1;
	// 	do
	// 	{
	// 		phi3_sew(m, ee, ff);
	// 		ee = phi1(m, ee);
	// 		ff = phi_1(m, ff);
	// 	} while (ee != f0);

	// 	phi2_unsew(m, d1);
	// 	phi2_sew(m, d1, f0);
	// 	phi2_sew(m, d2, f1);
	// 	f0 = phi1(m, f0);
	// 	f1 = phi_1(m, f1);

	// 	d3 = phi<1, 2, 1>(m, d1);
	// 	d4 = phi2(m, d3);

	// 	phi2_unsew(m, d3);
	// 	phi2_sew(m, d3, f0);
	// 	phi2_sew(m, d4, f1);
	// 	f0 = phi1(m, f0);
	// 	f1 = phi_1(m, f1);

	// 	d3 = phi<1, 2, 1>(m, d3);
	// 	d4 = phi2(m, d3);

	// 	phi2_unsew(m, d3);
	// 	phi2_sew(m, d3, f0);
	// 	phi2_sew(m, d4, f1);

	// 	d1 = d2;
	// } while(d1 != d0);

	// do{
	// 	d2 = phi1(m, d1);

	// 	d3 = phi<2, -1>(m, d1);
	// 	d4 = phi<2, 1>(m, d2);
	// 	phi2_sew(m, d3, d4);

	// 	d1 = phi<2, 3, 2>(m, d2);
	// }while(d1 != d0);

	// d0 = phi<1, 2, 1, 1, 2, -1>(m, w.dart);
	// d1 = d0;
	// d2 = phi2(m, d1);

	// if (is_indexed<CMap3::Vertex>(m))
	// {
	// 	do{
	// 		set_index<CMap3::Vertex>(m, d1, index_of(m, CMap3::Vertex(phi1(m, d2))));
	// 		set_index<CMap3::Vertex>(m, d2, index_of(m, CMap3::Vertex(phi1(m, d1))));
	// 		d3 = phi1(m, d1);
	// 		do{
	// 			d4 = phi3(m, d3);
	// 			set_index<CMap3::Vertex>(m, d3, index_of(m, CMap3::Vertex(phi<2, 1>(m, d3))));
	// 			set_index<CMap3::Vertex>(m, d4, index_of(m, CMap3::Vertex(phi<2, 1>(m, d4))));
	// 			d3 = phi1(m, d3);
	// 		} while(d3 != d1);
	// 		d1 = phi3(m, d2);
	// 		d2 = phi2(m, d1);
	// 	} while(d1 != d0);
	// }

	// if (is_indexed<CMap3::Vertex2>(m))
	// {
	// 	do {
	// 		set_index<CMap3::Vertex2>(m, d1, index_of(m, CMap3::Vertex2(phi1(m, d2))));
	// 		set_index<CMap3::Vertex2>(m, d2, index_of(m, CMap3::Vertex2(phi1(m, d1))));
	// 		d3 = phi1(m, d1);
	// 		do {
	// 			d4 = phi3(m, d3);
	// 			set_index<CMap3::Vertex2>(m, d3, index_of(m, CMap3::Vertex2(phi<2, 1>(m, d3))));
	// 			set_index<CMap3::Vertex2>(m, d4, index_of(m, CMap3::Vertex2(phi<2, 1>(m, d4))));
	// 			d3 = phi1(m, d3);
	// 		} while(d3 != d1);
	// 		d1 = phi3(m, d2);
	// 		d2 = phi2(m, d1);
	// 	} while(d1 != d0);
	// }

	// if(is_indexed<CMap3::Edge>(m))
	// {
	// 	set_index<CMap3::Edge>(m, CMap3::Edge(d0), new_index<CMap3::Edge>(m));
	// 	do{
	// 		d3 = phi1(m, d1);
	// 		do{
	// 			d4 = phi3(m, d3);
	// 			set_index<CMap3::Edge>(m, d3, index_of(m, CMap3::Edge(phi2(m, d3))));
	// 			set_index<CMap3::Edge>(m, d4, index_of(m, CMap3::Edge(phi2(m, d4))));
	// 			d3 = phi1(m, d3);
	// 		}while(d3 != d1);
	// 		d1 = phi<2, 3>(m, d1);
	// 	}while(d1 != d0);
	// }

	// if (is_indexed<CMap3::Edge2>(m))
	// {
	// 	do{
	// 	}while(d1 != d0);

	// 	do{
	// 		set_index<CMap3::Edge2>(m, CMap3::Edge2(d1), new_index<CMap3::Edge2>(m));
	// 		d3 = phi1(m, d1);
	// 		do{
	// 			d4 = phi3(m, d3);
	// 			set_index<CMap3::Edge2>(m, CMap3::Edge2(d3), new_index<CMap3::Edge2>(m));
	// 			set_index<CMap3::Edge2>(m, CMap3::Edge2(d4), new_index<CMap3::Edge2>(m));
	// 			d3 = phi1(m, d3);
	// 		}while(d3 != d1);
	// 		d1 = phi<2, 3>(m, d1);
	// 	}while(d1 != d0);
	// }

	// if (is_indexed<CMap3::Face>(m))
	// {
	// 	do{
	// 		set_index<CMap3::Face>(m, CMap3::Face(d1), new_index<CMap3::Face>(m));
	// 		d1 = phi<2, 3>(m, d1);
	// 	}while(d1 != d0);
	// }

	// if (is_indexed<CMap3::Face2>(m))
	// {
	// 	do{
	// 		set_index<CMap3::Face2>(m, CMap3::Face2(d1), new_index<CMap3::Face2>(m));
	// 		set_index<CMap3::Face2>(m, CMap3::Face2(d2), new_index<CMap3::Face2>(m));
	// 		d1 = phi3(m, d2);
	// 		d2 = phi2(m, d1);
	// 	}while(d1 != d0);
	// }

	// if (is_indexed<CMap3::Volume>(m))
	// {
	// 	do{
	// 		set_index<CMap3::Volume>(m, CMap3::Volume(d1), new_index<CMap3::Volume>(m));
	// 		d1 = phi3(m, d2);
	// 		d2 = phi2(m, d1);
	// 	}while(d1 != d0);
	// }
}

/// Pre cut hex, input: original vertex dart
inline CMap3::Vertex octosect_hex(CMap3& m, CMap3::Volume w)
{
	Dart d0, d10, d11, d20, d21, d22, d23;
	Dart d1, d2, d3, d4;
	std::vector<Dart> path0, path1, path2, path3;

	d0 = phi<1, 1>(m, w.dart);
	d10 = phi<2, 1>(m, d0);
	d11 = phi<1, 2>(m, w.dart);
	d20 = phi<1, 2, 1, 1>(m, d10);
	d21 = phi<2, 1, 2, 1>(m, d20);
	d22 = phi<1, 2, 1, 1>(m, d11);
	d23 = phi<2, 1, 2, 1>(m, d22);

	d1 = d0;
	do
	{
		path0.push_back(d1);
		d1 = phi<1, 2, 1>(m, d1);
	} while (d1 != d0);

	cut_volume(m, path0);
	path0.clear();
	CMap3::Vertex v = quadrangulate_face(m, CMap3::Face(phi2(m, d0)));

	d1 = d10;
	d2 = d11;
	do
	{
		path0.push_back(d1);
		path1.push_back(d2);
		d1 = phi<1, 2, 1>(m, d1);
		d2 = phi<1, 2, 1>(m, d2);
	} while (d1 != d10);
	cut_volume(m, path0);
	cut_volume(m, path1);
	path0.clear();
	path1.clear();

	d1 = phi<-1, 2>(m, d20);
	cut_face(m, CMap3::Vertex(d1), CMap3::Vertex(phi<1, 1, 1>(m, d1)));
	d1 = phi<-1, 2>(m, d22);
	cut_face(m, CMap3::Vertex(d1), CMap3::Vertex(phi<1, 1, 1>(m, d1)));

	d1 = d20;
	d2 = d21;
	d3 = d22;
	d4 = d23;
	do
	{
		path0.push_back(d1);
		path1.push_back(d2);
		path2.push_back(d3);
		path3.push_back(d4);
		d1 = phi<1, 2, 1>(m, d1);
		d2 = phi<1, 2, 1>(m, d2);
		d3 = phi<1, 2, 1>(m, d3);
		d4 = phi<1, 2, 1>(m, d4);
	} while (d1 != d20);
	cut_volume(m, path0);
	cut_volume(m, path1);
	cut_volume(m, path2);
	cut_volume(m, path3);

	return v;
}

template <typename MESH, typename FUNC1, typename FUNC2, typename FUNC3>
auto primal_cut_all_volumes(MESH& m, const FUNC1& on_edge_cut, const FUNC2& on_face_cut, const FUNC3& on_vol_cut)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMap3&>>
{
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Vertex2 = typename mesh_traits<MESH>::Vertex2;
	using Edge2 = typename mesh_traits<MESH>::Edge2;
	using Face2 = typename mesh_traits<MESH>::Face2;
	using Volume = typename mesh_traits<MESH>::Volume;
	static_assert(is_func_parameter_same<FUNC1, Vertex>::value, "Given function should take a Vertex");
	static_assert(is_func_parameter_same<FUNC2, Vertex>::value, "Given function should take a Vertex");
	static_assert(is_func_parameter_same<FUNC3, Vertex>::value, "Given function should take a Vertex");

	// CellCache<MESH> vol_cache(m);
	// vol_cache.template build<Volume>();

	CellCache<MESH> edge_vert_cache(m);
	CellCache<MESH> face_vert_cache(m);

	CellMarker<MESH, Volume> cm(m);
	foreach_cell(m, [&](Volume w) -> bool {
		cm.mark(w);
		return true;
	});

	quadrangulate_all_faces(
		m,
		[&](Vertex v) -> void {
			edge_vert_cache.add(v);
			on_edge_cut(v);
		},
		[&](Vertex v) -> void {
			face_vert_cache.add(v);
			on_face_cut(v);
		});

	foreach_cell(edge_vert_cache, [&](Vertex ve) -> bool {
		Dart d = ve.dart;
		do
		{
			if (!is_boundary(m, d) && cm.is_marked(Volume(d)))
			{
				Dart d01 = phi_1(m, d);
				Dart d02 = phi2(m, d01);
				Dart d21 = phi<2, 1>(m, d);
				Dart d22 = phi2(m, d21);

				Dart f0 = add_face(static_cast<CMap1&>(m), 4, false).dart;
				Dart f1 = add_face(static_cast<CMap1&>(m), 4, false).dart;
				Dart ee = f0;
				Dart ff = f1;
				do
				{
					phi3_sew(m, ee, ff);
					ee = phi1(m, ee);
					ff = phi_1(m, ff);
				} while (ee != f0);

				phi2_unsew(m, d01);
				phi2_unsew(m, d21);

				phi2_sew(m, d01, f0);
				phi2_sew(m, d02, f1);

				phi2_sew(m, d21, phi_1(m, f0));
				phi2_sew(m, d22, phi1(m, f1));
			}
			d = phi<2, 3>(m, d);
		} while (d != ve.dart);
		return true;
	});

	parallel_foreach_cell(face_vert_cache, [&](Vertex vf) -> bool {
		Dart d0 = vf.dart;
		Dart d1 = phi<2, 3, 2, 3>(m, vf.dart);
		Dart d;
		if (cm.is_marked(Volume(d0)))
		{
			d = d0;
			do
			{
				phi2_sew(m, phi<2, 1>(m, d), phi<-1, 2, -1>(m, d));
				d = phi<2, 3, 2, 1>(m, d);
			} while (d != d0);
		}

		if (!is_boundary(m, d1) && cm.is_marked(Volume(d1)))
		{
			d = d1;
			do
			{
				phi2_sew(m, phi<2, 1>(m, d), phi<-1, 2, -1>(m, d));
				d = phi<2, 3, 2, 1>(m, d);
			} while (d != d1);
		}
		return true;
	});

	CellCache<MESH> vol_vert_cache(m);
	foreach_cell(m, [&](Volume w) -> bool {
		Dart d0 = w.dart;
		Vertex vw = Vertex(phi<1, 2, -1>(m, d0));
		vol_vert_cache.add(vw);
		return true;
	});

	foreach_cell(vol_vert_cache, [&](Vertex v) -> bool {
		if (is_indexed<Vertex>(m))
			set_index(m, v, new_index<Vertex>(m));

		foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
			if (is_indexed<Edge>(m))
			{
				if (index_of(m, Edge(d)) == INVALID_INDEX)
					set_index(m, Edge(d), new_index<Edge>(m));
			}
			if (is_indexed<Face>(m))
			{
				if (index_of(m, Face(d)) == INVALID_INDEX)
					set_index(m, Face(d), new_index<Face>(m));
			}
			if (is_indexed<Volume>(m))
			{
				if (index_of(m, Volume(d)) == INVALID_INDEX)
					set_index(m, Volume(d), new_index<Volume>(m));
			}
			return true;
		});

		foreach_incident_volume(m, v, [&](Volume w) -> bool {
			foreach_dart_of_orbit(m, w, [&](Dart d) -> bool {
				if (is_indexed<HalfEdge>(m))
				{
					if (index_of(m, HalfEdge(d)) == INVALID_INDEX)
						set_index<HalfEdge>(m, HalfEdge(d), new_index<HalfEdge>(m));
				}
				if (is_indexed<Vertex2>(m))
				{
					if (index_of(m, Vertex2(d)) == INVALID_INDEX)
						set_index(m, Vertex2(d), new_index<Vertex2>(m));
				}
				if (is_indexed<Edge2>(m))
				{
					if (index_of(m, Edge2(d)) == INVALID_INDEX)
						set_index(m, Edge2(d), new_index<Edge2>(m));
				}
				if (is_indexed<Face2>(m))
				{
					if (index_of(m, Face2(d)) == INVALID_INDEX)
						set_index(m, Face2(d), new_index<Face2>(m));
				}
				return true;
			});
			return true;
		});
		return true;
	});

	if (is_indexed<Vertex>(m))
	{
		parallel_foreach_cell(face_vert_cache, [&](Vertex vf) -> bool {
			set_index<Vertex>(m, vf, index_of(m, vf));
			return true;
		});
		parallel_foreach_cell(edge_vert_cache, [&](Vertex ve) -> bool {
			set_index<Vertex>(m, ve, index_of(m, ve));
			return true;
		});
	}

	if (is_indexed<Edge>(m))
	{
		parallel_foreach_cell(face_vert_cache, [&](Vertex vf) -> bool {
			Dart d0 = vf.dart;
			Dart d = d0;
			do
			{
				set_index<Edge>(m, Edge(d), index_of(m, Edge(d)));
				d = phi<2, 3, 2, 1>(m, d);
			} while (d != d0);
			return true;
		});
	}

	foreach_cell(vol_vert_cache, [&](Vertex v) -> bool {
		on_vol_cut(v);
		return true;
	});
}

/* -------------------------- BUTTERFLY VOLUME MASKS -------------------------- */

// Remplit les vecteurs de p/q-points pour le volume incident à d
// Si le masque est mal défini (proche du bord) on ne garde que les p-points
// et le vecteur de q-points est vide après l'appel
template <typename MESH>
auto volumePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	// p-points : ajouter les sommets de Face(d) + Face(phi<2, 1, 1, 2>(m,d))
	Dart b = d;
	Dart t = b;
	do
	{
		p_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);
	b = phi<2, 1, 1, 2>(m, d);
	t = b;
	do
	{
		p_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);

	// q-points : pour chaque face f du volume,
	// ajouter Face(Inherit::phi<3, 2, 1, 1, 2>(f)) (face opposée du volume adjacent)
	// si le volume adjacent n'existe pas (bord) alors le masque est mal défini,
	// donc on ignore les q-points
	Dart f = d;
	unsigned i = 0;
	do
	{
		// S'il n'y a pas de q-points pour cette face,
		// on vide le vecteur de q-points et on quitte la recherche de q-points
		if (is_boundary(m, phi3(m, f)))
		{
			q_point.clear();
			break;
		}
		// Sinon, récupération des brins de la face opposée du volume adjacent
		b = phi<3, 2, 1, 1, 2>(m, f);
		t = b;
		do
		{
			q_point.push_back(t);
			t = phi1(m, t);
		} while (t != b);
		// face suivante sur l'hexaèdre
		f = phi2(m, (i % 2 == 0 ? phi1(m, f) : phi_1(m, f)));
		++i;
	} while (f != d);
}

// Remplit les vecteurs de p/q/r/s/t-points pour la face incidente à d
// Si le masque est mal défini (proche du bord) on ne garde que les p-points
// et les autres vecteurs sont vides après l'appel
template <typename MESH>
auto facePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point,
				   std::vector<Dart>& r_point, std::vector<Dart>& s_point, std::vector<Dart>& t_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	// p-points : ajouter les sommets de Face(d)
	Dart t = d;
	foreach_incident_vertex(m, typename MESH::Face(d), [&](typename MESH::Vertex dd) -> bool {
		p_point.push_back(dd.dart);
		return true;
	});

	// q-points : ajouter les sommets de Face(Inherit::phi<2, 1, 1, 2>(d)) + Face(Inherit::phi<3, 2, 1, 1, 2>(d))
	Dart b = phi<2, 1, 1, 2>(m, d);
	t = b;
	do
	{
		q_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);
	b = phi<3, 2, 1, 1, 2>(m, d);
	t = b;
	do
	{
		q_point.push_back(t);
		t = phi1(m, t);
	} while (t != b);

	// r/t-points
	// Question : si un volume adjacent au volume sur ou sous la face n'existe pas, doit-on chercher
	// les r/s/t-points par d'autres chemins ? (nombre de config énorme...)

	// Pour chaque brin de p-point
	for (Dart p : p_point)
	{
		// Est-ce que les volumes adjacents aux volumes au dessus et au dessous existent ou est-ce en surface
		bool face_phi2_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi2(m, p))),
			 face_phi32_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi<3, 2>(m, p)));

		Dart phi323211 = phi<3, 2, 3, 2, 1, 1>(m, p), phi3232111 = phi<3232111>(m, p),
			 phi23211 = phi<2, 3, 2, 1, 1>(m, p), phi232111 = phi<2, 3, 2, 1, 1, 1>(m, p);

		// Si Face(phi2(p)) est sur la surface et Face(phi<3, 2>(p)) aussi
		// alors il n'y a pas de r/t-points
		if (face_phi2_is_surface && face_phi32_is_surface)
		{
			q_point.clear();
			r_point.clear();
			s_point.clear();
			t_point.clear();
			break;
		}

		// Si Face(phi2(p)) est sur la surface et pas Face(phi<3, 2>(p))
		// alors il y a deux r-points : Edge(phi<3, 2, 3, 2, 1 ,1>(p))
		if (face_phi2_is_surface && !face_phi32_is_surface)
		{
			r_point.push_back(phi323211);
			r_point.push_back(phi3232111);
			continue;
		}

		// Si Face(phi2(p)) n'est pas sur la surface mais que Face(phi<3, 2>(p)) l'est
		// alors il y a deux r-points : Edge(phi<2, 3, 2, 1, 1>(p))
		if (!face_phi2_is_surface && face_phi32_is_surface)
		{
			r_point.push_back(phi23211);
			r_point.push_back(phi232111);
			continue;
		}

		// Ni Face(phi2(p)) ni Face(phi<3, 2>(p)) ne sont sur la surface
		// alors soit Edge(phi<2, 3, 2, 1, 1>(p)) == Edge(phi<3, 2, 3, 2, 1, 1>(p)) (topologie régulière)
		// soit Edge(phi<2, 3, 2, 1, 1>(p)) != Edge(phi<3, 2, 3, 2, 1, 1>(p)) (arete extraordinaire)

		// Si Edge(phi<2, 3, 2, 1, 1>(p)) != Edge(phi<3, 2, 3, 2, 1, 1>(p)) (arete extraordinaire)
		// alors il y a 4 t-points : Edge(phi<2, 3, 2, 1, 1>(p)) + Edge(phi<3, 2, 3, 2, 1, 1>(p))
		if (index_of(m, typename MESH::Edge(phi23211)) != index_of(m, typename MESH::Edge(phi323211)))
		{
			t_point.push_back(phi23211);
			t_point.push_back(phi232111);
			t_point.push_back(phi323211);
			t_point.push_back(phi3232111);
		}
		// Sinon Edge(phi<2, 3, 2, 1, 1>(p)) == Edge(phi<3, 2, 3, 2, 1, 1>(p))
		// et donc il y a 2 r-points
		else
		{
			r_point.push_back(phi23211);
			r_point.push_back(phi232111);
		}
	}

	// s-points : pour chaque q, ajouter Edge(phi<2, 3, 2, 1, 1>(q))
	for (Dart q : q_point)
	{
		// S'il n'y a pas de s-points
		if (is_boundary(m, phi<2, 3>(m, q)))
		{
			q_point.clear();
			r_point.clear();
			s_point.clear();
			t_point.clear();
			break;
		}
		s_point.push_back(phi<2, 3, 2, 1, 1>(m, q));
		s_point.push_back(phi<2, 3, 2, 1, 1, 1>(m, q));
	}
}

// Remplit les vecteurs de p/q/r/s-points pour l'arete incidente à d
// Si le masque est mal défini (proche du bord) on ne garde que les p-points
// et les autres vecteurs sont vides après l'appel
// Attention, un appel récursif est fait (un seul)
template <typename MESH>
auto edgePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point,
				   std::vector<Dart>& r_point, std::vector<Dart>& s_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	if (is_incident_to_boundary(m, typename MESH::Edge(d)))
	{
		std::cerr << "Fatal error : edgePointMask called on dart which is surface" << std::endl;
		exit(EXIT_FAILURE);
	}

	// pour savoir si c'est le premier appel ou le rappel pour la deuxieme partie des r-points
	bool first_call = p_point.size() == 0;

	// p-points : ajouter les extremités de l'arete (sauf si c'est le deuxième appel
	if (first_call)
	{
		p_point.push_back(d);
		p_point.push_back(phi1(m, d));
	}

	// Pour toutes les faces incidentes à l'arete (au moins 4)
	// (arret garanti car l'arete n'est pas en surface)
	Dart t = d;
	do
	{
		// q-points : ajouter l'extremité de l'arete opposée (autre extremité : appel récursif)
		q_point.push_back(phi<1, 1>(m, t));

		// r/s-points du coté sens direct
		// Question : s'il n'y a pas de volume adjacent à celui de la face courante,
		// doit-on chercher le point r/s par d'autres chemins (nombre de config!)

		bool face_phi12_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi<1, 2>(m, t))),
			 face_phi132_is_surface = is_incident_to_boundary(m, typename MESH::Face(phi<1, 3, 2>(m, t)));

		Dart phi1323211 = phi<1, 3, 2, 3, 2, 1, 1>(m, t), phi1232_1 = phi<1, 2, 3, 2, -1>(m, t);

		// Si Face(phi<1, 2>(t)) et Face(phi<1, 3, 2>(t)) sont sur la surface
		// alors pas de r/s-points coté direct et on peux donc arreter la recherche le points
		if (face_phi12_is_surface && face_phi132_is_surface)
		{
			q_point.clear();
			r_point.clear();
			s_point.clear();
			break;
		}
		// Au moins l'une des face n'est pas sur la surface
		// Si ni Face(phi<1, 2>(t)) ni Face(phi<1, 3, 2>(t)) ne sont sur la surface
		// alors soit Vertex(phi<1, 3, 2, 3, 2, 1, 1>(t)) == Vertex(phi<1, 2, 3, 2, -1>(t)) (topologie régulière)
		// soit Vertex(phi<1, 3, 2, 3, 2, 1, 1>(t)) != Vertex(phi<1, 2, 3, 2, -1>(t)) (arete extraordinaire)
		if (!face_phi12_is_surface && !face_phi132_is_surface)
		{
			// Si Vertex(phi<1, 3, 2, 3, 2, 1, 1>(t)) != Vertex(phi<1, 2, 3, 2, -1>(t)) (arete extraordinaire)
			// alors 2 s-points Vertex(phi<1, 3, 2, 3, 2, 1, 1>(t)) et Vertex(phi<1, 2, 3, 2, -1>(t))
			if (index_of(m, typename MESH::Vertex(phi1323211)) != index_of(m, typename MESH::Vertex(phi1232_1)))
			{
				s_point.push_back(phi1323211);
				s_point.push_back(phi1232_1);
			}
			// Sinon Vertex(phi<1, 3, 2, 3, 2, 1, 1>(t)) == Vertex(phi<1, 2, 3, 2, -1>(t)) (topologie régulière)
			// alors 1 r-point
			else
			{
				r_point.push_back(phi1323211);
			}
		}
		// sinon si Face(phi<1, 2>(t)) est sur la surface (et pas Face(phi<1, 3, 2>(t)))
		// alors 1 r-point : Vertex(phi<1, 3, 2, 3, 2, 1, 1>(t))
		else if (face_phi12_is_surface)
		{
			r_point.push_back(phi1323211);
		}
		// sinon (si Face(phi<1, 3, 2>(t)) est sur la surface et pas Face(phi<1, 2>(t))))
		// alors 1 r-point : Vertex(phi<1, 2, 3, 2, -1>(t))
		else
		{
			r_point.push_back(phi1232_1);
		}

		// Face incidente suivante
		t = phi<2, 3>(m, t);
	} while (t != d);

	// Sinon, faire la meme chose dans l'autre sens (si ce n'est pas déjà fait)
	// pour avoir les r/s-points de l'autre coté de l'arete
	if (q_point.size() != 0 && first_call)
	{
		edgePointMask(m, phi2(m, d), p_point, q_point, r_point, s_point);
	}
}

/* ------------------------- BUTTERFLY SURFACE MASKS ------------------------- */

template <typename MESH>
auto surfaceAdjacentFace(const MESH& m, Dart d) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Dart>
{
	if (!is_boundary(m, phi3(m, d)))
	{
		std::cerr << "Fatal error : surfaceAdjacentFace called on dart which is not surface" << std::endl;
		exit(EXIT_FAILURE);
	}
	Dart t = phi2(m, d);
	while (!is_boundary(m, phi3(m, t)))
	{
		t = phi<3, 2>(m, t);
	}
	return t;
}

template <typename MESH>
auto surfaceFacePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	Dart t = d;
	do
	{
		p_point.push_back(t);
		Dart e = surfaceAdjacentFace(m, t);
		q_point.push_back(phi<1, 1>(m, e));
		q_point.push_back(phi_1(m, e));
		t = phi1(m, t);
	} while (t != d);
}

template <typename MESH>
auto surfaceEdgePointMask(const MESH& m, Dart d, std::vector<Dart>& p_point, std::vector<Dart>& q_point,
						  std::vector<Dart>& r_point) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	p_point.push_back(d);
	p_point.push_back(phi1(m, d));

	q_point.push_back(phi<1, 1>(m, d));
	q_point.push_back(phi_1(m, d));
	Dart t = surfaceAdjacentFace(m, phi1(m, d));
	r_point.push_back(phi_1(m, t));
	t = surfaceAdjacentFace(m, phi_1(m, d));
	r_point.push_back(phi<1, 1>(m, t));

	Dart e = surfaceAdjacentFace(m, d);
	q_point.push_back(phi<1, 1>(m, e));
	q_point.push_back(phi_1(m, e));
	t = surfaceAdjacentFace(m, phi1(m, e));
	r_point.push_back(phi_1(m, t));
	t = surfaceAdjacentFace(m, phi_1(m, e));
	r_point.push_back(phi<1, 1>(m, t));
}

/* ----------------------------- BUTTERFLY RULES ----------------------------- */

template <typename MESH, typename T>
auto sum(const MESH& m, const std::vector<Dart>& v, typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	Vec3 s = Vec3(0.f, 0.f, 0.f);
	for (Dart i : v)
		s += value<T>(m, attribute, typename MESH::Vertex(i));
	return s;
}

template <typename T, typename MESH>
auto volumePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
					 typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	if (q_point.size() == 0)
	{
		return (1.f / 8) * sum<MESH>(m, p_point, attribute);
	}
	else
	{
		return ((6 * _W_ + 1) / 8) * sum<MESH>(m, p_point, attribute) - (_W_ / 4) * sum<MESH>(m, q_point, attribute);
	}
}

template <typename T, typename MESH>
auto facePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
				   const std::vector<Dart>& r_point, const std::vector<Dart>& s_point, const std::vector<Dart>& t_point,
				   typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	if (q_point.size() == 0)
	{
		return (1.f / 4) * sum<MESH>(m, p_point, attribute);
	}
	else
	{
		float w1 = (2 * _W_ + 1) / 4, w2 = _W_ / 4, w3 = _W_ / 8;
		return w1 * sum<MESH>(m, p_point, attribute) + w2 * sum<MESH>(m, q_point, attribute) -
			   w2 * sum<MESH>(m, r_point, attribute) - w3 * sum<MESH>(m, s_point, attribute) -
			   w3 * sum<MESH>(m, t_point, attribute);
	}
}

template <typename T, typename MESH>
auto edgePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
				   const std::vector<Dart>& r_point, const std::vector<Dart>& s_point,
				   typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	int32 N = uint32(q_point.size()) / 2;
	float32 w1 = 1.f / 2;
	if (N == 0)
	{
		return w1 * sum<MESH>(m, p_point, attribute);
	}
	else
	{
		float w2 = _W_ / N, w3 = _W_ / (2 * N);
		return w1 * sum<MESH>(m, p_point, attribute) + w2 * sum<MESH>(m, q_point, attribute) -
			   w2 * sum<MESH>(m, r_point, attribute) - w3 * sum<MESH>(m, s_point, attribute);
	}
}

template <typename T, typename MESH>
auto surfaceFacePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
						  typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	Vec3 sum_p = sum<MESH>(m, p_point, attribute), sum_q = sum<MESH>(m, q_point, attribute);
	float wp = 1.f / 4.f + _W_, wq = -_W_ / 2.f;
	return wp * sum_p + wq * sum_q;
}

template <typename T, typename MESH>
auto surfaceEdgePointRule(const MESH& m, const std::vector<Dart>& p_point, const std::vector<Dart>& q_point,
						  const std::vector<Dart>& r_point,
						  typename mesh_traits<MESH>::template Attribute<T>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Vec3>
{
	static const float _W_ = 1.f / 16;
	Vec3 sum_p = sum<MESH>(m, p_point, attribute), sum_q = sum<MESH>(m, q_point, attribute),
		 sum_r = sum<MESH>(m, r_point, attribute);
	float wp = 1.f / 2.f, wq = _W_ / 2.f, wr = -_W_ / 2.f;
	return wp * sum_p + wq * sum_q + wr * sum_r;
}

template <typename MESH>
auto subdivideEdge(MESH& m, Dart d, Vec3& p, typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	cut_edge(m, typename MESH::Edge(d));
	value<Vec3>(m, attribute, typename MESH::Vertex(phi1(m, d))) = p;
}

// /!\ The edges of the face are already cut
template <typename MESH>
auto subdivideFace(MESH& m, Dart d, Vec3& p, typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, Dart>
{
	using Vertex = typename MESH::Vertex;
	Dart e = phi1(m, d);
	Dart f = phi<1, 1, 1, 1>(m, e);
	cut_face(m, Vertex(e), Vertex(f));
	subdivideEdge(m, phi_1(m, e), p, attribute);

	Dart result = phi_1(m, e);

	cut_face(m, Vertex(phi_1(m, e)), Vertex(phi<1, 1>(m, e)));
	cut_face(m, Vertex(phi_1(m, f)), Vertex(phi<1, 1>(m, f)));

	return result;
}

// /!\ The faces of the volume are already cut
template <typename MESH>
auto subdivideVolume(MESH& m, Dart d, Vec3& p, typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Vertex = typename MESH::Vertex;
	Dart first_cut_dir = phi1(m, d);
	Dart left_cut_dir = phi1(m, first_cut_dir), right_cut_dir = phi<2, 1, 1, 1, 2>(m, first_cut_dir);

	Dart cut_direction[4];
	cut_direction[0] = phi<1, 2, 1, 1>(m, left_cut_dir);
	cut_direction[1] = phi<1, 2, 1, 2, 1, 1, 1>(m, left_cut_dir);
	cut_direction[2] = phi<1, 2, 1, 1>(m, right_cut_dir);
	cut_direction[3] = phi<1, 2, 1, 2, 1, 1, 1>(m, right_cut_dir);

	auto cutVolume = [&](Dart d) -> typename MESH::Face {
		std::vector<Dart> cut_path;
		Dart t = d;
		do
		{
			cut_path.push_back(t);
			t = phi<1, 2, 1>(m, t);
		} while (t != d);
		return cut_volume(m, cut_path);
	};

	typename MESH::Face F = cutVolume(first_cut_dir);
	subdivideFace(m, F.dart, p, attribute);

	F = cutVolume(left_cut_dir);
	cut_face(m, Vertex(F.dart), Vertex(phi<1, 1, 1>(m, F.dart)));
	F = cutVolume(right_cut_dir);
	cut_face(m, Vertex(F.dart), Vertex(phi<1, 1, 1>(m, F.dart)));
	for (int i = 0; i < 4; ++i)
	{
		cutVolume(cut_direction[i]);
	}
}

template <typename MESH>
auto subdivideListEdges(MESH& m, std::vector<Dart>& edges, std::queue<Vec3>& edge_points,
						typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	for (Dart d : edges)
	{
		subdivideEdge(m, d, edge_points.front(), attribute);
		edge_points.pop();
	}
}

inline void subdivideListEdges(CPH3& m, std::vector<Dart>& edges, std::queue<Vec3>& edge_points,
							   typename mesh_traits<CPH3>::template Attribute<Vec3>* attribute)
{
	CPH3 m2(m);
	for (Dart d : edges)
	{
		m2.current_level_ = m.edge_level(d) + 1;
		subdivideEdge(m2, d, edge_points.front(), attribute);
		edge_points.pop();
	}
}

template <typename MESH>
auto subdivideListFaces(MESH& m, std::vector<Dart>& faces, std::queue<Vec3>& face_points,
						typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	for (Dart d : faces)
	{
		subdivideFace(m, d, face_points.front(), attribute);
		face_points.pop();
	}
}

inline void subdivideListFaces(CPH3& m, std::vector<Dart>& faces, std::queue<Vec3>& face_points,
							   typename mesh_traits<CPH3>::template Attribute<Vec3>* attribute)
{
	CPH3 m2(m);
	for (Dart d : faces)
	{
		m2.current_level_ = m.face_level(d) + 1;
		subdivideFace(m2, d, face_points.front(), attribute);
		face_points.pop();
	}
}

template <typename MESH>
auto subdivideListVolumes(MESH& m, std::vector<Dart>& volumes, std::queue<Vec3>& volume_points,
						  typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	for (Dart d : volumes)
	{
		subdivideVolume(m, d, volume_points.front(), attribute);
		volume_points.pop();
	}
}

inline void subdivideListVolumes(CPH3& m, std::vector<Dart>& volumes, std::queue<Vec3>& volume_points,
								 typename mesh_traits<CPH3>::template Attribute<Vec3>* attribute)
{
	CPH3 m2(m);
	for (Dart d : volumes)
	{
		m2.current_level_ = m.volume_level(d) + 1;
		subdivideVolume(m2, d, volume_points.front(), attribute);
		volume_points.pop();
	}
}

template <typename MESH>
auto butterflySubdivisionVolumeAdaptative(MESH& m, double angle_threshold,
										  typename mesh_traits<MESH>::template Attribute<Vec3>* attribute)
	-> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>
{
	using Volume = typename MESH::Volume;
	using Face = typename MESH::Face;
	using Edge = typename MESH::Edge;

	CellMarker<MESH, Edge> cm_surface(m);
	CellMarker<MESH, Volume> cm_cell(m);

	for (Dart d = m.begin(); d != m.end(); d = m.next(d))
	{
		if (is_boundary(m, phi3(m, d)) && !cm_surface.is_marked(Edge(d)))
		{
			Dart d2 = phi2(m, d);
			while (!is_boundary(m, phi3(m, d2)))
			{
				d2 = phi<3, 2>(m, d2);
			}
			auto edge_angle = geometry::angle(m, typename MESH::Face2(d), typename MESH::Face2(d2), attribute);
			if (std::abs(edge_angle) > angle_threshold)
			{
				if (!cm_cell.is_marked(Volume(d)))
					cm_cell.mark(Volume(d));
				if (!cm_cell.is_marked(Volume(d2)))
					cm_cell.mark(Volume(d2));
			}
			cm_surface.mark(Edge(d));
		}
	}

	CellMarker<MESH, Edge> cm_edge(m);
	CellMarker<MESH, Face> cm_face(m);
	CellMarker<MESH, Volume> cm_volume(m);
	std::queue<Vec3> volume_points, face_points, edge_points;
	std::vector<Dart> edges, faces, volumes;
	std::vector<Dart> p_point, q_point, r_point, s_point, t_point;

	// computing the new vertices's embedding
	for (Dart d = m.begin(); d != m.end(); d = m.next(d))
	{
		if (!is_boundary(m, d) && cm_cell.is_marked(Volume(d)))
		{
			foreach_dart_of_orbit(m, Volume(d), [&](Dart t) -> bool {
				// edges vertices
				if (!cm_edge.is_marked(Edge(t)))
				{
					if (is_incident_to_boundary(m, Edge(t)) && !is_boundary(m, phi3(m, t)))
					{
						// we ignore the surfaces darts which are in junction of two volumes
					}
					// surface case
					else if (is_boundary(m, phi3(m, t)))
					{
						p_point.clear();
						q_point.clear();
						r_point.clear();
						surfaceEdgePointMask(m, t, p_point, q_point, r_point);

						Vec3 vec = surfaceEdgePointRule<Vec3>(m, p_point, q_point, r_point, attribute);

						edge_points.push(vec);
						edges.push_back(t);
						cm_edge.mark(Edge(t));
					}
					// volume case
					else
					{
						p_point.clear();
						q_point.clear();
						r_point.clear();
						s_point.clear();
						edgePointMask(m, t, p_point, q_point, r_point, s_point);

						Vec3 vec = edgePointRule<Vec3>(m, p_point, q_point, r_point, s_point, attribute);

						edge_points.push(vec);

						edges.push_back(t);
						cm_edge.mark(Edge(t));
					}
				}
				// faces vertices
				if (!cm_face.is_marked(Face(t)))
				{
					// cas surface
					if (is_incident_to_boundary(m, Face(t)))
					{
						p_point.clear();
						q_point.clear();
						surfaceFacePointMask(m, t, p_point, q_point);

						Vec3 vec = surfaceFacePointRule<Vec3>(m, p_point, q_point, attribute);

						face_points.push(vec);
					}
					// volume case
					else
					{
						p_point.clear();
						q_point.clear();
						r_point.clear();
						s_point.clear();
						t_point.clear();
						facePointMask(m, t, p_point, q_point, r_point, s_point, t_point);

						Vec3 vec = facePointRule<Vec3>(m, p_point, q_point, r_point, s_point, t_point, attribute);

						face_points.push(vec);
					}
					faces.push_back(t);
					cm_face.mark(Face(t));
				}
				// volumes vertices
				if (!cm_volume.is_marked(Volume(t)))
				{
					p_point.clear();
					q_point.clear();
					volumePointMask(m, t, p_point, q_point);

					Vec3 vec = volumePointRule<Vec3>(m, p_point, q_point, attribute);

					volume_points.push(vec);
					volumes.push_back(t);
					cm_volume.mark(Volume(t));
				}
				return true;
			});
			cm_cell.unmark(Volume(d));
		}
	}

	subdivideListEdges(m, edges, edge_points, attribute);
	subdivideListFaces(m, faces, face_points, attribute);
	subdivideListVolumes(m, volumes, volume_points, attribute);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_H_
