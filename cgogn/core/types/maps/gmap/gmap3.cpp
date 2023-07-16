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
#include <cgogn/core/types/maps/gmap/gmap3.h>

#include <cgogn/core/types/cell_marker.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

GMap3::Vertex cut_edge(GMap3& m, GMap3::Edge e, bool set_indices)
{
	Dart d0 = e.dart;
	Dart d23 = beta<2, 3>(m, d0);

	GMap3::Vertex v(cut_edge(static_cast<GMap2&>(m), GMap2::Edge(d0), false).dart);

	while (d23 != e.dart)
	{
		d0 = d23;
		d23 = beta<2, 3>(m, d23);
		Dart d1 = beta0(m, d0);
		cut_edge(static_cast<GMap2&>(m), GMap2::Edge(d0), false);
		beta3_sew(m, beta0(m, d0), beta0(m, beta3(m, d0)));
		beta3_sew(m, beta0(m, d1), beta0(m, beta3(m, d1)));
	}

	Dart e1 = beta0(m, e.dart);
	Dart e2 = beta0(m, beta3(m, e.dart));
	beta3_sew(m, e1, e2);
	beta3_sew(m, beta1(m, e1), beta1(m, e2));

	if (set_indices)
	{
		if (is_indexed<GMap3::Vertex>(m))
			set_index(m, v, new_index<GMap3::Vertex>(m));

		if (is_indexed<GMap3::Vertex2>(m))
		{
			Dart di = v.dart;
			do
			{
				if (!is_boundary(m, di))
					set_index(m, GMap3::Vertex2(di), new_index<GMap3::Vertex2>(m));
				di = beta<2, 3>(m, di);
			} while (di != v.dart);
		}

		if (is_indexed<GMap3::Edge>(m))
		{
			uint32 emb1 = index_of(m, v);
			uint32 emb2 = new_index<GMap3::Edge>(m);
			set_index(m, GMap3::Edge(v.dart), emb1);
			set_index(m, GMap3::Edge(beta1(m, v.dart)), emb2);
		}

		if (is_indexed<GMap3::Face2>(m))
		{
			foreach_dart_of_orbit(m, v, [&](Dart di) -> bool {
				if (is_boundary(m, di))
					return true;
				copy_index<GMap3::Face>(m, di, beta0(m, di));
				return true;
			});
		}

		if (is_indexed<GMap3::Face>(m))
		{
			foreach_dart_of_orbit(m, v, [&](Dart di) -> bool {
				if (is_boundary(m, di))
					return true;
				copy_index<GMap3::Face>(m, di, beta0(m, di));
				return true;
			});
		}

		if (is_indexed<GMap3::Volume>(m))
		{
			foreach_dart_of_orbit(m, v, [&](Dart di) -> bool {
				if (is_boundary(m, di))
					return true;
				std::cout << "Dart of v" << di << " : beta0 " << beta0(m, di) << " -> "
						  << index_of(m, GMap3::Volume(beta0(m, di))) << std::endl;
				copy_index<GMap3::Volume>(m, di, beta0(m, di));
				return true;
			});
		}
	}

	return v;
}

GMap3::Edge cut_face(GMap3& m, GMap3::Vertex v1, GMap3::Vertex v2, bool set_indices)
{
	Dart d = v1.dart;
	Dart e = v2.dart;

	Dart dd = beta3(m, d);
	Dart ee = beta3(m, e);

	GMap2::Edge e1 = cut_face(static_cast<GMap2&>(m), GMap2::Vertex(d), GMap2::Vertex(e), false);
	GMap2::Edge e2 = cut_face(static_cast<GMap2&>(m), GMap2::Vertex(dd), GMap2::Vertex(ee), false);

	Dart e10 = beta0(m, e1.dart);
	Dart e20 = beta0(m, e2.dart);

	// beta3_sew(m, e1.dart, e20);
	// beta3_sew(m, e2.dart, e10);
	// beta3_sew(m, beta2(m, e1.dart), beta2(m, e20));
	// beta3_sew(m, beta2(m, e2.dart), beta2(m, e10));
	beta3_sew(m, e1.dart, e2.dart);
	beta3_sew(m, e10, e20);
	beta3_sew(m, beta2(m, e1.dart), beta2(m, e2.dart));
	beta3_sew(m, beta2(m, e10), beta2(m, e20));

	if (set_indices)
	{
		if (is_indexed<GMap3::Vertex>(m))
		{
			uint32 emb1 = index_of(m, v1);
			Dart ei = e1.dart;
			if (!is_boundary(m, ei))
			{
				set_index<GMap3::Vertex>(m, ei, emb1);
				ei = beta2(m, ei);
				set_index<GMap3::Vertex>(m, ei, emb1);
				ei = beta3(m, ei);
			}
			else
				ei = beta<2, 3>(m, ei);
			if (!is_boundary(m, ei))
			{
				set_index<GMap3::Vertex>(m, ei, emb1);
				ei = beta2(m, ei);
				set_index<GMap3::Vertex>(m, ei, emb1);
			}

			uint32 emb2 = index_of(m, v2);
			Dart ej = beta0(m, e1.dart);
			if (!is_boundary(m, ej))
			{
				set_index<GMap3::Vertex>(m, ej, emb2);
				ej = beta2(m, ej);
				set_index<GMap3::Vertex>(m, ej, emb2);
				ej = beta3(m, ej);
			}
			else
				ej = beta<2, 3>(m, ej);
			if (!is_boundary(m, ej))
			{
				set_index<GMap3::Vertex>(m, ej, emb2);
				ej = beta2(m, ej);
				set_index<GMap3::Vertex>(m, ej, emb2);
			}
		}

		if (is_indexed<GMap3::Vertex2>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<GMap3::Vertex2>(m, phi_1(m, d), e);
				copy_index<GMap3::Vertex2>(m, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<GMap3::Vertex2>(m, phi_1(m, dd), ee);
				copy_index<GMap3::Vertex2>(m, phi_1(m, ee), dd);
			}
		}
		if (is_indexed<GMap3::Edge>(m))
			set_index(m, GMap3::Edge(e10), new_index<GMap3::Edge>(m));

		if (is_indexed<GMap3::Face2>(m))
		{
			if (!is_boundary(m, e1.dart))
			{
				copy_index<GMap3::Face2>(m, e1.dart, beta1(m, e1.dart));
				copy_index<GMap3::Face2>(m, e10, beta1(m, e10));
				set_index(m, GMap3::Face2(beta2(m, e10)), new_index<GMap3::Face2>(m));
			}
			if (!is_boundary(m, e2.dart))
			{
				copy_index<GMap3::Face2>(m, e2.dart, beta1(m, e2.dart));
				copy_index<GMap3::Face2>(m, e20, beta1(m, e20));
				set_index(m, GMap3::Face2(beta2(m, e20)), new_index<GMap3::Face2>(m));
			}
		}

		if (is_indexed<GMap3::Face>(m))
		{
			uint32 emb = index_of(m, GMap3::Face(v1.dart));
			set_index<GMap3::Face>(m, e1.dart, emb);
			set_index<GMap3::Face>(m, e10, emb);
			set_index<GMap3::Face>(m, e2.dart, emb);
			set_index<GMap3::Face>(m, e20, emb);
			set_index(m, GMap3::Face(beta2(m, e10)), new_index<GMap3::Face>(m));
		}

		if (is_indexed<GMap3::Volume>(m))
		{
			if (!is_boundary(m, d))
			{
				uint32 emb = index_of(m, GMap3::Volume(d));
				foreach_dart_of_orbit(m, GMap3::Edge2(e1.dart), [&](Dart di) {
					set_index<GMap3::Volume>(m, di, emb);
					return true;
				});
			}
			if (!is_boundary(m, dd))
			{
				uint32 emb = index_of(m, GMap3::Volume(dd));
				foreach_dart_of_orbit(m, GMap3::Edge2(e2.dart), [&](Dart di) {
					set_index<GMap3::Volume>(m, di, emb);
					return true;
				});
			}
		}
	}

	return GMap3::Edge(e1.dart);
}

GMap3::Face cut_volume(GMap3& m, const std::vector<Dart>& path, bool set_indices)
{
	Dart f0 = add_face(static_cast<GMap1&>(m), uint32(path.size()), false).dart;
	Dart f1 = add_face(static_cast<GMap1&>(m), uint32(path.size()), false).dart;

	for (Dart d0 : path)
	{
		Dart d1 = phi2(m, d0);
		phi2_unsew(m, d0);

		phi2_sew(m, d0, f0);
		phi2_sew(m, d1, f1);

		phi3_sew(m, f0, f1);

		f0 = phi_1(m, f0);
		f1 = phi1(m, f1);
	}
	f0 = phi_1(m, f0);

	if (set_indices)
	{
		if (is_indexed<GMap3::Vertex>(m))
		{
			foreach_dart_of_orbit(m, GMap3::Face(f0), [&](Dart d) -> bool {
				Dart dd = beta2(m, d);
				copy_index<GMap3::Vertex>(m, d, dd);
				return true;
			});
		}
		if (is_indexed<GMap3::Vertex2>(m))
		{
			foreach_dart_of_orbit(m, GMap3::Face(f0), [&](Dart d) -> bool {
				copy_index<GMap3::Vertex2>(m, d, beta2(m, d));
				return true;
			});
		}
		if (is_indexed<GMap3::Edge>(m))
		{
			foreach_dart_of_orbit(m, GMap3::Face(f0), [&](Dart d) -> bool {
				copy_index<GMap3::Edge>(m, d, beta2(m, d));
				return true;
			});
		}
		if (is_indexed<GMap3::Face2>(m))
		{
			set_index(m, GMap3::Face2(f0), new_index<GMap3::Face2>(m));
			set_index(m, GMap3::Face2(beta3(m, f0)), new_index<GMap3::Face2>(m));
		}

		if (is_indexed<GMap3::Face>(m))
		{
			set_index(m, GMap3::Face(f0), new_index<GMap3::Face>(m));
		}

		if (is_indexed<GMap3::Volume>(m))
		{
			foreach_dart_of_orbit(m, GMap3::Face2(f0), [&](Dart d) -> bool {
				copy_index<GMap3::Volume>(m, d, beta2(m, d));
				return true;
			});
			set_index(m, GMap3::Volume(f1), new_index<GMap3::Volume>(m));
		}
	}

	return GMap3::Face(f0);
}

GMap3::Volume close_hole(GMap3& m, Dart d, bool set_indices)
{
	cgogn_message_assert(phi3(m, d) == d, "GMap3: close hole called on a dart that is not a phi3 fix point");

	DartMarkerStore marker(m);
	DartMarkerStore hole_volume_marker(m);

	std::vector<Dart> visited_faces;
	visited_faces.reserve(1024u);

	visited_faces.push_back(d);
	foreach_dart_of_orbit(m, GMap3::Face2(d), [&](Dart fd) -> bool {
		marker.mark(fd);
		return true;
	});

	uint32 count = 0u;

	for (uint32 i = 0u; i < uint32(visited_faces.size()); ++i)
	{
		const Dart it = visited_faces[i];
		Dart f = it;

		GMap3::Face2 hf = add_face(static_cast<GMap1&>(m), codegree(m, GMap3::Face(f)));
		foreach_dart_of_orbit(m, hf, [&](Dart fd) -> bool {
			hole_volume_marker.mark(fd);
			return true;
		});

		++count;

		Dart bit = hf.dart;
		do
		{
			Dart e = phi3(m, phi2(m, f));
			bool found = false;
			do
			{
				if (phi3(m, e) == e)
				{
					found = true;
					if (!marker.is_marked(e))
					{
						visited_faces.push_back(e);
						foreach_dart_of_orbit(m, GMap3::Face2(e), [&](Dart fd) -> bool {
							marker.mark(fd);
							return true;
						});
					}
				}
				else
				{
					if (hole_volume_marker.is_marked(e))
					{
						found = true;
						phi2_sew(m, e, bit);
					}
					else
						e = phi3(m, phi2(m, e));
				}
			} while (!found);

			phi3_sew(m, f, bit);
			bit = phi_1(m, bit);
			f = phi1(m, f);
		} while (f != it);
	}

	GMap3::Volume hole(phi3(m, d));

	if (set_indices)
	{
		foreach_dart_of_orbit(m, hole, [&](Dart hd) -> bool {
			Dart hd3 = phi3(m, hd);
			if (is_indexed<GMap3::Vertex>(m))
				copy_index<GMap3::Vertex>(m, hd, phi1(m, hd3));
			if (is_indexed<GMap3::Edge>(m))
				copy_index<GMap3::Edge>(m, hd, hd3);
			if (is_indexed<GMap3::Face>(m))
				copy_index<GMap3::Face>(m, hd, hd3);
			return true;
		});
	}

	return hole;
}

uint32 close(GMap3& m, bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		if (phi3(m, d) == d)
			fix_point_darts.push_back(d);

	for (Dart d : fix_point_darts)
	{
		if (phi3(m, d) == d)
		{
			GMap3::Volume h = close_hole(m, d, set_indices);
			foreach_dart_of_orbit(m, h, [&](Dart hd) -> bool {
				set_boundary(m, hd, true);
				return true;
			});
			++nb_holes;
		}
	}

	return nb_holes;
}

} // namespace cgogn
