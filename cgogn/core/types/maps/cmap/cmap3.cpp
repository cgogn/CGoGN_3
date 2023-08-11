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
#include <cgogn/core/types/maps/cmap/cmap3.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/cell_marker.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

CMap3::Vertex cut_edge(CMap3& m, CMap3::Edge e, bool set_indices)
{
	Dart d0 = e.dart;
	Dart d23 = phi<2, 3>(m, d0);

	CMap3::Vertex v(cut_edge(static_cast<CMap2&>(m), CMap2::Edge(d0), false).dart);

	while (d23 != e.dart)
	{
		d0 = d23;
		d23 = phi<2, 3>(m, d23);

		cut_edge(static_cast<CMap2&>(m), CMap2::Edge(d0), false);

		const Dart d3 = phi3(m, d0);
		phi3_unsew(m, d0);

		phi3_sew(m, d0, phi1(m, d3));
		phi3_sew(m, d3, phi1(m, d0));
	}

	const Dart d3 = phi3(m, e.dart);
	phi3_unsew(m, e.dart);

	phi3_sew(m, e.dart, phi1(m, d3));
	phi3_sew(m, d3, phi1(m, e.dart));

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
			set_index(m, v, new_index<CMap3::Vertex>(m));
		if (is_indexed<CMap3::Vertex2>(m))
		{
			Dart d = v.dart;
			do
			{
				if (!is_boundary(m, d))
					set_index(m, CMap3::Vertex2(d), new_index<CMap3::Vertex2>(m));
				d = phi<2, 3>(m, d);
			} while (d != v.dart);
		}
		if (is_indexed<CMap3::Edge>(m))
		{
			set_index(m, CMap3::Edge(v.dart), new_index<CMap3::Edge>(m));
			set_index(m, e, index_of(m, e));
		}
		if (is_indexed<CMap3::Face2>(m))
		{
			Dart d = e.dart;
			do
			{
				// if (!is_boundary(m, d))
				copy_index<CMap3::Face2>(m, phi1(m, d), d);
				// if (!is_boundary(m, phi3(m, d)))
				copy_index<CMap3::Face2>(m, phi3(m, d), phi<1, 3>(m, d));
				d = phi<2, 3>(m, d);
			} while (d != e.dart);
		}
		if (is_indexed<CMap3::Face>(m))
		{
			Dart d = e.dart;
			do
			{
				copy_index<CMap3::Face>(m, phi1(m, d), d);
				copy_index<CMap3::Face>(m, phi3(m, d), d);
				// copy_index<CMap3::Face>(m, phi2(m, d), phi<1, 2>(m, d));
				d = phi<2, 3>(m, d);
			} while (d != e.dart);
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			foreach_dart_of_orbit(m, e, [&](Dart d) -> bool {
				if (is_boundary(m, d))
					return true;
				copy_index<CMap3::Volume>(m, phi1(m, d), d);
				copy_index<CMap3::Volume>(m, phi2(m, d), d);
				return true;
			});
		}
	}

	return v;
}

CMap3::Edge cut_face(CMap3& m, CMap3::Vertex v1, CMap3::Vertex v2, bool set_indices)
{
	Dart d = v1.dart;
	Dart e = v2.dart;

	Dart dd = phi<3, 1>(m, d);
	Dart ee = phi<3, 1>(m, e);

	cut_face(static_cast<CMap2&>(m), CMap2::Vertex(v1.dart), CMap2::Vertex(e), false);
	cut_face(static_cast<CMap2&>(m), CMap2::Vertex(dd), CMap2::Vertex(ee), false);

	phi3_sew(m, phi_1(m, v1.dart), phi_1(m, ee));
	phi3_sew(m, phi_1(m, dd), phi_1(m, e));

	CMap3::Edge edge(phi_1(m, e));

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
		{
			copy_index<CMap3::Vertex>(m, phi_1(m, e), v1.dart);
			copy_index<CMap3::Vertex>(m, phi_1(m, ee), v1.dart);
			copy_index<CMap3::Vertex>(m, phi_1(m, d), e);
			copy_index<CMap3::Vertex>(m, phi_1(m, dd), e);
		}
		if (is_indexed<CMap3::Vertex2>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CMap3::Vertex2>(m, phi_1(m, d), e);
				copy_index<CMap3::Vertex2>(m, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CMap3::Vertex2>(m, phi_1(m, dd), ee);
				copy_index<CMap3::Vertex2>(m, phi_1(m, ee), dd);
			}
		}
		if (is_indexed<CMap3::Edge>(m))
			set_index(m, CMap3::Edge(phi_1(m, v1.dart)), new_index<CMap3::Edge>(m));
		if (is_indexed<CMap3::Face2>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CMap3::Face2>(m, phi_1(m, d), d);
				set_index(m, CMap3::Face2(e), new_index<CMap3::Face2>(m));
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CMap3::Face2>(m, phi_1(m, dd), dd);
				set_index(m, CMap3::Face2(ee), new_index<CMap3::Face2>(m));
			}
		}
		if (is_indexed<CMap3::Face>(m))
		{
			copy_index<CMap3::Face>(m, phi_1(m, ee), d);
			copy_index<CMap3::Face>(m, phi_1(m, d), d);
			set_index(m, CMap3::Face(e), new_index<CMap3::Face>(m));
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CMap3::Volume>(m, phi_1(m, d), d);
				copy_index<CMap3::Volume>(m, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CMap3::Volume>(m, phi_1(m, dd), dd);
				copy_index<CMap3::Volume>(m, phi_1(m, ee), dd);
			}
		}
	}

	return edge;
}

CMap3::Face cut_volume(CMap3& m, const std::vector<Dart>& path, bool set_indices)
{
	Dart f0 = add_face(static_cast<CMap1&>(m), uint32(path.size()), false).dart;
	Dart f1 = add_face(static_cast<CMap1&>(m), uint32(path.size()), false).dart;

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
		if (is_indexed<CMap3::Vertex>(m))
		{
			foreach_dart_of_orbit(m, CMap3::Face(f0), [&](Dart d) -> bool {
				copy_index<CMap3::Vertex>(m, d, phi<2, 1>(m, d));
				return true;
			});
		}
		if (is_indexed<CMap3::Vertex2>(m))
		{
			foreach_dart_of_orbit(m, CMap3::Face(f0), [&](Dart d) -> bool {
				copy_index<CMap3::Vertex2>(m, d, phi<2, 1>(m, d));
				return true;
			});
		}
		if (is_indexed<CMap3::Edge>(m))
		{
			foreach_dart_of_orbit(m, CMap3::Face2(f0), [&](Dart d) -> bool {
				copy_index<CMap3::Edge>(m, d, phi2(m, d));
				copy_index<CMap3::Edge>(m, phi3(m, d), d);
				return true;
			});
		}
		if (is_indexed<CMap3::Face2>(m))
		{
			set_index(m, CMap3::Face2(f0), new_index<CMap3::Face2>(m));
			set_index(m, CMap3::Face2(phi3(m, f0)), new_index<CMap3::Face2>(m));
		}
		if (is_indexed<CMap3::Face>(m))
		{
			set_index(m, CMap3::Face(f0), new_index<CMap3::Face>(m));
		}
		if (is_indexed<CMap3::Volume>(m))
		{
			foreach_dart_of_orbit(m, CMap3::Face2(f0), [&](Dart d) -> bool {
				copy_index<CMap3::Volume>(m, d, phi2(m, d));
				return true;
			});
			set_index(m, CMap3::Volume(f1), new_index<CMap3::Volume>(m));
		}
	}

	return CMap3::Face(f0);
}

CMap3::Volume close_hole(CMap3& m, Dart d, bool set_indices)
{
	cgogn_message_assert(phi3(m, d) == d, "CMap3: close hole called on a dart that is not a phi3 fix point");

	DartMarkerStore marker(m);
	DartMarkerStore hole_volume_marker(m);

	std::vector<Dart> visited_faces;
	visited_faces.reserve(1024u);

	visited_faces.push_back(d);
	foreach_dart_of_orbit(m, CMap3::Face2(d), [&](Dart fd) -> bool {
		marker.mark(fd);
		return true;
	});

	uint32 count = 0u;

	for (uint32 i = 0u; i < uint32(visited_faces.size()); ++i)
	{
		const Dart it = visited_faces[i];
		Dart f = it;

		CMap3::Face2 hf = add_face(static_cast<CMap1&>(m), codegree(m, CMap3::Face(f)));
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
						foreach_dart_of_orbit(m, CMap3::Face2(e), [&](Dart fd) -> bool {
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

	CMap3::Volume hole(phi3(m, d));

	if (set_indices)
	{
		foreach_dart_of_orbit(m, hole, [&](Dart hd) -> bool {
			Dart hd3 = phi3(m, hd);
			if (is_indexed<CMap3::Vertex>(m))
				copy_index<CMap3::Vertex>(m, hd, phi1(m, hd3));
			if (is_indexed<CMap3::Edge>(m))
				copy_index<CMap3::Edge>(m, hd, hd3);
			if (is_indexed<CMap3::Face>(m))
				copy_index<CMap3::Face>(m, hd, hd3);
			return true;
		});
	}

	return hole;
}

uint32 close(CMap3& m, bool set_indices)
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
			CMap3::Volume h = close_hole(m, d, set_indices);
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
// Cells information
/*************************************************************************/

std::vector<uint32> hexahedra_vertex_indices(const CMap3& m, CMap3::Attribute<uint32>* vertex_id, CMap3::Volume v)
{
	using Vertex = CMap3::Vertex;

	Dart d1 = v.dart;
	Dart d2 = phi<2, 1, 1, 2>(m, d1);

	return {value<uint32>(m, vertex_id, Vertex(d1)),
			value<uint32>(m, vertex_id, Vertex(phi1(m, d2))),
			value<uint32>(m, vertex_id, Vertex(phi<1, 1>(m, d2))),
			value<uint32>(m, vertex_id, Vertex(phi_1(m, d1))),
			value<uint32>(m, vertex_id, Vertex(phi1(m, d1))),
			value<uint32>(m, vertex_id, Vertex(d2)),
			value<uint32>(m, vertex_id, Vertex(phi_1(m, d2))),
			value<uint32>(m, vertex_id, Vertex(phi<1, 1>(m, d1)))

	};
}

/*************************************************************************/
// Debugging helper functions
/*************************************************************************/

bool check_integrity(CMap3& m, bool verbose)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = true;
		relations &= phi3(m, d) != d && phi<3, 3>(m, d) == d && phi<3, 1, 3, 1>(m, d) == d;
		relations &= phi2(m, d) != d && phi<2, 2>(m, d) == d;
		relations &= phi<-1, 1>(m, d) == d && phi<1, -1>(m, d) == d;
		if (verbose && !relations)
			std::cerr << "Dart " << d << " has bad relations" << std::endl;

		bool boundary = is_boundary(m, d) == is_boundary(m, phi1(m, d)) &&
						is_boundary(m, d) == is_boundary(m, phi2(m, d)) &&
						(!is_boundary(m, d) || !is_boundary(m, phi3(m, d)));
		if (verbose && !boundary)
			std::cerr << "Dart " << d << " has bad boundary" << std::endl;

		result &= relations && boundary;
	}
	result &= check_indexing<CMap3::Vertex>(m);
	result &= check_indexing<CMap3::Vertex2>(m);
	result &= check_indexing<CMap3::HalfEdge>(m);
	result &= check_indexing<CMap3::Edge>(m);
	result &= check_indexing<CMap3::Edge2>(m);
	result &= check_indexing<CMap3::Face>(m);
	result &= check_indexing<CMap3::Face2>(m);
	result &= check_indexing<CMap3::Volume>(m);
	return result;
}

} // namespace cgogn
