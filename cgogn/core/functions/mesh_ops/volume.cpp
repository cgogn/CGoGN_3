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

#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/cmap/orbit_traversal.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// add_pyramid(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

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
			foreach_incident_vertex(m, vol, [&](CMap2::Vertex v) -> bool {
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(m, vol, [&](CMap2::Edge e) -> bool {
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				set_index(m, CMap2::HalfEdge(phi2(m, e.dart)), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(m, vol, [&](CMap2::Edge e) -> bool {
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Face>(m))
		{
			foreach_incident_face(m, vol, [&](CMap2::Face f) -> bool {
				set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// add_prism(MESH& m, uint32 size, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap2 //
///////////

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
	close_hole(m, phi<11>(m, first.dart), false);		 // Add the top face

	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(m, vol, [&](CMap2::Vertex v) -> bool {
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(m, vol, [&](CMap2::Edge e) -> bool {
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				set_index(m, CMap2::HalfEdge(phi2(m, e.dart)), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(m, vol, [&](CMap2::Edge e) -> bool {
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Face>(m))
		{
			foreach_incident_face(m, vol, [&](CMap2::Face f) -> bool {
				set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Face
// cut_volume(MESH& m, const std::vector<Dart>& path, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap3 //
///////////

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
				copy_index<CMap3::Vertex>(m, d, phi<21>(m, d));
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

//////////
// CPH3 //
//////////

CPH3::CMAP::Face cut_volume(CPH3& m, const std::vector<Dart>& path, bool set_indices)
{
	CPH3::CMAP& map = static_cast<CPH3::CMAP&>(m);

	uint32 vid = m.refinement_face_id(path);
	uint32 vlevel = m.volume_level(path[0]);

	CPH3::CMAP::Face result = cut_volume(map, path, false);

	Dart f0 = result.dart;
	Dart f1 = phi3(m, f0);

	foreach_dart_of_orbit(m, result, [&](Dart d) -> bool {
		m.set_edge_id(d, m.edge_id(phi2(m, d)));
		m.set_face_id(d, vid);
		m.set_dart_level(d, m.current_level_);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CPH3::CMAP::Vertex>(m))
		{
			foreach_dart_of_orbit(m, CPH3::CMAP::Face(f0), [&](Dart d) -> bool {
				copy_index<CPH3::CMAP::Vertex>(m, d, phi<21>(m, d));
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Edge>(m))
		{
			foreach_dart_of_orbit(m, CPH3::CMAP::Face2(f0), [&](Dart d) -> bool {
				copy_index<CPH3::CMAP::Edge>(map, d, phi2(m, d));
				copy_index<CPH3::CMAP::Edge>(map, phi3(m, d), d);
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Face>(m))
			set_index(m, CPH3::CMAP::Face(f0), new_index<CPH3::CMAP::Face>(m));
		if (is_indexed<CPH3::CMAP::Volume>(m))
		{
			if (vlevel == m.current_level_)
			{
				foreach_dart_of_orbit(m, CPH3::CMAP::Face2(f0), [&](Dart d) -> bool {
					copy_index<CPH3::CMAP::Volume>(map, d, phi2(m, d));
					return true;
				});
				set_index(m, CPH3::CMAP::Volume(f1), new_index<CPH3::CMAP::Volume>(m));
			}
			else
			{
				uint32 ved1 = new_index<CPH3::CMAP::Volume>(m);
				uint32 ved2 = new_index<CPH3::CMAP::Volume>(m);
				foreach_dart_of_orbit(m, CPH3::CMAP::Volume(f0), [&](Dart d) -> bool {
					if (m.dart_level(d) == m.current_level_)
						set_index<CPH3::CMAP::Volume>(m, d, ved1);
					return true;
				});
				foreach_dart_of_orbit(m, CPH3::CMAP::Volume(f1), [&](Dart d) -> bool {
					if (m.dart_level(d) == m.current_level_)
						set_index<CPH3::CMAP::Volume>(m, d, ved2);
					return true;
				});
			}
		}
	}

	return result;
}

/*****************************************************************************/

// template <typename MESH>
// typename mesh_traits<MESH>::Volume
// close_hole(MESH& m, Dart d, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap3 //
///////////

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

/*****************************************************************************/

// template <typename MESH>
// uint32
// close(MESH& m, bool set_indices = true);

/*****************************************************************************/

///////////
// CMap3 //
///////////

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

} // namespace cgogn
