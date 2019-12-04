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

#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/cells.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/dart.h>
#include <cgogn/core/functions/phi.h>
#include <cgogn/core/functions/mesh_info.h>

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

CMap2::Volume
add_pyramid(CMap2& m, uint32 size, bool set_indices)
{
	CMap1::Face first = add_face(static_cast<CMap1&>(m), 3u, false); // First triangle
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next triangles
	{
		CMap1::Face next = add_face(static_cast<CMap1&>(m), 3u, false);
		phi2_sew(m,phi_1(m,current), phi1(m,next.dart));
		current = next.dart;
	}
	phi2_sew(m,phi_1(m,current), phi1(m,first.dart)); // Finish the umbrella
	CMap2::Face base = close_hole(m,first.dart, false); // Add the base face
	
	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(m, vol, [&] (CMap2::Vertex v) -> bool
			{
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				set_index(m, CMap2::HalfEdge(phi2(m,e.dart)), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Face>(m))
		{
			foreach_incident_face(m, vol, [&] (CMap2::Face f) -> bool
			{
				set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

CMap2::Volume
add_prism(CMap2& m, uint32 size, bool set_indices)
{
	CMap1::Face first = add_face(static_cast<CMap1&>(m), 4u, false); // first quad
	Dart current = first.dart;
	for (uint32 i = 1u; i < size; ++i) // Next quads
	{
		CMap1::Face next = add_face(static_cast<CMap1&>(m), 4u, false);
		phi2_sew(m,phi_1(m,current), phi1(m,next.dart));
		current = next.dart;
	}
	phi2_sew(m,phi_1(m,current), phi1(m,first.dart)); // Finish the sides
	CMap2::Face base = close_hole(m,first.dart, false); // Add the base face
	close_hole(m,phi<11>(m,first.dart), false); // Add the top face
	
	CMap2::Volume vol(base.dart);

	if (set_indices)
	{
		if (is_indexed<CMap2::Vertex>(m))
		{
			foreach_incident_vertex(m, vol, [&] (CMap2::Vertex v) -> bool
			{
				set_index(m, v, new_index<CMap2::Vertex>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, CMap2::HalfEdge(e.dart), new_index<CMap2::HalfEdge>(m));
				set_index(m, CMap2::HalfEdge(phi2(m,e.dart)), new_index<CMap2::HalfEdge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			foreach_incident_edge(m, vol, [&] (CMap2::Edge e) -> bool
			{
				set_index(m, e, new_index<CMap2::Edge>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Face>(m))
		{
			foreach_incident_face(m, vol, [&] (CMap2::Face f) -> bool
			{
				set_index(m, f, new_index<CMap2::Face>(m));
				return true;
			});
		}
		if (is_indexed<CMap2::Volume>(m))
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}

	return vol;
}

///////////
// CMap3 //
///////////

CMap3::Face
cut_volume(CMap3& m, const std::vector<Dart>& path, bool set_indices)
{
	Dart f0 = add_face(static_cast<CMap1&>(m), path.size(), false).dart;
	Dart f1 = add_face(static_cast<CMap1&>(m), path.size(), false).dart;
	
	for(Dart d0: path)
	{
		Dart d1 = phi2(m,d0);
		phi2_unsew(m,d0);

		phi2_sew(m,d0, f0);
		phi2_sew(m,d1, f1);

		phi3_sew(m,f0, f1);

		f0 = phi_1(m,f0);
		f1 = phi1(m,f1);
	}

	if (set_indices)
	{
		if (is_indexed<CMap3::Vertex>(m))
		{
			foreach_dart_of_orbit(m,CMap3::Face(f0), [&](Dart d) -> bool
			{
				copy_index<CMap3::Vertex>(m,d, phi<21>(m,d));
				return true;
			});
		}
		if (is_indexed<CMap3::Edge>(m))
		{
			foreach_dart_of_orbit(m,CMap3::Face2(f0), [&](Dart d) -> bool
			{
				copy_index<CMap3::Edge>(m,d, phi2(m,d));
				copy_index<CMap3::Edge>(m,phi3(m,d), d);
				return true;
			});
		}
		if (is_indexed<CMap3::Face>(m))
		{
			set_index(m, CMap3::Face(f0), new_index<CMap3::Face>(m));
		}
	}

	return CMap3::Face(phi_1(m,f0));
}

CMap3::Volume close_hole(CMap3& m,Dart d, bool set_indices)
{
	cgogn_message_assert(phi3(m,d) == d, "CMap3: close hole called on a dart that is not a phi3 fix point");

	DartMarkerStore marker(m.mesh());
	DartMarkerStore hole_volume_marker(m.mesh());

	std::vector<Dart> visited_faces;
	visited_faces.reserve(1024u);

	visited_faces.push_back(d);
	foreach_dart_of_orbit(m,CMap3::Face2(d), [&] (Dart fd) -> bool { marker.mark(fd); return true; });
	
	uint32 count = 0u;

	for (uint32 i = 0u; i < visited_faces.size(); ++i)
	{
		const Dart it = visited_faces[i];
		Dart f = it;

		CMap3::Face2 hf = add_face(static_cast<CMap1&>(m), codegree(m, CMap3::Face(f)));
		foreach_dart_of_orbit(m,hf, [&] (Dart fd) -> bool { hole_volume_marker.mark(fd); return true; });
		
		++count;

		Dart bit = hf.dart;
		do
		{
			Dart e = phi3(m,phi2(m,f));
			bool found = false;
			do
			{
				if (phi3(m,e) == e)
				{
					found = true;
					if (!marker.is_marked(e))
					{
						visited_faces.push_back(e);
						foreach_dart_of_orbit(m,CMap3::Face2(e), [&] (Dart fd) -> bool { marker.mark(fd); return true; });
					}
				}
				else
				{
					if (hole_volume_marker.is_marked(e))
					{
						found = true;
						phi2_sew(m,e, bit);
					}
					else
						e = phi3(m,phi2(m,e));
				}
			} while(!found);

			phi3_sew(m,f, bit);
			bit = phi_1(m,bit);
			f = phi1(m,f);
		} while (f != it);
	}

	CMap3::Volume hole(phi3(m,d));

	if (set_indices)
	{
		foreach_dart_of_orbit(m,hole, [&] (Dart hd) -> bool
		{
			Dart hd3 = phi3(m,hd);
			if (is_indexed<CMap3::Vertex>(m))
				copy_index<CMap3::Vertex>(m,hd, phi1(m,hd3));
			if (is_indexed<CMap3::Edge>(m))
				copy_index<CMap3::Edge>(m,hd, hd3);
			if (is_indexed<CMap3::Face>(m))
				copy_index<CMap3::Face>(m,hd, hd3);
			return true;
		});
	}

	return hole;
}

uint32 close(CMap3& m,bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	foreach_dart(m,[&] (Dart d) -> bool
	{
		if (phi3(m,d) == d)
			fix_point_darts.push_back(d);
		return true;
	});

	for (Dart d : fix_point_darts)
	{
		if (phi3(m,d) == d)
		{
			CMap3::Volume h = close_hole(m,d, set_indices);
			foreach_dart_of_orbit(m,h, [&] (Dart hd) -> bool { set_boundary(m,hd, true); return true; });
			++nb_holes;
		}
	}

	return nb_holes;
}

void sew_volumes(CMap3& m,Dart d0, Dart d1)
{
	cgogn_message_assert(codegree(m, CMap3::Face(d0)) == codegree(m, CMap3::Face(d1)), "The faces to sew do not have the same codegree");
	Dart it0 = d0;
	Dart it1 = d1;
	do
	{
		cgogn_message_assert(phi3(m,it0) == it0 && phi3(m,it1) == it1, "The faces to sew are already sewn");
		phi3_sew(m,it0, it1);
		it0 = phi1(m,it0);
		it1 = phi_1(m,it1);
	} while (it0 != d0);
}

} // namespace cgogn
