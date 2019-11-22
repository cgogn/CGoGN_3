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

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{

CMap3::Volume CMap3::close_hole(Dart d, bool set_indices)
{
	cgogn_message_assert(phi3(d) == d, "CMap3: close hole called on a dart that is not a phi3 fix point");

	DartMarkerStore marker(*this);
	DartMarkerStore hole_volume_marker(*this);

	std::vector<Dart> visited_faces;
	visited_faces.reserve(1024u);

	visited_faces.push_back(d);
	CMap2::foreach_dart_of_orbit(CMap2::Face(d), [&] (Dart fd) -> bool { marker.mark(fd); return true; });
	
	uint32 count = 0u;

	for (uint32 i = 0u; i < visited_faces.size(); ++i)
	{
		const Dart it = visited_faces[i];
		Dart f = it;

		CMap1::Face hf = add_face(static_cast<CMap1&>(*this), codegree(*this, CMap3::Face(f)));
		CMap1::foreach_dart_of_orbit(hf, [&] (Dart fd) -> bool { hole_volume_marker.mark(fd); return true; });
		++count;

		Dart bit = hf.dart;
		do
		{
			Dart e = phi3(phi2(f));
			bool found = false;
			do
			{
				if (phi3(e) == e)
				{
					found = true;
					if (!marker.is_marked(e))
					{
						visited_faces.push_back(e);
						CMap2::foreach_dart_of_orbit(CMap2::Face(e), [&] (Dart fd) -> bool { marker.mark(fd); return true; });
					}
				}
				else
				{
					if (hole_volume_marker.is_marked(e))
					{
						found = true;
						phi2_sew(e, bit);
					}
					else
						e = phi3(phi2(e));
				}
			} while(!found);

			phi3_sew(f, bit);
			bit = phi_1(bit);
			f = phi1(f);
		} while (f != it);
	}

	Volume hole(phi3(d));

	if (set_indices)
	{
		foreach_dart_of_orbit(hole, [&] (Dart hd) -> bool
		{
			Dart hd3 = phi3(hd);
			if (is_indexed<Vertex>())
				copy_index<Vertex>(hd, phi1(hd3));
			if (is_indexed<Edge>())
				copy_index<Edge>(hd, hd3);
			if (is_indexed<Face>())
				copy_index<Face>(hd, hd3);
			return true;
		});
	}

	return hole;
}

uint32 CMap3::close(bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	foreach_dart([&] (Dart d) -> bool
	{
		if (phi3(d) == d)
			fix_point_darts.push_back(d);
		return true;
	});

	for (Dart d : fix_point_darts)
	{
		if (phi3(d) == d)
		{
			Volume h = close_hole(d, set_indices);
			foreach_dart_of_orbit(h, [&] (Dart hd) -> bool { set_boundary(hd, true); return true; });
			++nb_holes;
		}
	}

	return nb_holes;
}

void CMap3::sew_volumes(Dart d0, Dart d1)
{
	cgogn_message_assert(codegree(*this, CMap3::Face(d0)) == codegree(*this, CMap3::Face(d1)), "The faces to sew do not have the same codegree");
	Dart it0 = d0;
	Dart it1 = d1;
	do
	{
		cgogn_message_assert(phi3(it0) == it0 && phi3(it1) == it1, "The faces to sew are already sewn");
		phi3_sew(it0, it1);
		it0 = phi1(it0);
		it1 = phi_1(it1);
	} while (it0 != d0);
}

} // namespace cgogn
