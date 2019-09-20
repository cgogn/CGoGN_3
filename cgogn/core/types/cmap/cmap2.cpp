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

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/core/functions/mesh_ops/face.h>

namespace cgogn
{

Dart CMap2::close_hole(Dart d, bool set_indices)
{
	cgogn_message_assert(phi2(d) == d, "CMap2: close hole called on a dart that is not a phi2 fix point");

	Dart first = add_dart();	// First edge of the face that will fill the hole
	phi2_sew(d, first);			// 2-sew the new edge to the hole

	Dart d_next = d;			// Turn around the hole
	Dart d_phi1;				// to complete the face
	do
	{
		do
		{
			d_phi1 = phi1(d_next); // Search and put in d_next
			d_next = phi2(d_phi1); // the next dart of the hole
		} while (d_next != d_phi1 && d_phi1 != d);

		if (d_phi1 != d)
		{
			Dart next = add_dart();	// Add a vertex into the built face
			phi1_sew(first, next);
			phi2_sew(d_next, next);	// and 2-sew the face to the hole
		}
	} while (d_phi1 != d);

	if (set_indices)
	{
		Dart it = first;
		do
		{
			Dart it2 = phi2(it);
			if (is_embedded<Vertex>())
				copy_embedding<Vertex>(it, phi1(it2));
			if (is_embedded<Edge>())
				copy_embedding<Edge>(it, it2);
			if (is_embedded<Volume>())
				copy_embedding<Volume>(it, it2);
			it = phi1(it);
		} while (it != first);
	}

	return first;
}

uint32 CMap2::close(bool set_indices)
{
	uint32 nb_holes = 0u;

	std::vector<Dart> fix_point_darts;
	foreach_dart([&] (Dart d) -> bool
	{
		if (phi2(d) == d)
			fix_point_darts.push_back(d);
		return true;
	});

	for (Dart d : fix_point_darts)
	{
		if (phi2(d) == d)
		{
			Dart h = close_hole(d, set_indices);
			foreach_dart_of_orbit(CMap2::Face(h), [&] (Dart hd) -> bool { set_boundary(hd, true); return true; });
			++nb_holes;
		}
	}

	return nb_holes;
}

} // namespace cgogn
