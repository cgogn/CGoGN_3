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
#include <cgogn/core/types/maps/gmap/gmap1.h>
#include <cgogn/core/functions/mesh_info.h>

namespace cgogn
{

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
	beta0_unsew(m, e0);

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
	Dart ee = beta0(m, e.dart);
	Dart d2 = beta1(m, ee);

	beta1_unsew(m, d1);
	beta1_unsew(m, d2);

	beta1_sew(m, d1, d2);

	remove_edge(m, e);

	Vertex v(d2);

	if (set_indices)
		copy_index<Vertex>(m, d2, d1);

	return v;
}

} // namespace cgogn
                    
