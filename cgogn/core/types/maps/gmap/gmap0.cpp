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

#include <cgogn/core/types/maps/gmap/gmap0.h>

#include <cgogn/core/types/cell_marker.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

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
	Dart ee = beta0(m, e.dart);
	remove_dart(m, e.dart);
	remove_dart(m, ee);
}

} // namespace cgogn
