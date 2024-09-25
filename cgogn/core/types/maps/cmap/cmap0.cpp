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

#include <cgogn/core/types/maps/cmap/cmap0.h>

namespace cgogn
{

/*************************************************************************/
// Operators
/*************************************************************************/

CMap0::Vertex add_vertex(CMap0& m, bool set_indices)
{
	using Vertex = CMap0::Vertex;

	Dart d = add_dart(m);
	Vertex v(d);
	if (set_indices)
	{
		if (is_indexed<Vertex>(m))
			set_index(m, v, new_index<Vertex>(m));
	}
	return v;
}

void remove_vertex(CMap0& m, CMap0::Vertex v)
{
	remove_dart(m, v.dart_);
}

} // namespace cgogn
