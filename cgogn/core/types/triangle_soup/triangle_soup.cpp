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

#include <cgogn/core/types/triangle_soup/triangle_soup.h>

namespace cgogn
{

/*************************************************************************/
// Clear mesh
/*************************************************************************/

void clear(TriangleSoup& ts, bool keep_attributes)
{
	for (TriangleSoup::AttributeContainer& container : ts.attribute_containers_)
	{
		container.clear_attributes();
		if (!keep_attributes)
			container.remove_attributes();
	}
}

/*************************************************************************/
// Copy mesh
/*************************************************************************/

void copy(TriangleSoup& dst, const TriangleSoup& src)
{
	clear(dst, false);
	for (uint32 i = 0; i < dst.attribute_containers_.size(); ++i)
		dst.attribute_containers_[i].copy(src.attribute_containers_[i]);
}

/*************************************************************************/
// Operators
/*************************************************************************/

TriangleSoup::Face add_face(TriangleSoup& ts, TriangleSoup::Vertex v1, TriangleSoup::Vertex v2, TriangleSoup::Vertex v3)
{
	TriangleSoup::Face f = new_index<TriangleSoup::Face>(ts);
	(*ts.face_incident_vertices_)[f] = {v1, v2, v3};
	return f;
}

} // namespace cgogn
