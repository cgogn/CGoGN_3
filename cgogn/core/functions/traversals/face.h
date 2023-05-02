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

#ifndef CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_
#define CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_

// #include <cgogn/core/cgogn_core_export.h>

// #include <cgogn/core/utils/tuples.h>
// #include <cgogn/core/utils/type_traits.h>

namespace cgogn
{

template <typename MESH, typename CELL>
std::vector<typename mesh_traits<MESH>::Face> incident_faces(const MESH& m, CELL c)
{
	using Face = typename mesh_traits<MESH>::Face;
	std::vector<Face> faces;
	faces.reserve(32u);
	foreach_incident_face(m, c, [&](Face f) -> bool {
		faces.push_back(f);
		return true;
	});
	return faces;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_TRAVERSALS_FACE_H_
