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

#ifndef CGOGN_MODELING_ALGOS_VOLUME_UTILS_H_
#define CGOGN_MODELING_ALGOS_VOLUME_UTILS_H_

#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

using geometry::Vec3;

//////////
// CMap //
//////////

void extract_volume_surface(CMap3& m3, CMap3::Attribute<Vec3>* m3_vertex_position, CMap2& m2,
							CMap2::Attribute<Vec3>* m2_vertex_position,
							CMap2::Attribute<CMap3::Vertex>* m2_vertex_m3_vertex = nullptr,
							CMap3::Attribute<CMap2::Vertex>* m3_vertex_m2_vertex = nullptr);


//void extract_volume_surface(GMap3& m3, GMap3::Attribute<Vec3>* m3_vertex_position, GMap2& m2,
//							GMap2::Attribute<Vec3>* m2_vertex_position,
//							GMap2::Attribute<GMap3::Vertex>* m2_vertex_m3_vertex = nullptr,
//							GMap3::Attribute<GMap2::Vertex>* m3_vertex_m2_vertex = nullptr);

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_VOLUME_UTILS_H_
