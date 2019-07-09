/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps *
 * Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published by the
 ** Free Software Foundation; either version 2.1 of the License, or (at your
 ** option) any later version.
 **
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details. *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License *
 * along with this library; if not, write to the Free Software Foundation, *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA. *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/ * Contact information: cgogn@unistra.fr
 **
 *                                                                              *
 *******************************************************************************/

#include <cgogn/core/cmap/cmap2.h>
#include <cgogn/geometry/algos/bounding_box.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/io/map_import.h>
#include <cgogn/rendering/shaders/vbo.h>

#ifndef __INSTANCIED_TEMPLATES__
#define EXTERN extern
#else
#define EXTERN
#endif

namespace cgogn {
using Vec3 = Eigen::Vector3d;
using VA = CMap2::VertexAttribute<Vec3>;
EXTERN template void io::import_surface<Vec3>(CMap2&, const std::string&);
EXTERN template void geometry::compute_normal(const CMap2&, const VA&, VA&);
EXTERN template void geometry::compute_AABB(const VA&, AABB<Vec3>&);
EXTERN template void rendering::update_vbo(const VA&, VBO*);
}  // namespace cgogn
