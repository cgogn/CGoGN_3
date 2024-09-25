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

#ifndef CGOGN_IO_POINT_IMPORT_H_
#define CGOGN_IO_POINT_IMPORT_H_

#include <cgogn/io/cgogn_io_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

namespace cgogn
{

struct CMap0;

namespace io
{

using geometry::Vec3;

struct PointImportData
{
	uint32 nb_vertices_ = 0;

	std::vector<Vec3> vertex_position_;
	std::string vertex_position_attribute_name_ = "position";

	inline void reserve(uint32 nb_vertices)
	{
		vertex_position_.reserve(nb_vertices);
		nb_vertices_ = nb_vertices;
	}
};

void CGOGN_IO_EXPORT import_point_data(CMap0& m, const PointImportData& point_data);

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_POINT_IMPORT_H_
