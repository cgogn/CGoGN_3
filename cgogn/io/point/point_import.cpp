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

#include <cgogn/io/point/point_import.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/core/types/cmap/cmap_ops.h>


#include <algorithm>
#include <set>
#include <vector>

namespace cgogn
{

namespace io
{

void import_point_data(CMap0& m, const PointImportData& point_data)
{
	using Vertex = CMap0::Vertex;

	auto position = get_or_add_attribute<geometry::Vec3, Vertex>(m, point_data.vertex_position_attribute_name_);

	for (uint32 i = 0u; i < point_data.nb_vertices_; ++i)
	{
		uint32 vertex_id = new_index<Vertex>(m);
		Vertex v = add_vertex(m);
		set_index(m, v, vertex_id);
		(*position)[vertex_id] = point_data.vertex_position_[i];
	}
}

} // namespace io

} // namespace cgogn
