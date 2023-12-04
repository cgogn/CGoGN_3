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

#ifndef CGOGN_IO_SURFACE_IMPORT_H_
#define CGOGN_IO_SURFACE_IMPORT_H_

#include <cgogn/io/cgogn_io_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <vector>

namespace cgogn
{

struct CMap2;
struct GMap2;
struct IncidenceGraph;
struct TriangleSoup;

namespace io
{


using geometry::Vec3;
using geometry::Vec2;

template <typename VEC>
struct SurfaceImportDataTGen
{
	using VEC_TYPE = VEC;

	uint32 nb_vertices_ = 0;
	uint32 nb_faces_ = 0;

	std::vector<VEC> vertex_position_;
	std::string vertex_position_attribute_name_ = "position";

	std::vector<uint32> faces_nb_vertices_;
	std::vector<uint32> faces_vertex_indices_;

	std::vector<uint32> vertex_id_after_import_;

	inline void reserve(uint32 nb_vertices, uint32 nb_faces)
	{
		nb_vertices_ = nb_vertices;
		nb_faces_ = nb_faces;
		vertex_position_.reserve(nb_vertices);
		faces_nb_vertices_.reserve(nb_faces);
		faces_vertex_indices_.reserve(nb_faces * 4u);
		vertex_id_after_import_.reserve(nb_vertices);
	}
};

using SurfaceImportData = SurfaceImportDataTGen<Vec3>;
using SurfaceImportData2D = SurfaceImportDataTGen<Vec2>;


void CGOGN_IO_EXPORT import_surface_data(CMap2& m, SurfaceImportData& surface_data, bool reconstruct_phi2 = true);
void CGOGN_IO_EXPORT import_surface_data(GMap2& m, SurfaceImportData& surface_data, bool reconstruct_phi2 = true);
void CGOGN_IO_EXPORT import_surface_data(IncidenceGraph& m, SurfaceImportData& surface_data);
void CGOGN_IO_EXPORT import_surface_data(TriangleSoup& m, SurfaceImportData& surface_data);

void CGOGN_IO_EXPORT import_surface_data(CMap2& m, SurfaceImportData2D& surface_data, bool reconstruct_phi2 = true);
void CGOGN_IO_EXPORT import_surface_data(GMap2& m, SurfaceImportData2D& surface_data, bool reconstruct_phi2 = true);

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_IMPORT_H_
