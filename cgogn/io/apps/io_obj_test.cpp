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

#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/types/cell_marker.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/io/surface/obj.h>

using namespace cgogn;

using geometry::Vec3;

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

uint32 count_bdart(const CMap2& m)
{
	uint32 nb = 0u;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
		if (is_boundary(m, d))
			++nb;
	return nb;
}

int main()
{
	CMap2 map_pos;
	CMap2 map_tc;
	CMap2 map_no;

	auto filename = std::string(DEFAULT_MESH_PATH) + std::string("wavefront_obj/cube.obj");

	io::import_OBJ_tn(map_pos, map_tc, map_no, filename);

	std::cout << "nb darts pos: " << nb_darts(map_pos);
	std::cout << "   with in boundary : " << count_bdart(map_pos) << std::endl;
	std::cout << "nb darts tc: " << nb_darts(map_tc);
	std::cout << "   with in boundary : " << count_bdart(map_tc) << std::endl;
	std::cout << "nb darts no: " << nb_darts(map_no);
	std::cout << "   with in boundary : " << count_bdart(map_no) << std::endl;

	std::cout << "nb Faces pos: " << nb_cells<CMap2::Face>(map_pos) << std::endl;
	std::cout << "nb Faces tc: " << nb_cells<CMap2::Face>(map_tc) << std::endl;
	std::cout << "nb Faces no: " << nb_cells<CMap2::Face>(map_no) << std::endl;

	std::cout << "nb Vertices pos: " << nb_cells<CMap2::Vertex>(map_pos) << std::endl;
	std::cout << "nb Vertices tc: " << nb_cells<CMap2::Vertex>(map_tc) << std::endl;
	std::cout << "nb Vertices no: " << nb_cells<CMap2::Vertex>(map_no) << std::endl;

	std::cout << "nb Edge pos: " << nb_cells<CMap2::Edge>(map_pos) << std::endl;
	std::cout << "nb Edge tc: " << nb_cells<CMap2::Edge>(map_tc) << std::endl;
	std::cout << "nb Edge no: " << nb_cells<CMap2::Edge>(map_no) << std::endl;

	foreach_cell(map_pos, [&](CMap2::Edge e)->bool 
	{
		Dart d = e.dart_;
		Dart d2 = phi2(map_pos,d);
		if (is_boundary(map_pos, phi2(map_pos, d)) != is_boundary(map_tc, phi2(map_tc, d)))
			std::cout << "Edge of dart " << d.index_ << " is split in TC"<< std::endl;
		return true;
	});
}
