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

#ifndef CGOGN_CORE_CMAP_CMAP_BASE_CHECK_IND_HPP_
#define CGOGN_CORE_CMAP_CMAP_BASE_CHECK_IND_HPP_

#include <sstream>
#include <iostream>
#include <iomanip>


namespace cgogn
{

template <typename CELL, typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>* = nullptr>
bool check_indexing(MESH& m, bool verbose = true)
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");

	if (!is_indexed<CELL>(m))
		return true;

	bool result = true;

	auto counter = add_attribute<uint32, CELL>(m, "__cell_counter");
	counter->fill(0);

	foreach_cell(
		m,
		[&](CELL c) -> bool {
			const uint32 index = index_of(m, c);

			++(*counter)[index];

			bool valid_index = index != INVALID_INDEX;
			if (verbose && !valid_index)
				std::cerr << "Cell " << c << " (" << cell_name<CELL>(m) << ") has invalid index" << std::endl;

			bool all_darts_same_index = true;
			foreach_dart_of_orbit(m, c, [&](Dart d) -> bool {
				const uint32 index_d = index_of(m, CELL(d));
				if (index_d != index)
				{
					if (verbose)
						std::cerr << "Cell " << c << " (" << cell_name<CELL>(m) << ") has darts with different indices"
								  << std::endl;
					all_darts_same_index = false;
				}
				return true;
			});

			result &= valid_index && all_darts_same_index;
			return true;
		},
		CMapBase::TraversalPolicy::DART_MARKING);

	// check that all lines of the attribute container are used
	for (uint32 i = m.attribute_containers_[CELL::ORBIT].first_index(),
				end = m.attribute_containers_[CELL::ORBIT].last_index();
		 i != end; i = m.attribute_containers_[CELL::ORBIT].next_index(i))
	{
		if ((*counter)[i] == 0)
		{
			if (verbose)
				std::cerr << "Cell index " << i << " is not used in container " << cell_name<CELL>(m) << std::endl;
			result = false;
		}
		else
		{
			if ((*counter)[i] >= 2ul)
			{
				if (verbose)
					std::cerr << "Multiple cells with same index " << i << " in container " << cell_name<CELL>(m)
							  << std::endl;
				result = false;
			}
		}
	}

	remove_attribute<CELL>(m, counter);

	return result;
}




void clear(CMapBase& m, bool keep_attributes = true);


template <typename MESH, typename std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>>* = nullptr>
void copy(MESH& dst, const MESH& src)
{
	clear(dst, false);
	for (uint32 orbit = 0; orbit < NB_ORBITS; ++orbit)
	{
		if (src.cells_indices_[orbit] != nullptr)
			init_cells_indexing(dst, Orbit(orbit));
	}
	dst.darts_.copy(src.darts_);
	for (uint32 i = 0; i < NB_ORBITS; ++i)
		dst.attribute_containers_[i].copy(src.attribute_containers_[i]);
	dst.boundary_marker_ = dst.darts_.get_mark_attribute();
	dst.boundary_marker_->copy(*src.boundary_marker_);
}


template <typename MESH, typename CELL>
auto is_incident_to_boundary(const MESH& m, CELL c) -> std::enable_if_t<std::is_convertible_v<MESH&, CMapBase&>, bool>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	bool result = false;
	foreach_dart_of_orbit(m, c, [&m, &result](Dart d) -> bool {
		if (is_boundary(m, d))
		{
			result = true;
			return false;
		}
		return true;
	});
	return result;
}


} // namespace cgogn


#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_
