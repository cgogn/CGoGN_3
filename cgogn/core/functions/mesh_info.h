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

#ifndef CGOGN_CORE_FUNCTIONS_MESH_INFO_H_
#define CGOGN_CORE_FUNCTIONS_MESH_INFO_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/cmap/cmap3.h>

#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/volume.h>

namespace cgogn
{

/*****************************************************************************/

// template <typename CELL, typename MESH>
// std::string cell_name(const MESH& m);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename MESH>
std::string cell_name(const MESH& m)
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	return mesh_traits<MESH>::cell_names[tuple_type_index<CELL, typename mesh_traits<MESH>::Cells>::value];
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// uint32 nb_cells(const MESH& m);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename CELL, typename MESH>
uint32 nb_cells(const MESH& m)
{
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	uint32 result = 0;
	foreach_cell(m, [&](CELL) -> bool {
		++result;
		return true;
	});
	return result;
}

/*****************************************************************************/

// template <typename MESH>
// bool is_simplicial(const MESH& m);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH>
bool is_simplicial(const MESH& m)
{
	bool res = true;
	if constexpr (mesh_traits<MESH>::dimension == 2)
	{
		foreach_cell(m, [&](typename mesh_traits<MESH>::Face f) -> bool {
			res = incident_vertices(m, f).size() == 3;
			return res;
		});
	}
	if constexpr (mesh_traits<MESH>::dimension == 3)
	{
		foreach_cell(m, [&](typename mesh_traits<MESH>::Volume v) -> bool {
			res = incident_vertices(m, v).size() == 4;
			return res;
		});
	}
	return res;
}

/*****************************************************************************/

// template <typename CELL, typename MESH>
// void check_indexing(MESH& m);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

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

/*****************************************************************************/

// template <typename MESH>
// void check_integrity(MESH& m);

/*****************************************************************************/

///////////
// CMap3 //
///////////

inline bool check_integrity(CMap3& m, bool verbose = true)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = true;
		relations &= phi3(m, d) != d && phi<33>(m, d) == d && phi<3131>(m, d) == d;
		relations &= phi2(m, d) != d && phi<22>(m, d) == d;
		relations &= phi1(m, phi_1(m, d)) == d && phi_1(m, phi1(m, d)) == d;
		if (verbose && !relations)
			std::cerr << "Dart " << d << " has bad relations" << std::endl;

		bool boundary = is_boundary(m, d) == is_boundary(m, phi1(m, d)) &&
						is_boundary(m, d) == is_boundary(m, phi2(m, d)) &&
						(!is_boundary(m, d) || !is_boundary(m, phi3(m, d)));
		if (verbose && !boundary)
			std::cerr << "Dart " << d << " has bad boundary" << std::endl;

		result &= relations && boundary;
	}
	result &= check_indexing<CMap3::Vertex>(m);
	result &= check_indexing<CMap3::Vertex2>(m);
	result &= check_indexing<CMap3::HalfEdge>(m);
	result &= check_indexing<CMap3::Edge>(m);
	result &= check_indexing<CMap3::Edge2>(m);
	result &= check_indexing<CMap3::Face>(m);
	result &= check_indexing<CMap3::Face2>(m);
	result &= check_indexing<CMap3::Volume>(m);
	return result;
}

///////////
// CMap2 //
///////////

inline bool check_integrity(CMap2& m, bool verbose = true)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = true;
		relations &= phi2(m, d) != d && phi<22>(m, d) == d;
		relations &= phi1(m, phi_1(m, d)) == d && phi_1(m, phi1(m, d)) == d;
		if (verbose && !relations)
		{
			std::cerr << "Dart " << d << " has bad relations" << std::endl;
			if (phi2(m, d) == d)
				std::cerr << "  phi2 fixed point" << std::endl;
			if (phi<22>(m, d) != d)
				std::cerr << "  phi2 not involution" << std::endl;
		}

		bool boundary =
			is_boundary(m, d) == is_boundary(m, phi1(m, d)) && (!is_boundary(m, d) || !is_boundary(m, phi2(m, d)));
		if (verbose && !boundary)
			std::cerr << "Dart " << d << " has bad boundary" << std::endl;

		result &= relations && boundary;
	}
	result &= check_indexing<CMap2::Vertex>(m);
	result &= check_indexing<CMap2::HalfEdge>(m);
	result &= check_indexing<CMap2::Edge>(m);
	result &= check_indexing<CMap2::Face>(m);
	result &= check_indexing<CMap2::Volume>(m);
	return result;
}

///////////
// CMap1 //
///////////

inline bool check_integrity(CMap1& m, bool verbose = true)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = phi1(m, phi_1(m, d)) == d && phi_1(m, phi1(m, d)) == d;
		if (verbose && !relations)
			std::cerr << "Dart " << d << " has bad relations" << std::endl;

		result &= relations;
	}
	result &= check_indexing<CMap1::Vertex>(m);
	result &= check_indexing<CMap1::Edge>(m);
	result &= check_indexing<CMap1::Face>(m);
	return result;
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// uint32 degree(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH>
uint32 degree(const MESH& m, typename mesh_traits<MESH>::Vertex v)
{
	uint32 result = 0;
	foreach_incident_edge(m, v, [&](typename mesh_traits<MESH>::Edge) -> bool {
		++result;
		return true;
	});
	return result;
}

template <typename MESH>
uint32 degree(const MESH& m, typename mesh_traits<MESH>::Edge e)
{
	uint32 result = 0;
	foreach_incident_face(m, e, [&](typename mesh_traits<MESH>::Face) -> bool {
		++result;
		return true;
	});
	return result;
}

template <typename MESH>
uint32 degree(const MESH& m, typename mesh_traits<MESH>::Face f)
{
	uint32 result = 0;
	foreach_incident_volume(m, f, [&](typename mesh_traits<MESH>::Volume) -> bool {
		++result;
		return true;
	});
	return result;
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// uint32 codegree(const MESH& m, CELL c);

/*****************************************************************************/

/////////////
// GENERIC //
/////////////

template <typename MESH>
uint32 codegree(const MESH& m, typename mesh_traits<MESH>::Edge e)
{
	uint32 result = 0;
	foreach_incident_vertex(m, e, [&](typename mesh_traits<MESH>::Vertex) -> bool {
		++result;
		return true;
	});
	return result;
}

template <typename MESH>
uint32 codegree(const MESH& m, typename mesh_traits<MESH>::Face f)
{
	uint32 result = 0;
	foreach_incident_edge(m, f, [&](typename mesh_traits<MESH>::Edge) -> bool {
		++result;
		return true;
	});
	return result;
}

template <typename MESH>
uint32 codegree(const MESH& m, typename mesh_traits<MESH>::Volume v)
{
	uint32 result = 0;
	foreach_incident_face(m, v, [&](typename mesh_traits<MESH>::Face) -> bool {
		++result;
		return true;
	});
	return result;
}

/*****************************************************************************/

// template <typename MESH, typename CELL>
// bool is_incident_to_boundary(const MESH& m, CELL c);

/*****************************************************************************/

//////////////
// CMapBase //
//////////////

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

/*****************************************************************************/

// template <typename MESH>
// bool edge_can_collapse(const MESH& m, typename mesh_traits<MESH>::Edge e);

/*****************************************************************************/

///////////
// CMap2 //
///////////

inline bool edge_can_collapse(const CMap2& m, CMap2::Edge e)
{
	using Vertex = CMap2::Vertex;
	using Face = CMap2::Face;

	auto vertices = incident_vertices(m, e);

	if (is_incident_to_boundary(m, vertices[0]) || is_incident_to_boundary(m, vertices[1]))
		return false;

	uint32 val_v1 = degree(m, vertices[0]);
	uint32 val_v2 = degree(m, vertices[1]);

	if (val_v1 + val_v2 < 8 || val_v1 + val_v2 > 14)
		return false;

	Dart e1 = e.dart;
	Dart e2 = phi2(m, e.dart);
	if (codegree(m, Face(e1)) == 3)
	{
		if (degree(m, Vertex(phi_1(m, e1))) < 4)
			return false;
	}
	if (codegree(m, Face(e2)) == 3)
	{
		if (degree(m, Vertex(phi_1(m, e2))) < 4)
			return false;
	}

	auto next_edge = [&m](Dart d) { return phi2(m, phi_1(m, d)); };

	// Check vertex sharing condition
	std::vector<uint32> vn1;
	Dart it = next_edge(next_edge(e1));
	Dart end = phi1(m, e2);
	do
	{
		vn1.push_back(index_of(m, Vertex(phi1(m, it))));
		it = next_edge(it);
	} while (it != end);
	it = next_edge(next_edge(e2));
	end = phi1(m, e1);
	do
	{
		auto vn1it = std::find(vn1.begin(), vn1.end(), index_of(m, Vertex(phi1(m, it))));
		if (vn1it != vn1.end())
			return false;
		it = next_edge(it);
	} while (it != end);

	return true;
}

/*****************************************************************************/

// template <typename MESH>
// bool edge_can_flip(const MESH& m, typename mesh_traits<MESH>::Edge e);

/*****************************************************************************/

///////////
// CMap2 //
///////////

inline bool edge_can_flip(const CMap2& m, CMap2::Edge e)
{
	if (is_incident_to_boundary(m, e))
		return false;

	Dart e1 = e.dart;
	Dart e2 = phi2(m, e1);

	auto next_edge = [&m](Dart d) { return phi2(m, phi_1(m, d)); };

	if (codegree(m, CMap2::Face(e1)) == 3 && codegree(m, CMap2::Face(e2)) == 3)
	{
		uint32 idxv2 = index_of(m, CMap2::Vertex(phi_1(m, e2)));
		Dart d = phi_1(m, e1);
		Dart it = d;
		do
		{
			if (index_of(m, CMap2::Vertex(phi1(m, it))) == idxv2)
				return false;
			it = next_edge(it);
		} while (it != d);
	}

	return true;
}

} // namespace cgogn

#endif // CGOGN_CORE_FUNCTIONS_MESH_INFO_H_
