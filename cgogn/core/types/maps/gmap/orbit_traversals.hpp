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

#ifndef CGOGN_CORE_TYPES_GMAP_ORBIT_TRAVERSAL_H_
#define CGOGN_CORE_TYPES_GMAP_ORBIT_TRAVERSAL_H_

#include <cgogn/core/types/maps/gmap/beta.h>
#include <cgogn/core/types/maps/dart_marker.h>

namespace cgogn
{

/*****************************************************************************/
// orbits traversals
/*****************************************************************************/

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 0)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (f(d))
		f(beta0(m, d));
}


template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA1(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if (f(d))
		f(beta1(m, d));
}


template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension>=1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it =d;
	do
	{
		if (!f(it))
			return;
		it = beta0(m, it);
		if (!f(it))
			return;
		it = beta1(m, it);
	} while ((it != d));
}

//template <typename MESH, typename BETA, typename FUNC>
//auto foreach_dart_of_BETAs(const MESH& m, Dart d, const BETA& betas, const FUNC& f)
//	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
//{
//	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
//	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
//
//	Dart it = d;
//	do
//	{
//		if (!f(it))
//			return;
//		it = betas(m, it);
//	} while (it != d);
//}


template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA01(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	Dart it = d;
	do
	{
		if (!f(it))
			return;

		it = beta<0,1>(m, it);
	} while (it != d);
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA21(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 1)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		it = beta<2,1>(m, it);
	} while ((it != d));
}

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA1_BETA2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 2)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		it = beta1(m, it);
		if (!f(it))
			return;
		it = beta2(m, it);
	} while ((it != d));
}


template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 2)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	if (!f(it))
		return;
	it = beta0(m, it);

	if (!f(it))
		return;
	it = beta2(m, it);

	if (!f(it))
		return;
	it = beta0(m, it);

	f(it);
}

// DIM 3

template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA2_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension == 3)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	Dart it = d;
	do
	{
		if (!f(it))
			return;
		if (!f(beta0(m,it)))
			return;
		it = beta2(m, it);
		if (!f(it))
			return;
		if (!f(beta0(m,it)))
			return;
		it = beta3(m, it);
	} while ((it != d));
}


template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1_BETA2(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension >= 2)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");


		// NON OPTIMIZED (BUT SO SIMPLE TO CODE!)

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		if (!f(curr_dart))
			break;

		const Dart d0 = beta0(m, curr_dart);
		const Dart d1 = beta1(m, curr_dart);
		const Dart d2 = beta2(m, curr_dart);

		if (!marker.is_marked(d0))
			marker.mark(d0);
		if (!marker.is_marked(d1))
			marker.mark(d1);
		if (!marker.is_marked(d2))
			marker.mark(d2);
	}
	//DartMarkerStore<MESH> marker(m);

	//std::vector<Dart> visited_faces;
	//visited_faces.reserve(32);
	//visited_faces.push_back(d); // Start with the face of d

	//// For every face added to the list
	//for (uint32 i = 0; i < uint32(visited_faces.size()); ++i)
	//{
	//	const Dart e = visited_faces[i];
	//	if (!marker.is_marked(e)) // Face has not been visited yet
	//	{
	//		foreach_dart_of_orbit(m, GMap1::Face(d), [&](Dart fd) -> bool 
	//		{
	//			if (!f(fd))
	//				return false;
	//			marker.mark(fd);				  // Mark
	//			const Dart adj = beta2(m, fd);	  // Get adjacent face
	//			if (!marker.is_marked(adj))		  // no need to test FP because of marker
	//				visited_faces.push_back(adj); // Add it
	//			return true;
	//		});
	//	}
	//}
}


template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&> && (mesh_traits<MESH>::dimension == 3)>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");
	foreach_dart_of_BETA0_BETA1(m, d, [&](Dart fd) -> bool {
		if (f(fd))
		{
			Dart e = beta3(m, fd);
			if (e != fd)
				return f(e);
		}
		return false;
	});
}



template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA1_BETA2_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		if (!f(curr_dart))
			break;

		const Dart d1 = beta1(m, curr_dart);
		const Dart d2 = beta2(m, curr_dart);
		const Dart d3 = beta3(m, curr_dart);

		if (!marker.is_marked(d1))
			marker.mark(d1);
		if (!marker.is_marked(d2))
			marker.mark(d2);
		if (!marker.is_marked(d3))
			marker.mark(d3);
	}
}


// TODO optimized version
template <typename MESH, typename FUNC>
auto foreach_dart_of_BETA0_BETA1_BETA2_BETA3(const MESH& m, Dart d, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	// NON OPTIMIZED (BUT SO SIMPLE TO CODE!)

	DartMarkerStore<MESH> marker(m);
	const std::vector<Dart>& marked_darts = marker.marked_darts();

	marker.mark(d);
	for (uint32 i = 0; i < uint32(marked_darts.size()); ++i)
	{
		const Dart curr_dart = marked_darts[i];
		if (!f(curr_dart))
			break;

		const Dart d0 = beta0(m, curr_dart);
		const Dart d1 = beta1(m, curr_dart);
		const Dart d2 = beta2(m, curr_dart);
		const Dart d3 = beta3(m, curr_dart);

		if (!marker.is_marked(d0))
			marker.mark(d0);
		if (!marker.is_marked(d1))
			marker.mark(d1);
		if (!marker.is_marked(d2))
			marker.mark(d2);
		if (!marker.is_marked(d3))
			marker.mark(d3);
	}
}


template <typename MESH, typename CELL, typename FUNC>
auto foreach_dart_of_orbit(const MESH& m, CELL c, const FUNC& f)
	-> std::enable_if_t<std::is_convertible_v<MESH&, GMapBase&>>
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Dart>::value, "Given function should take a Dart as parameter");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	static const Orbit orbit = CELL::ORBIT;

	if constexpr (orbit == DART)
	{
		unused_parameters(m);
		f(c.dart);
		return;
	}
	if constexpr (orbit == BETA0)
	{
		foreach_dart_of_BETA0(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA1)
	{
		foreach_dart_of_BETA1(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1)
	{
		foreach_dart_of_BETA0_BETA1(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA2)
	{
		foreach_dart_of_BETA0_BETA2(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA1_BETA2)
	{
		foreach_dart_of_BETA1_BETA2(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1_BETA2)
	{
		foreach_dart_of_BETA0_BETA1_BETA2(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1_BETA3)
	{
		foreach_dart_of_BETA0_BETA1_BETA3(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA2_BETA3)
	{
		foreach_dart_of_BETA0_BETA2_BETA3(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA1_BETA2_BETA3)
	{
		foreach_dart_of_BETA1_BETA2_BETA3(m, c.dart, f);
		return;
	}
	if constexpr (orbit == BETA0_BETA1_BETA2_BETA3)
	{
		foreach_dart_of_BETA0_BETA1_BETA2_BETA3(m, c.dart, f);
		return;
	}
}

} // namespace cgogn

#endif // CGOGN_CORE_GMAP_ORBIT_TRAVERSAL_H_
