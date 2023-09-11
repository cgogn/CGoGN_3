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

#ifndef CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_UTILS_H_
#define CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_UTILS_H_

#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/modeling/algos/subdivision_utils.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH, typename FUNC>
void cut_all_edges(MESH& m, const FUNC& on_edge_cut)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Given function should take a Vertex");

	CellCache<MESH> cache(m);
	cache.template build<Edge>();

	foreach_cell(cache, [&](Edge e) -> bool {
		on_edge_cut(cut_edge(m, e));
		return true;
	});
}

template <typename MESH, typename FUNC1, typename FUNC2>
void quadrangulate_all_faces(MESH& m, const FUNC1& on_edge_cut, const FUNC2& on_face_cut)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;
	static_assert(is_func_parameter_same<FUNC1, Vertex>::value, "Given function should take a Vertex");
	static_assert(is_func_parameter_same<FUNC2, Vertex>::value, "Given function should take a Vertex");

	CellCache<MESH> cache(m);
	cache.template build<Face>();

	CellMarker<MESH, Edge> em(m);
	foreach_cell(cache, [&](Face f) -> bool {
		foreach_incident_edge(m, f, [&](Edge ie) -> bool {
			if (!em.is_marked(ie))
			{
				em.mark(ie);
				cache.add(ie);
			}
			return true;
		});
		return true;
	});

	CellMarker<MESH, Vertex> vm(m);
	foreach_cell(cache, [&](Edge e) -> bool {
		Vertex v = cut_edge(m, e);
		vm.mark(v);
		on_edge_cut(v);
		return true;
	});

	foreach_cell(cache, [&](Face f) -> bool {
		on_face_cut(quadrangulate_face(m, f, vm));
		return true;
	});
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_SUBDIVISION_SURFACE_UTILS_H_
