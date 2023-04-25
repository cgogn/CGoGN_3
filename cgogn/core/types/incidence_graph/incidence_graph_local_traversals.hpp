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

#ifndef CGOGN_CORE_INCIDENCE_GRAPH_LOCAL_TRAVERSAL_H_
#define CGOGN_CORE_INCIDENCE_GRAPH_LOCAL_TRAVERSAL_H_

#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/functions/traversals/face.h>

namespace cgogn
{
// VERTEX

template <typename CELL, typename FUNC>
auto foreach_incident_vertex(const IncidenceGraph& ig, CELL c, const FUNC& func)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	static_assert(is_in_tuple<CELL, mesh_traits<IncidenceGraph>::Cells>::value, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, Edge>)
	{
		const std::pair<Vertex, Vertex>& evs = (*ig.edge_incident_vertices_)[c.index_];
		if (func(evs.first))
			func(evs.second);
	}
	else if constexpr (std::is_same_v<CELL, Face>)
	{
		// strong precondition: edges are sorted in the face & edges dirs are computed
		const std::vector<Edge>& edges = (*ig.face_incident_edges_)[c.index_];
		const std::vector<uint8>& edges_dir = (*ig.face_incident_edges_dir_)[c.index_];
		for (uint32 i = 0, end = edges.size(); i < end - 1; ++i)
		{
			const std::pair<Vertex, Vertex>& evs = (*ig.edge_incident_vertices_)[edges[i].index_];
			if (i == 0)
			{
				if (edges_dir[i] == 0)
				{
					if (!func(evs.first))
						break;
					if (!func(evs.second))
						break;
				}
				else
				{
					if (!func(evs.second))
						break;
					if (!func(evs.first))
						break;
				}
			}
			else
			{
				if (edges_dir[i] == 0)
				{
					if (!func(evs.second))
						break;
				}
				else
				{
					if (!func(evs.first))
						break;
				}
			}
		}
	}
}

template <typename FUNC>
auto foreach_adjacent_vertex_through_edge(const IncidenceGraph& ig, IncidenceGraph::Vertex v, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	for (IncidenceGraph::Edge e : (*ig.vertex_incident_edges_)[v.index_])
	{
		const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& ev = (*ig.edge_incident_vertices_)[e.index_];
		if (ev.first.index_ != v.index_)
		{
			if (!func(ev.first))
				break;
		}
		else
		{
			if (!func(ev.second))
				break;
		}
	}
}




template <typename CELL, typename FUNC>
auto foreach_incident_edge(const IncidenceGraph& ig, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this IncidenceGraph");
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, IncidenceGraph::Vertex>)
	{
		for (auto& ep : (*ig.vertex_incident_edges_)[c.index_])
		{
			if (!func(ep))
				break;
		}
	}
	else if constexpr (std::is_same_v<CELL, IncidenceGraph::Face>)
	{
		for (auto& ep : (*ig.face_incident_edges_)[c.index_])
		{
			if (!func(ep))
				break;
		}
	}
}

template <typename FUNC>
void foreach_adjacent_edge_through_face(const IncidenceGraph& ig, IncidenceGraph::Edge e, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	bool stop = false;
	CellMarkerStore<IncidenceGraph, IncidenceGraph::Edge> marker(ig);
	marker.mark(e);
	foreach_incident_face(ig, e, [&](IncidenceGraph::Face f0) -> bool {
		foreach_incident_edge(ig, f0, [&](IncidenceGraph::Edge e1) -> bool {
			if (!marker.is_marked(e1))
			{
				marker.mark(e1);
				stop = !func(e1);
			}
			return !stop;
		});
		return !stop;
	});
}

// FACE
template <typename CELL, typename FUNC>
auto foreach_incident_face(const IncidenceGraph& ig, CELL c, const FUNC& func)
{
	static_assert(is_in_tuple<CELL, mesh_traits<IncidenceGraph>::Cells>::value,
				  "CELL not supported in this IncidenceGraph");
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, IncidenceGraph::Vertex>)
	{
		CellMarkerStore<IncidenceGraph, IncidenceGraph::Face> marker(ig);
		for (auto& ep : (*ig.vertex_incident_edges_)[c.index_])
		{
			bool stop = false;
			for (auto& fp : (*ig.edge_incident_faces_)[ep.index_])
			{
				if (!marker.is_marked(fp))
				{
					marker.mark(fp);
					stop = !func(fp);
					if (stop)
						break;
				}
			}
			if (stop)
				break;
		}
	}
	else if constexpr (std::is_same_v<CELL, IncidenceGraph::Edge>)
	{
		for (auto& fp : (*ig.edge_incident_faces_)[c.index_])
		{
			if (!func(fp))
				break;
		}
	}
}

template <typename FUNC>
auto foreach_adjacent_face_through_edge(const IncidenceGraph& ig, IncidenceGraph::Face f, const FUNC& func)
{
	static_assert(is_func_parameter_same<FUNC, IncidenceGraph::Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	bool stop = false;
	CellMarkerStore<IncidenceGraph, IncidenceGraph::Face> marker(ig);
	marker.mark(f);
	for (auto& ie : (*ig.face_incident_edges_)[f.index_])
	{
		for (auto& iface : (*ig.edge_incident_faces_)[ie.index_])
		{
			if (!marker.is_marked(iface))
			{
				marker.mark(iface);
				stop = !func(iface);
			}
			if (stop)
				break;
		}
		if (stop)
			break;
	}
}



} // namespace cgogn

#endif // CGOGN_CORE_INCIDENCE_GRAPH_LOCAL_TRAVERSAL_H_
