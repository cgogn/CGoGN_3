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

#ifndef CGOGN_CORE_TYPES_INCIDENCE_GRAPH_H_
#define CGOGN_CORE_TYPES_INCIDENCE_GRAPH_H_

#include <cgogn/core/types/incidence_graph/incidence_graph_base.h>

namespace cgogn
{

struct IncidenceGraph : public IncidenceGraphBase
{
	static const uint8 dimension = 2;

	using Vertex = IncidenceGraphBase::Vertex;
	using Edge = IncidenceGraphBase::Edge;
	using Face = IncidenceGraphBase::Face;

	using Cells = std::tuple<Vertex, Edge, Face>;

	IncidenceGraph() : IncidenceGraphBase()
	{
	}
};

template <>
struct mesh_traits<IncidenceGraph>
{
	static constexpr const char* name = "IncidenceGraph";
	static constexpr const uint8 dimension = 2;

	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	using Cells = std::tuple<Vertex, Edge, Face>;
	static constexpr const char* cell_names[] = {"Vertex", "Edge", "Face"};

	template <typename T>
	using Attribute = IncidenceGraphBase::Attribute<T>;
	using AttributeGen = IncidenceGraphBase::AttributeGen;
	using MarkAttribute = IncidenceGraphBase::MarkAttribute;
};

/*************************************************************************/
// Local vertex traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_vertex(const MESH& ig, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraph&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
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

template <typename MESH, typename FUNC>
auto foreach_adjacent_vertex_through_edge(const MESH& ig, typename mesh_traits<MESH>::Vertex v, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraph&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	static_assert(is_func_parameter_same<FUNC, Vertex>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	for (Edge e : (*ig.vertex_incident_edges_)[v.index_])
	{
		const std::pair<Vertex, Vertex>& ev = (*ig.edge_incident_vertices_)[e.index_];
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

/*************************************************************************/
// Local edge traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_edge(const MESH& ig, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraph&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(has_cell_type_v<MESH, CELL>, "CELL not supported in this MESH");
	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, Vertex>)
	{
		for (auto& ep : (*ig.vertex_incident_edges_)[c.index_])
		{
			if (!func(ep))
				break;
		}
	}
	else if constexpr (std::is_same_v<CELL, Face>)
	{
		for (auto& ep : (*ig.face_incident_edges_)[c.index_])
		{
			if (!func(ep))
				break;
		}
	}
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_edge_through_face(const MESH& ig, typename mesh_traits<MESH>::Edge e, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraph&>>
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_func_parameter_same<FUNC, Edge>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	bool stop = false;
	CellMarkerStore<MESH, Edge> marker(ig);
	marker.mark(e);
	foreach_incident_face(ig, e, [&](Face f0) -> bool {
		foreach_incident_edge(ig, f0, [&](Edge e1) -> bool {
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

/*************************************************************************/
// Local face traversals
/*************************************************************************/

template <typename MESH, typename CELL, typename FUNC>
auto foreach_incident_face(const MESH& ig, CELL c, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraph&>>
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_in_tuple<CELL, typename mesh_traits<MESH>::Cells>::value,
				  "CELL not supported in this IncidenceGraph");
	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	if constexpr (std::is_same_v<CELL, Vertex>)
	{
		CellMarkerStore<MESH, Face> marker(ig);
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
	else if constexpr (std::is_same_v<CELL, Edge>)
	{
		for (auto& fp : (*ig.edge_incident_faces_)[c.index_])
		{
			if (!func(fp))
				break;
		}
	}
}

template <typename MESH, typename FUNC>
auto foreach_adjacent_face_through_edge(const MESH& ig, typename mesh_traits<MESH>::Face f, const FUNC& func)
	-> std::enable_if_t<std::is_convertible_v<MESH&, IncidenceGraph&>>
{
	using Face = typename mesh_traits<MESH>::Face;

	static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function cell parameter type");
	static_assert(is_func_return_same<FUNC, bool>::value, "Given function should return a bool");

	bool stop = false;
	CellMarkerStore<MESH, Face> marker(ig);
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

/*************************************************************************/
// Operators
/*************************************************************************/

IncidenceGraph::Vertex add_vertex(IncidenceGraph& ig);
IncidenceGraph::Edge connect_vertices(IncidenceGraph& g, IncidenceGraph::Vertex v1, IncidenceGraph::Vertex v2);
void remove_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v);

IncidenceGraph::Edge add_edge(IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1);
IncidenceGraph::Vertex cut_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, bool set_indices = true);
std::pair<IncidenceGraph::Vertex, std::vector<IncidenceGraph::Edge>> collapse_edge(IncidenceGraph& ig,
																				   IncidenceGraph::Edge e,
																				   bool set_indices = true);
void remove_edge(IncidenceGraph& ig, IncidenceGraph::Edge e);

IncidenceGraph::Face add_face(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge>& edges);
IncidenceGraph::Edge cut_face(IncidenceGraph& m, IncidenceGraph::Vertex v1, IncidenceGraph::Vertex v2);
void remove_face(IncidenceGraph& ig, IncidenceGraph::Face f);

/*************************************************************************/
// Helper functions
/*************************************************************************/

bool same_edge(IncidenceGraph& ig, IncidenceGraph::Edge e1, IncidenceGraph::Edge e2);
IncidenceGraph::Vertex common_vertex(IncidenceGraph& ig, IncidenceGraph::Edge e0, IncidenceGraph::Edge e1);

void remove_edge_in_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v, IncidenceGraph::Edge edge_to_remove);

void remove_face_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Face face_to_remove);
void replace_vertex_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Vertex old_vertex,
							IncidenceGraph::Vertex new_vertex);

void remove_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge edge_to_remove);
void replace_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge old_edge,
						  IncidenceGraph::Edge new_edge);
std::vector<IncidenceGraph::Vertex> sorted_face_vertices(IncidenceGraph& ig, IncidenceGraph::Face f);
bool sort_face_edges(IncidenceGraph& ig, IncidenceGraph::Face f);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_INCIDENCE_GRAPH_H_
