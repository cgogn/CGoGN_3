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

#ifndef CGOGN_CORE_INCIDENCE_GRAPH_H_
#define CGOGN_CORE_INCIDENCE_GRAPH_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/container/attribute_container.h>
#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/types/container/vector.h>
#include <cgogn/core/types/mesh_traits.h>

#include <any>
#include <array>

namespace cgogn
{

struct CGOGN_CORE_EXPORT IncidenceGraph
{
	// using AttributeContainer = AttributeContainerT<Vector>;
	using AttributeContainer = AttributeContainerT<ChunkArray>;

	template <typename T>
	using Attribute = AttributeContainer::Attribute<T>;
	using AttributeGen = AttributeContainer::AttributeGen;
	using MarkAttribute = AttributeContainer::MarkAttribute;

	/*************************************************************************/
	// Graph attributes container
	/*************************************************************************/

	// std::unordered_map<std::string, std::any> attributes_;

	struct Vertex
	{
		static const uint32 CELL_INDEX = 0;
		uint32 index_;
		inline Vertex() : index_(INVALID_INDEX)
		{
		}
		inline Vertex(uint32 id) : index_(id)
		{
		}
		bool operator<(Vertex v) const
		{
			return index_ < v.index_;
		}
		bool operator==(Vertex v) const
		{
			return index_ == v.index_;
		}
		bool operator!=(Vertex v) const
		{
			return index_ != v.index_;
		}
		inline bool is_valid() const
		{
			return index_ != INVALID_INDEX;
		}
	};

	struct Edge
	{
		static const uint32 CELL_INDEX = 1;
		uint32 index_;
		inline Edge() : index_(INVALID_INDEX)
		{
		}
		inline Edge(uint32 id) : index_(id)
		{
		}
		bool operator<(Edge e) const
		{
			return index_ < e.index_;
		}
		bool operator==(Edge e) const
		{
			return index_ == e.index_;
		}
		bool operator!=(Edge e) const
		{
			return index_ != e.index_;
		}
		inline bool is_valid() const
		{
			return index_ != INVALID_INDEX;
		}
	};

	struct Face
	{
		static const uint32 CELL_INDEX = 2;
		uint32 index_;
		inline Face() : index_(INVALID_INDEX)
		{
		}
		inline Face(uint32 id) : index_(id)
		{
		}
		bool operator<(Face f) const
		{
			return index_ < f.index_;
		}
		bool operator==(Face f) const
		{
			return index_ == f.index_;
		}
		bool operator!=(Face f) const
		{
			return index_ != f.index_;
		}
		inline bool is_valid() const
		{
			return index_ != INVALID_INDEX;
		}
	};

	mutable std::array<AttributeContainer, 3> attribute_containers_;

	std::shared_ptr<Attribute<std::vector<Edge>>> vertex_incident_edges_;
	std::shared_ptr<Attribute<std::pair<Vertex, Vertex>>> edge_incident_vertices_;
	std::shared_ptr<Attribute<std::vector<Face>>> edge_incident_faces_;
	std::shared_ptr<Attribute<std::vector<Edge>>> face_incident_edges_;
	std::shared_ptr<Attribute<std::vector<uint8>>> face_incident_edges_dir_;

	IncidenceGraph()
	{
		vertex_incident_edges_ =
			attribute_containers_[Vertex::CELL_INDEX].add_attribute<std::vector<Edge>>("incident_edges");
		edge_incident_vertices_ =
			attribute_containers_[Edge::CELL_INDEX].add_attribute<std::pair<Vertex, Vertex>>("incident_vertices");
		edge_incident_faces_ =
			attribute_containers_[Edge::CELL_INDEX].add_attribute<std::vector<Face>>("incident_faces");
		face_incident_edges_ =
			attribute_containers_[Face::CELL_INDEX].add_attribute<std::vector<Edge>>("incident_edges");
		face_incident_edges_dir_ =
			attribute_containers_[Face::CELL_INDEX].add_attribute<std::vector<uint8>>("incident_edges_dir");
	};
	// ~IncidenceGraph();
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
	using Attribute = IncidenceGraph::Attribute<T>;
	using AttributeGen = IncidenceGraph::AttributeGen;
	using MarkAttribute = IncidenceGraph::MarkAttribute;
};



bool sort_face_edges(IncidenceGraph& ig, IncidenceGraph::Face f);

bool same_edge(IncidenceGraph& ig, IncidenceGraph::Edge e1, IncidenceGraph::Edge e2);

void remove_edge_in_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v, IncidenceGraph::Edge edge_to_remove);


void remove_face_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Face face_to_remove);

void remove_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge edge_to_remove);

void replace_vertex_in_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, IncidenceGraph::Vertex old_vertex,
								   IncidenceGraph::Vertex new_vertex);

void replace_edge_in_face(IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge old_edge,
								 IncidenceGraph::Edge new_edge);


IncidenceGraph::Vertex common_vertex(IncidenceGraph& ig, IncidenceGraph::Edge e0, IncidenceGraph::Edge e1);

std::vector<IncidenceGraph::Vertex> sorted_face_vertices(IncidenceGraph& ig, IncidenceGraph::Face f);


IncidenceGraph::Vertex CGOGN_CORE_EXPORT add_vertex(IncidenceGraph& ig);

void CGOGN_CORE_EXPORT remove_vertex(IncidenceGraph& ig, IncidenceGraph::Vertex v);

IncidenceGraph::Edge CGOGN_CORE_EXPORT connect_vertices(IncidenceGraph& g, IncidenceGraph::Vertex v1, IncidenceGraph::Vertex v2);


IncidenceGraph::Edge CGOGN_CORE_EXPORT add_edge(IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1);

IncidenceGraph::Vertex CGOGN_CORE_EXPORT cut_edge(IncidenceGraph& ig, IncidenceGraph::Edge e, bool set_indices = true);

std::pair<IncidenceGraph::Vertex, std::vector<IncidenceGraph::Edge>> collapse_edge(IncidenceGraph& ig,IncidenceGraph::Edge e, bool set_indices = true);

void CGOGN_CORE_EXPORT remove_edge(IncidenceGraph& ig, IncidenceGraph::Edge e);

IncidenceGraph::Face CGOGN_CORE_EXPORT add_face(IncidenceGraph& ig, std::vector<IncidenceGraph::Edge>& edges);

void CGOGN_CORE_EXPORT remove_face(IncidenceGraph& ig, IncidenceGraph::Face f);

IncidenceGraph::Edge CGOGN_CORE_EXPORT cut_face(IncidenceGraph& m, IncidenceGraph::Vertex v1,
												IncidenceGraph::Vertex v2);

inline void copy(IncidenceGraph& /*dst*/, const IncidenceGraph& /*src*/)
{
	// TODO
}


template <typename CELL>
CELL add_cell(IncidenceGraph& ig)
{
	return CELL(ig.attribute_containers_[CELL::CELL_INDEX].new_index());
}

template <typename CELL>
void remove_cell(IncidenceGraph& ig, CELL c)
{
	ig.attribute_containers_[CELL::CELL_INDEX].release_index(c.index_);
}




template <typename CELL>
bool is_indexed(const IncidenceGraph& /*m*/)
{
	return true;
}



template <typename CELL>
uint32 index_of(const IncidenceGraph& /*m*/, CELL c)
{
	return c.index_;
}


template <typename CELL>
uint32 new_index(const IncidenceGraph& ig)
{
	uint32 id = ig.attribute_containers_[CELL::CELL_INDEX].new_index();
	// (*ig.cells_indices_[CELL::CELL_INDEX])[id] = id;
	return id;
}

} // namespace cgogn

#include <cgogn/core/types/incidence_graph/incidence_graph_attribute.hpp>
#include <cgogn/core/types/incidence_graph/incidence_graph_local_traversals.hpp>
#include <cgogn/core/types/incidence_graph/incidence_graph_global_traversals.hpp>
#endif // CGOGN_CORE_INCIDENCE_GRAPH_H_
