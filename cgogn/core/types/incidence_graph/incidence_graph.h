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

// #include <cgogn/core/cgogn_core_export.h>

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

	MarkAttribute* vertex_boundary_marker_;
	MarkAttribute* edge_boundary_marker_;
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
		vertex_boundary_marker_ = attribute_containers_[Vertex::CELL_INDEX].get_mark_attribute();
		edge_boundary_marker_ = attribute_containers_[Edge::CELL_INDEX].get_mark_attribute();
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

} // namespace cgogn

#endif // CGOGN_CORE_INCIDENCE_GRAPH_H_
