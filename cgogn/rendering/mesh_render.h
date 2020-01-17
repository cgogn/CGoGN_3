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

#ifndef CGOGN_RENDERING_MESH_RENDER_H_
#define CGOGN_RENDERING_MESH_RENDER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>

#include <cgogn/core/utils/numerics.h>

 #include <cgogn/geometry/algos/ear_triangulation.h>

#include <memory>

#include <cgogn/rendering/drawer.h>
#include <cgogn/rendering/ebo.h>
#include <cgogn/rendering/vbo.h>

namespace cgogn
{

namespace rendering
{

enum DrawingType : uint8
{
	POINTS = 0,
	LINES,
	TRIANGLES,
	BOUNDARY,
	VOLUMES_VERTICES,
	VOLUMES_EDGES,
	VOLUMES_FACES,
	SIZE_BUFFER,
	LINES_TB,
	TRIANGLES_TB,
};

class CGOGN_RENDERING_EXPORT MeshRender
{
protected:
	std::array<std::unique_ptr<EBO>, SIZE_BUFFER> indices_buffers_;
	std::array<bool, SIZE_BUFFER> indices_buffers_uptodate_;
	std::array<uint32, SIZE_BUFFER> nb_indices_;

public:
	MeshRender();
	~MeshRender();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(MeshRender);

	inline bool is_primitive_uptodate(DrawingType prim)
	{
		return indices_buffers_uptodate_[prim];
	}
	inline void set_primitive_dirty(DrawingType prim)
	{
		indices_buffers_uptodate_[prim] = false;
	}

protected:
	template <typename MESH>
	inline void init_points(const MESH& m, std::vector<uint32>& table_indices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		foreach_cell(m, [&](Vertex v) -> bool {
			table_indices.push_back(index_of(m, v));
			return true;
		});
	}

	template <typename MESH>
	inline void init_lines(const MESH& m, std::vector<uint32>& table_indices)
	{
		if constexpr (mesh_traits<MESH>::dimension > 0)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			foreach_cell(m, [&](Edge e) -> bool {
				foreach_incident_vertex(m, e, [&](Vertex v) -> bool {
					table_indices.push_back(index_of(m, v));
					return true;
				});
				return true;
			});
		}
	}

	template <typename MESH>
	inline void init_lines(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_indices_edge)
	{
		if constexpr (mesh_traits<MESH>::dimension > 0)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			foreach_cell(m, [&](Edge e) -> bool {
				foreach_incident_vertex(m, e, [&](Vertex v) -> bool {
					table_indices.push_back(index_of(m, v));
					return true;
				});
				table_indices.push_back(index_of(m, e));
				return true;
			});
		}
	}

	template <typename MESH>
	inline void init_triangles(const MESH& m, std::vector<uint32>& table_indices)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			foreach_cell(m, [&](Face f) -> bool {
				std::vector<Vertex> vertices = incident_vertices(m, f);
				for (uint32 i = 1; i < vertices.size() - 1; ++i)
				{
					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[i]));
					table_indices.push_back(index_of(m, vertices[i + 1]));
				}
				return true;
			});
		}
	}


	template <typename MESH>
	inline void init_triangles(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_indices_face)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			foreach_cell(m, [&](Face f) -> bool {
				std::vector<Vertex> vertices = incident_vertices(m, f);
				for (uint32 i = 1; i < vertices.size() - 1; ++i)
				{
					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[i]));
					table_indices.push_back(index_of(m, vertices[i + 1]));
				}
				table_indices_face.push_back(index_of(m,f));
				return true;
			});
		}
	}

	template <typename MESH>
	inline void init_ear_triangles(const MESH& m, std::vector<uint32>& table_indices,
				   const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			foreach_cell(m, [&](Face f) -> bool {
				if (codegree(m,f)==3)
				{
					std::vector<Vertex> vertices = incident_vertices(m, f);
					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[1]));
					table_indices.push_back(index_of(m, vertices[2]));
				}
				else
				{
					cgogn::geometry::append_ear_triangulation(m, f, position, table_indices);
				}
				return true;
			});
		}
	}

	template <typename MESH>
	inline void init_ear_triangles(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_indices_face,
				   const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			foreach_cell(m, [&](Face f) -> bool {
				if (codegree(m,f)==3)
				{
					std::vector<Vertex> vertices = incident_vertices(m,f);
					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[1]));
					table_indices.push_back(index_of(m, vertices[2]));
				}
				else
				{
					cgogn::geometry::append_ear_triangulation(m, f, position, table_indices);
				}
				table_indices_face.push_back(index_of(m,f));
				return true;
			});
		}
	}
	template <typename MESH>
	inline void init_volumes_faces(const MESH& m, std::vector<uint32>& table_indices,
			const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{
		if constexpr (mesh_traits<MESH>::dimension > 2)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			using Volume = typename mesh_traits<MESH>::Volume;
			foreach_cell(m, [&](Volume vol) -> bool {
				foreach_incident_face(m, vol, [&](Face f)-> bool {
					if (codegree(m,f)==3)
					{
						std::vector<Vertex> vertices = incident_vertices(m, f);
						table_indices.push_back(index_of(m, vertices[0]));
						table_indices.push_back(index_of(m, vertices[1]));
						table_indices.push_back(index_of(m, vertices[2]));
					}
					else
					{
						cgogn::geometry::append_ear_triangulation(m, f, position, table_indices);
					}
					table_indices.push_back(index_of(m, vol));
					return true;
				});
				return true;
			});
		}
	}

	template <typename MESH>
	inline void init_volumes_edges(const MESH& m, std::vector<uint32>& table_indices)
	{
		if constexpr (mesh_traits<MESH>::dimension > 2)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			using Volume = typename mesh_traits<MESH>::Volume;
			foreach_cell(m, [&](Volume vol) -> bool {
				foreach_incident_edge(m, vol, [&](Edge e)-> bool {
					std::vector<Vertex> vertices = incident_vertices(m, e);
					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[1]));
					table_indices.push_back(index_of(m, vol));
					return true;
				});
				return true;
			});
		}
	}

	template <typename MESH>
	inline void init_volumes_vertices(const MESH& m, std::vector<uint32>& table_indices)
	{
		if constexpr (mesh_traits<MESH>::dimension > 2)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Volume = typename mesh_traits<MESH>::Volume;
			foreach_cell(m, [&](Volume vol) -> bool {
				foreach_incident_vertex(m, vol, [&](Vertex v)-> bool {
					table_indices.push_back(index_of(m, v));
					table_indices.push_back(index_of(m, vol));
					return true;
				});
				return true;
			});
		}
	}


public:
	inline uint32 nb_indices(DrawingType prim) const
	{
		return nb_indices_[prim];
	}

	template <typename MESH>
	inline void init_primitives(const MESH& m, DrawingType prim,
								const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position = nullptr)
	{
		std::vector<uint32> table_indices;
		table_indices.reserve(1024u);
		std::vector<uint32> table_indices_cell;
		table_indices_cell.reserve(1024u);

		 using Edge = typename mesh_traits<MESH>::Edge;
		 using Face = typename mesh_traits<MESH>::Face;

		switch (prim)
		{
		case POINTS:
			init_points(m, table_indices);
			break;
		case LINES:
			if (is_indexed<Face>(m))
				init_lines(m, table_indices,table_indices_cell);
			else
				init_lines(m, table_indices);
			break;
		case TRIANGLES:
			if (is_indexed<Face>(m))
				init_ear_triangles(m, table_indices, position);
			else
				init_ear_triangles(m, table_indices, table_indices_cell, position);
			break;
		case BOUNDARY:
			break;
		case VOLUMES_EDGES:
			init_volumes_edges(m, table_indices);
			break;
		case VOLUMES_FACES:
			init_volumes_faces(m, table_indices, position);
			break;
		case VOLUMES_VERTICES:
			init_volumes_vertices(m, table_indices);
			break;
		default:
			break;
		}

		indices_buffers_uptodate_[prim] = true;
		nb_indices_[prim] = uint32(table_indices.size());

		if (table_indices.empty())
			return;

		if (!indices_buffers_[prim]->is_created())
			indices_buffers_[prim]->create();

		indices_buffers_[prim]->bind();
		indices_buffers_[prim]->allocate(table_indices.data(), nb_indices_[prim]);
		indices_buffers_[prim]->release();
	}

	void draw(DrawingType prim);
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_MESH_RENDER_H_
