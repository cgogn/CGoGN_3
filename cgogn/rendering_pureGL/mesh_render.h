/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

// #include <cgogn/geometry/algos/ear_triangulation.h>

#include <memory>

#include <cgogn/rendering_pureGL/drawer.h>
#include <cgogn/rendering_pureGL/vbo.h>
#include <cgogn/rendering_pureGL/ebo.h>

namespace cgogn
{

namespace rendering_pgl
{

enum DrawingType : uint8
{
	POINTS = 0,
	LINES,
	TRIANGLES,
	BOUNDARY,
	SIZE_BUFFER
};

class CGOGN_RENDERING_PUREGL_EXPORT MeshRender
{
protected:

	std::array<std::unique_ptr<EBO>, SIZE_BUFFER> indices_buffers_;
	std::array<bool, SIZE_BUFFER> indices_buffers_uptodate_;
	std::array<uint32, SIZE_BUFFER> nb_indices_;

public:

	MeshRender();
	~MeshRender();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(MeshRender);

	inline bool is_primitive_uptodate(DrawingType prim) { return indices_buffers_uptodate_[prim]; }
	inline void set_primitive_dirty(DrawingType prim) { indices_buffers_uptodate_[prim] = false; }

protected:

	template <typename MESH>
	inline void init_points(const MESH& m, std::vector<uint32>& table_indices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		foreach_cell(m, [&] (Vertex v) -> bool
		{
			table_indices.push_back(index_of(m, v));
			return true;
		});
	}

	template <typename MESH>
	inline void init_lines(const MESH& m, std::vector<uint32>& table_indices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Edge = typename mesh_traits<MESH>::Edge;
		foreach_cell(m, [&] (Edge e) -> bool
		{
			foreach_incident_vertex(m, e, [&] (Vertex v) -> bool { table_indices.push_back(index_of(m, v)); return true; });
			return true;
		});
	}

	template <typename MESH>
	inline void init_triangles(const MESH& m, std::vector<uint32>& table_indices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Face = typename mesh_traits<MESH>::Face;
		foreach_cell(m, [&] (Face f) -> bool
		{
			std::vector<Vertex> vertices = incident_vertices(m, f);
			for (uint32 i = 1; i < vertices.size() - 1; ++i)
			{
				table_indices.push_back(index_of(m, vertices[0]));
				table_indices.push_back(index_of(m, vertices[i]));
				table_indices.push_back(index_of(m, vertices[i+1]));
			}
			return true;
		});
	}

public:

	inline uint32 nb_indices(DrawingType prim) const
	{
		return nb_indices_[prim];
	}

	template <typename MESH>
	inline void init_primitives(const MESH& m, DrawingType prim)
	{
		std::vector<uint32> table_indices;

		switch (prim)
		{
			case POINTS:
				init_points(m, table_indices);
				break;
			case LINES:
				init_lines(m, table_indices);
				break;
			case TRIANGLES:
				init_triangles(m, table_indices);
				break;
			case BOUNDARY:
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
