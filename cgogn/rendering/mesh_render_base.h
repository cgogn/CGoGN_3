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

#ifndef CGOGN_RENDERING_MESH_RENDER_BASE_H_
#define CGOGN_RENDERING_MESH_RENDER_BASE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>


#include <cgogn/rendering/ebo.h>
#include <cgogn/rendering/vbo.h>

#include <memory>

namespace cgogn
{

namespace rendering
{

using TablesIndices = std::vector<std::vector<uint32>>;

enum DrawingType : uint32
{
	POINTS = 0,
	LINES,
	TRIANGLES,

	VOLUMES_FACES,
	VOLUMES_EDGES,
	VOLUMES_VERTICES,

	INDEX_EDGES,
	INDEX_FACES,
	INDEX_VOLUMES,

	POINTS_TB,
	LINES_TB,
	TRIANGLES_TB,

	VOLUMES_FACES_TB,
	VOLUMES_EDGES_TB,
	VOLUMES_VERTICES_TB,

	INDEX_EDGES_TB,
	INDEX_FACES_TB,
	INDEX_VOLUMES_TB
};

static const uint32 SIZE_BUFFER = uint32(POINTS_TB);

inline DrawingType& operator++(DrawingType& d)
{
	++*reinterpret_cast<int32*>(&d);
	return d;
}

inline int32* operator&(DrawingType& d)
{
	return reinterpret_cast<int*>(&d);
}

static std::vector<std::string> primitives_names = {
	"POINTS",			"LINES",		  "TRIANGLES",		  "VOLUMES_FACES",	  "VOLUMES_EDGES",
	"VOLUMES_VERTICES", "INDEX_EDGES",	  "INDEX_FACES",	  "INDEX_VOLUMES",	  "POINTS_TB",
	"LINES_TB",			"TRIANGLES_TB",	  "VOLUMES_FACES_TB", "VOLUMES_EDGES_TB", "VOLUMES_VERTICES_TB",
	"INDEX_EDGES_TB",	"INDEX_FACES_TB", "INDEX_VOLUMES_TB"};

class CGOGN_RENDERING_EXPORT MeshRender
{
protected:
	std::array<std::unique_ptr<EBO>, SIZE_BUFFER> indices_buffers_;
	std::array<bool, SIZE_BUFFER> indices_buffers_uptodate_;

public:
	MeshRender();
	~MeshRender();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(MeshRender);

	inline bool is_primitive_uptodate(DrawingType prim)
	{
		return indices_buffers_uptodate_[prim % SIZE_BUFFER];
	}

	inline void set_primitive_dirty(DrawingType prim)
	{
		indices_buffers_uptodate_[prim % SIZE_BUFFER] = false;
	}

	inline void set_all_primitives_dirty()
	{
		for (DrawingType p = POINTS; p < SIZE_BUFFER; ++p)
			indices_buffers_uptodate_[p] = false;
	}

	inline EBO* get_EBO(DrawingType prim)
	{
		return indices_buffers_[prim % SIZE_BUFFER].get();
	}

protected:
	template <typename MESH>
	void init_points(const MESH& m, TablesIndices& table_indices);

	template <bool EMB, typename MESH>
	void init_lines(const MESH& m, TablesIndices& table_indices, TablesIndices& table_emb_edge);

	template <bool EMB, typename MESH>
	void init_triangles(const MESH& m, TablesIndices& table_indices, TablesIndices& table_emb_face);
    
	template <bool EMB, typename MESH>
	void init_ear_triangles(const MESH& m, TablesIndices& table_indices, TablesIndices& table_emb_face,
								   const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position);

	template <bool EMB, typename MESH>
	void init_volumes(const MESH& m, TablesIndices& table_indices_f, TablesIndices& table_indices_e,
							 TablesIndices& table_indices_v, TablesIndices& table_emb_vol,
							 const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position);

public:
	template <typename MESH>
	void init_primitives(
		const MESH& m, DrawingType prim,
		const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position = nullptr);

	void draw(DrawingType prim);
};

} // namespace rendering

} // namespace cgogn

#endif 
