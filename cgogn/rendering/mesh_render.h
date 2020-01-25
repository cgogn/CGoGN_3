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
	FACES_CENTERS,
	VOLUMES
};

enum DrawingBufferType : uint8
{
	BUFFER_POINTS = 0,
	BUFFER_LINES,
	BUFFER_TRIANGLES,

	BUFFER_VOLUMES_FACES,
	BUFFER_VOLUMES_EDGES,
	BUFFER_VOLUMES_VERTICES,

	BUFFER_EMB_EDGES,
	BUFFER_EMB_FACES,
	BUFFER_EMB_VOLUMES,

	SIZE_BUFFER,
	BUFFER_POINTS_TB,
	BUFFER_LINES_TB,
	BUFFER_TRIANGLES_TB,

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

	inline bool is_primitive_uptodate(DrawingBufferType prim)
	{
		return indices_buffers_uptodate_[prim];
	}

	inline void set_primitive_dirty(DrawingBufferType prim)
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

//	template <typename MESH>
//	inline void init_lines(const MESH& m, std::vector<uint32>& table_indices)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 0)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Edge = typename mesh_traits<MESH>::Edge;
//			foreach_cell(m, [&](Edge e) -> bool {
//				foreach_incident_vertex(m, e, [&](Vertex v) -> bool {
//					table_indices.push_back(index_of(m, v));
//					return true;
//				});
//				return true;
//			});
//		}
//	}

	template <bool EMB, typename MESH>
	inline void init_lines(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_emb_edge)
	{
		if constexpr (mesh_traits<MESH>::dimension > 0)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			foreach_cell(m, [&](Edge e) -> bool
			{
				foreach_incident_vertex(m, e, [&](Vertex v) -> bool
				{
					table_indices.push_back(index_of(m, v));
					return true;
				});

				if (EMB)
					table_emb_edge.push_back(index_of(m,e));
				else
					table_emb_edge.push_back(table_emb_edge.size());

				return true;
			});
		}
	}

//	template <typename MESH>
//	inline void init_triangles(const MESH& m, std::vector<uint32>& table_indices)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 1)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Face = typename mesh_traits<MESH>::Face;
//			foreach_cell(m, [&](Face f) -> bool {
//				std::vector<Vertex> vertices = incident_vertices(m, f);
//				for (uint32 i = 1; i < vertices.size() - 1; ++i)
//				{
//					table_indices.push_back(index_of(m, vertices[0]));
//					table_indices.push_back(index_of(m, vertices[i]));
//					table_indices.push_back(index_of(m, vertices[i + 1]));
//				}
//				return true;
//			});
//		}
//	}


	template <bool EMB, typename MESH>
	inline void init_triangles(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_emb_face)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			std::vector<Vertex> vertices;
			vertices.reserve(32u);
			foreach_cell(m, [&](Face f) -> bool
			{
				vertices.clear();
				incident_vertices(m, f, vertices);
				for (uint32 i = 1; i < vertices.size() - 1; ++i)
				{
					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[i]));
					table_indices.push_back(index_of(m, vertices[i + 1]));
				}
				if (EMB)
					table_emb_face.push_back(index_of(m,f));
				else
					table_emb_face.push_back(table_emb_face.size());
				return true;
			});
		}
	}

//	template <typename MESH>
//	inline void init_ear_triangles(const MESH& m, std::vector<uint32>& table_indices,
//				   const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 1)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Face = typename mesh_traits<MESH>::Face;
//			std::vector<Vertex> vertices;
//			vertices.reserve(32u);
//			foreach_cell(m, [&](Face f) -> bool
//			{
//				if (codegree(m,f)==3)
//				{
//					vertices.clear();
//					incident_vertices(m, f, vertices);
//					table_indices.push_back(index_of(m, vertices[0]));
//					table_indices.push_back(index_of(m, vertices[1]));
//					table_indices.push_back(index_of(m, vertices[2]));
//				}
//				else
//				{
//					cgogn::geometry::append_ear_triangulation(m, f, position, table_indices);
//				}
//				return true;
//			});
//		}
//	}

	template <bool EMB, typename MESH>
	inline void init_ear_triangles(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_emb_face,
				   const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			std::vector<Vertex> vertices;
			vertices.reserve(32u);
			foreach_cell(m, [&](Face f) -> bool
			{
				if (codegree(m,f)==3)
				{
					vertices.clear();
					incident_vertices(m, f, vertices);

					table_indices.push_back(index_of(m, vertices[0]));
					table_indices.push_back(index_of(m, vertices[1]));
					table_indices.push_back(index_of(m, vertices[2]));
					table_emb_face.push_back(EMB ? index_of(m, f) : table_emb_face.size());
				}
				else
				{
					cgogn::geometry::append_ear_triangulation(m, f, position, table_indices,
					  [&] ()
					  {
						  table_indices.push_back(EMB ? index_of(m, f) : table_emb_face.size());
					  });
				}
				return true;
			});
		}
	}


//	template <bool EMB, typename MESH>
//	inline void init_faces_centers(const MESH& m, std::vector<uint32>& table_indices, std::vector<uint32>& table_emb_edge)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 0)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Face = typename mesh_traits<MESH>::Edge;
//			uint32 last = 0;
//			foreach_cell(m, [&](Face f) -> bool
//			{
//				foreach_incident_vertex(m, f, [&](Vertex v) -> bool
//				{
//					table_indices.push_back(index_of(m, v));
//					return true;
//					if (EMB)
//						table_indices.push_back(index_of(m,f));
//					else
//						table_indices.push_back(last++);
//				});
//				return true;
//			});
//		}
//	}

	template <bool EMB, typename MESH>
	inline void init_volumes(const MESH& m,
							 std::vector<uint32>& table_indices_f,
							 std::vector<uint32>& table_indices_e,
							 std::vector<uint32>& table_indices_v,
							 std::vector<uint32>& table_emb_vol,
			const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{

		if constexpr (mesh_traits<MESH>::dimension > 2)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			using Face = typename mesh_traits<MESH>::Face;
			using Volume = typename mesh_traits<MESH>::Volume;

			std::vector<Vertex> vertices;
			vertices.reserve(256);
			foreach_cell(m, [&](Volume vol) -> bool
			{
				foreach_incident_face(m, vol, [&](Face f)-> bool
				{
					if (codegree(m,f)==3)
					{
						vertices.clear();
						incident_vertices(m, f, vertices);
						table_indices_f.push_back(index_of(m, vertices[0]));
						table_indices_f.push_back(index_of(m, vertices[1]));
						table_indices_f.push_back(index_of(m, vertices[2]));
						table_indices_f.push_back(EMB ? index_of(m, vol) : table_emb_vol.size());
					}
					else
						cgogn::geometry::append_ear_triangulation(m, f, position, table_indices_f,
							[&] ()
						{
							table_indices_f.push_back(EMB ? index_of(m, vol) : table_emb_vol.size());
						});

//					if (EMB)
//						table_indices_f.push_back(index_of(m, vol));
//					else
//						table_indices_f.push_back(table_emb_vol.size());
					return true;
				});

				foreach_incident_edge(m, vol, [&](Edge e)-> bool
				{
					vertices.clear();
					incident_vertices(m, e, vertices);

					table_indices_e.push_back(index_of(m, vertices[0]));
					table_indices_e.push_back(index_of(m, vertices[1]));

					if (EMB)
						table_indices_e.push_back(index_of(m, vol));
					else
						table_indices_e.push_back(table_emb_vol.size());
					return true;
				});

				foreach_incident_vertex(m, vol, [&](Vertex v)-> bool
				{
					table_indices_v.push_back(index_of(m, v));

					if (EMB)
						table_indices_v.push_back(index_of(m, vol));
					else
						table_indices_v.push_back(table_emb_vol.size());
					return true;
				});

				if (EMB)
					table_emb_vol.push_back(index_of(m, vol));
				else
					table_emb_vol.push_back(table_emb_vol.size());

				return true;
			});
		}
	}

//	template <typename MESH, bool EMB>
//	inline void init_volumes_faces(const MESH& m, std::vector<uint32>& table_indices,
//			const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 2)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Face = typename mesh_traits<MESH>::Face;
//			using Volume = typename mesh_traits<MESH>::Volume;
//			foreach_cell(m, [&](Volume vol) -> bool
//			{
//				foreach_incident_face(m, vol, [&](Face f)-> bool
//				{
//					if (codegree(m,f)==3)
//					{
//						std::vector<Vertex> vertices = incident_vertices(m, f);
//						table_indices.push_back(index_of(m, vertices[0]));
//						table_indices.push_back(index_of(m, vertices[1]));
//						table_indices.push_back(index_of(m, vertices[2]));
//					}
//					else
//						cgogn::geometry::append_ear_triangulation(m, f, position, table_indices);

//					if (EMB)
//						table_indices.push_back(index_of(m, vol));
//					else
//						table_indices.push_back(table_indices.size());

//					return true;
//				});
//				return true;
//			});
//		}
//	}

//	template <typename MESH, bool EMB>
//	inline void init_volumes_edges(const MESH& m, std::vector<uint32>& table_indices)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 2)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Edge = typename mesh_traits<MESH>::Edge;
//			using Volume = typename mesh_traits<MESH>::Volume;
//			foreach_cell(m, [&](Volume vol) -> bool
//			{
//				foreach_incident_edge(m, vol, [&](Edge e)-> bool
//				{
//					std::vector<Vertex> vertices = incident_vertices(m, e);
//					table_indices.push_back(index_of(m, vertices[0]));
//					table_indices.push_back(index_of(m, vertices[1]));
//					if (EMB)
//						table_indices.push_back(index_of(m, vol));
//					else
//						table_indices.push_back(table_indices.size());
//					return true;
//				});
//				return true;
//			});
//		}
//	}

//	template <typename MESH>
//	inline void init_volumes_vertices(const MESH& m, std::vector<uint32>& table_indices)
//	{
//		if constexpr (mesh_traits<MESH>::dimension > 2)
//		{
//			using Vertex = typename mesh_traits<MESH>::Vertex;
//			using Volume = typename mesh_traits<MESH>::Volume;
//			foreach_cell(m, [&](Volume vol) -> bool
//			{
//				foreach_incident_vertex(m, vol, [&](Vertex v)-> bool
//				{
//					table_indices.push_back(index_of(m, v));
//					table_indices.push_back(index_of(m, vol));
//					return true;
//				});
//				return true;
//			});
//		}
//	}


public:
	inline uint32 nb_indices(DrawingType prim) const
	{
		return nb_indices_[prim];
	}

	template <typename MESH>
	inline void init_primitives(const MESH& m, DrawingType prim,
								const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position = nullptr)
	{

		auto start_timer = std::chrono::high_resolution_clock::now();

		std::vector<uint32> table_indices;
		table_indices.reserve(1024u);
		std::vector<uint32> table_indices_emb;
		table_indices_emb.reserve(1024u);

		std::vector<uint32> table_indices_e;
		table_indices_emb.reserve(1024u);
		std::vector<uint32> table_indices_v;
		table_indices_emb.reserve(1024u);

		using Edge = typename mesh_traits<MESH>::Edge;
		using Face = typename mesh_traits<MESH>::Face;
		using Volume = typename mesh_traits<MESH>::Volume;

		switch (prim)
		{
		case POINTS:
			init_points(m, table_indices);
			break;
		case LINES:
			if (is_indexed<Edge>(m))
				init_lines<true>(m, table_indices,table_indices_emb);
			else
				init_lines<false>(m, table_indices,table_indices_emb);

			break;
		case TRIANGLES:
			if (position == nullptr)
			{
				if (is_indexed<Face>(m))
					init_triangles<true>(m, table_indices, table_indices_emb);
				else
					init_triangles<false>(m, table_indices, table_indices_emb);
			}
			else
			{
				if (is_indexed<Face>(m))
					init_ear_triangles<true>(m, table_indices, table_indices_emb, position);
				else
					init_ear_triangles<false>(m, table_indices, table_indices_emb, position);
			}
			break;
//		case FACES_CENTERS:
//				if (is_indexed<Face>(m))
//					init_faces_center<true>(m, table_indices, position);
//				else
//					init_faces_center<false>(m, table_indices, position);
//				break;
		case VOLUMES:
				if (is_indexed<Volume>(m))
					init_volumes<true>(m, table_indices, table_indices_e,table_indices_v, table_indices_emb, position);
				else
					init_volumes<false>(m, table_indices, table_indices_e,table_indices_v, table_indices_emb, position);
				break;
		default:
			break;
		}



		auto func_update_ebo = [&] (DrawingBufferType pr, const std::vector<uint32>& table ) -> void
		{
			indices_buffers_uptodate_[pr] = true;
			nb_indices_[pr] = uint32(table.size());

			if (!table.empty())
			{
				if (!indices_buffers_[pr]->is_created())
					indices_buffers_[pr]->create();

				indices_buffers_[pr]->bind();
				indices_buffers_[pr]->allocate(table.data(), table.size());
				indices_buffers_[pr]->release();
			}
		};

		switch (prim)
		{
		case POINTS:
				func_update_ebo(BUFFER_POINTS,table_indices);
				break;
		case LINES:
				func_update_ebo(BUFFER_LINES,table_indices);
				func_update_ebo(BUFFER_EMB_EDGES,table_indices_emb);
				break;
		case TRIANGLES:
				func_update_ebo(BUFFER_TRIANGLES,table_indices);
				func_update_ebo(BUFFER_EMB_FACES,table_indices_emb);
				break;
			break;
		case VOLUMES:
				func_update_ebo(BUFFER_VOLUMES_FACES,table_indices);
				func_update_ebo(BUFFER_VOLUMES_EDGES,table_indices_e);
				func_update_ebo(BUFFER_VOLUMES_VERTICES,table_indices_v);
				func_update_ebo(BUFFER_EMB_VOLUMES, table_indices_emb);
				break;
		default:
			break;
		}

		auto end_timer = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
		std::cout << "init primitive "<<prim<< " in "<< elapsed_seconds.count() << std::endl;

	}

	void draw(DrawingBufferType prim,GLint binding_point=10);

	template <typename MESH>
	inline void init_primitives(const MESH& m, DrawingBufferType prim,
								const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position = nullptr)
	{
		if (prim<= BUFFER_TRIANGLES)
			init_primitives(m, static_cast<DrawingType>(prim), position);
		else if (prim<= BUFFER_VOLUMES_VERTICES)
			init_primitives(m, VOLUMES, position);
	}


	inline void bind_ebo_tb(DrawingBufferType prim, GLint binding_point)
	{
		if (prim >= SIZE_BUFFER)
			indices_buffers_[prim - SIZE_BUFFER - 1]->bind_tb(binding_point);
		else
			indices_buffers_[prim]->bind_tb(binding_point);
	}



};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_MESH_RENDER_H_
