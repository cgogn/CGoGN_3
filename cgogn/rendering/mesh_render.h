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

	SIZE_BUFFER,
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

static std::vector<std::string> primitives_names = {
"POINTS",
"LINES",
"TRIANGLES",
"VOLUMES_FACES",
"VOLUMES_EDGES",
"VOLUMES_VERTICES",
"INDEX_EDGES",
"INDEX_FACES",
"INDEX_VOLUMES",
"SIZE_BUFFER",
"POINTS_TB",
"LINES_TB",
"TRIANGLES_TB",
"VOLUMES_FACES_TB",
"VOLUMES_EDGES_TB",
"VOLUMES_VERTICES_TB",
"INDEX_EDGES_TB",
"INDEX_FACES_TB",
"INDEX_VOLUMES_TB"
};


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
		return indices_buffers_uptodate_[prim%SIZE_BUFFER];
	}

	inline void set_primitive_dirty(DrawingType prim)
	{
		indices_buffers_uptodate_[prim%SIZE_BUFFER] = false;
	}

protected:
	template <typename MESH>
	inline void init_points(const MESH& m, TablesIndices& table_indices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		parallel_foreach_cell(m, [&](Vertex v) -> bool
		{
			table_indices[current_worker_index()].push_back(index_of(m, v));
			return true;
		});
	}


	template <bool EMB, typename MESH>
	inline void init_lines(const MESH& m, TablesIndices& table_indices, TablesIndices& table_emb_edge)
	{
		if constexpr (mesh_traits<MESH>::dimension > 0)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			std::vector<uint32> i_e(thread_pool()->nb_workers(),0);
			parallel_foreach_cell(m, [&](Edge e) -> bool
			{
				uint32 worker_index = current_worker_index();
				foreach_incident_vertex(m, e, [&](Vertex v) -> bool
				{
					table_indices[worker_index].push_back(index_of(m, v));
					return true;
				});
				if (EMB)
					table_emb_edge[worker_index].push_back(index_of(m,e));
				else
					table_emb_edge[worker_index].push_back(i_e[worker_index]++);
				return true;
			});
		}
	}


	template <bool EMB, typename MESH>
	inline void init_triangles(const MESH& m, TablesIndices& table_indices, TablesIndices& table_emb_face)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;
			std::vector<std::vector<Vertex>> vvertices(thread_pool()->nb_workers());
			for (auto& v:vvertices)
				v.reserve(32u);
			std::vector<uint32> i_f(thread_pool()->nb_workers(),0);
			parallel_foreach_cell(m, [&](Face f) -> bool
			{
				uint32 worker_index = current_worker_index();
				if (EMB)
					i_f[worker_index] = index_of(m,f);
				auto& vertices = vvertices[worker_index];
				vertices.clear();
				incident_vertices(m, f, vertices);
				for (uint32 i = 1; i < vertices.size() - 1; ++i)
				{
					auto& tif = table_indices[worker_index];
					tif.push_back(index_of(m, vertices[0]));
					tif.push_back(index_of(m, vertices[i]));
					tif.push_back(index_of(m, vertices[i + 1]));
					table_emb_face[worker_index].push_back(i_f[worker_index]);
				}
				if (!EMB)
					i_f[worker_index]++;
				return true;
			});
		}
	}

	template <bool EMB,typename MESH>
	inline void init_ear_triangles(const MESH& m, TablesIndices& table_indices, TablesIndices& table_emb_face,
				   const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{
		if constexpr (mesh_traits<MESH>::dimension > 1)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Face = typename mesh_traits<MESH>::Face;

			std::vector<std::vector<Vertex>> vvertices(thread_pool()->nb_workers());
			for (auto& v:vvertices)
				v.reserve(32u);
			std::vector<uint32> i_f(thread_pool()->nb_workers(),0);
			parallel_foreach_cell(m, [&](Face f) -> bool
			{
				uint32 worker_index = current_worker_index();
				if (EMB)
					i_f[worker_index] = index_of(m,f);
				auto& tif = table_indices[worker_index];
				if (codegree(m,f)==3)
				{
					auto& vertices = vvertices[worker_index];
					vertices.clear();
					incident_vertices(m, f, vertices);

					tif.push_back(index_of(m, vertices[0]));
					tif.push_back(index_of(m, vertices[1]));
					tif.push_back(index_of(m, vertices[2]));
					table_emb_face[worker_index].push_back(i_f[worker_index]);
				}
				else
				{
					cgogn::geometry::append_ear_triangulation(m, f, position, tif,
					  [&] ()
					  {
						  table_emb_face[worker_index].push_back(i_f[worker_index]);
					  });
				}
				if (!EMB)
					i_f[worker_index]++;
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
//			bool emb = is_indexed<Face>(m);
//			foreach_cell(m, [&](Face f) -> bool
//			{
//				foreach_incident_vertex(m, f, [&](Vertex v) -> bool
//				{
//					table_indices.push_back(index_of(m, v));
//					table_indices.push_back(emb ? index_of(m,f) : last);
//					return true;
//				});
//				return true;
//				++last;
//			});
//		}
//	}

	template <bool EMB,typename MESH>
	inline void init_volumes(const MESH& m,
							 TablesIndices& table_indices_f,
							 TablesIndices& table_indices_e,
							 TablesIndices& table_indices_v,
							 TablesIndices& table_emb_vol,
			const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position)
	{

		if constexpr (mesh_traits<MESH>::dimension >= 2)
		{
			using Vertex = typename mesh_traits<MESH>::Vertex;
			using Edge = typename mesh_traits<MESH>::Edge;
			using Face = typename mesh_traits<MESH>::Face;
			using Volume = typename mesh_traits<MESH>::Volume;

			std::vector<std::vector<Vertex>> vvertices(thread_pool()->nb_workers());
			for (auto& v:vvertices)
				v.reserve(32u);
			std::vector<uint32> i_vol(thread_pool()->nb_workers(),0);
			parallel_foreach_cell(m, [&](Volume vol) -> bool
			{
				uint32 worker_index = current_worker_index();
				auto& ivol = i_vol[worker_index];

				if (EMB)
					ivol = index_of(m,vol);
				auto& vertices = vvertices[worker_index];
				foreach_incident_face(m, vol, [&](Face f)-> bool
				{
					auto& tif = table_indices_f[worker_index];
					if (codegree(m,f)==3)
					{
						vertices.clear();
						incident_vertices(m, f, vertices);	
						tif.push_back(index_of(m, vertices[0]));
						tif.push_back(index_of(m, vertices[1]));
						tif.push_back(index_of(m, vertices[2]));
						tif.push_back(ivol);
					}
					else
					{
						cgogn::geometry::append_ear_triangulation(m, f, position, tif,
								[&] () { tif.push_back(ivol); });
					}
					return true;
				});

				foreach_incident_edge(m, vol, [&](Edge e)-> bool
				{
					vertices.clear();
					incident_vertices(m, e, vertices);
					auto& ted = table_indices_e[worker_index];
					ted.push_back(index_of(m, vertices[0]));
					ted.push_back(index_of(m, vertices[1]));
					ted.push_back(ivol);
					return true;
				});

				foreach_incident_vertex(m, vol, [&](Vertex v)-> bool
				{
					auto& tv = table_indices_v[worker_index];
					tv.push_back(index_of(m, v));
					tv.push_back(ivol);
					return true;
				});

				table_emb_vol[worker_index].push_back(ivol);

				if (!EMB)
					ivol++;

				return true;
			});
		}
	}

public:

	template <typename MESH>
	inline void init_primitives(const MESH& m, DrawingType prim,
								const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* position = nullptr)
	{
		if (prim>=SIZE_BUFFER)
			prim = DrawingType(prim + POINTS - POINTS_TB);

		auto func_update_ebo = [&] (DrawingType pr, const TablesIndices& table ) -> void
		{
			uint32 total_size = 0;
			for (const auto& t : table)
				total_size += t.size();

			indices_buffers_uptodate_[pr] = true;
			if (total_size>0)
			{
				if (!indices_buffers_[pr]->is_created())
					indices_buffers_[pr]->create();

				indices_buffers_[pr]->allocate(total_size);
				uint32 beg = 0;
				for (const auto& t : table)
				{
					indices_buffers_[pr]->copy_data(beg, t.size(), t.data());
					beg += t.size();
				}
				indices_buffers_[pr]->set_name("EBO_"+primitives_names[pr]);
			}
		};

		auto func_update_ebo2 = [&] (DrawingType pr1, const TablesIndices& table1) -> void
		{
			uint32 total_size1 = 0;
			for (const auto& t : table1)
				total_size1 += t.size();

			indices_buffers_uptodate_[pr1] = true;
			if (total_size1>0)
			{
				if (!indices_buffers_[pr1]->is_created())
					indices_buffers_[pr1]->create();

				indices_buffers_[pr1]->allocate(total_size1);
				indices_buffers_[pr1]->bind();
				uint32* ptr = indices_buffers_[pr1]->lock_pointer();
				uint32 beg = 0;
				for (const auto& t : table1)
				{
					for (uint32 i : t)
						*ptr++ = i+beg;
					beg += t.empty()?0:t.back()+1;
				}
				indices_buffers_[pr1]->set_name("EBO_"+primitives_names[pr1]);
				indices_buffers_[pr1]->release_pointer();
				indices_buffers_[pr1]->release();
			}
		};


		auto func_update_ebo3 = [&] (DrawingType pr1, const TablesIndices& table1,
				const TablesIndices& table2, uint32 interv) -> void
		{
			uint32 total_size1 = 0;
			for (const auto& t : table1)
				total_size1 += t.size();

			indices_buffers_uptodate_[pr1] = true;
			if (total_size1 > 0)
			{
				if (!indices_buffers_[pr1]->is_created())
					indices_buffers_[pr1]->create();

				indices_buffers_[pr1]->allocate(total_size1);

				indices_buffers_[pr1]->bind();
				uint32* ptr1 = indices_buffers_[pr1]->lock_pointer();

				uint32 beg = 0;
				uint32 nb = table1.size();
				for (int j = 0; j<nb; ++j)
				{
					const auto& t1 = table1[j];
					uint32 sz = t1.size();
					for (uint32 k = 0; k<sz; ++k)
						*ptr1++ = (k%interv == interv-1) ? t1[k]+beg : t1[k];

					beg += table2[j].empty() ? 0 : table2[j].back()+1;
				}

				indices_buffers_[pr1]->set_name("EBO_"+primitives_names[pr1]);
				indices_buffers_[pr1]->release_pointer();
				indices_buffers_[pr1]->release();

			}
		};



		auto start_timer = std::chrono::high_resolution_clock::now();

		uint32 nbw = thread_pool()->nb_workers();

		TablesIndices table_indices(nbw);
		for (auto& t : table_indices)
			t.reserve(1024u);

		TablesIndices table_indices_emb(nbw);
		for (auto& t : table_indices_emb)
			t.reserve(1024u);

		TablesIndices table_indices_e(nbw);
		for (auto& t : table_indices_e)
			t.reserve(1024u);

		TablesIndices table_indices_v(nbw);
		for (auto& t : table_indices_v)
			t.reserve(1024u);



		switch (prim)
		{
		case POINTS:
			init_points(m, table_indices);
			func_update_ebo(POINTS,table_indices);
			break;
		case LINES:
		case INDEX_EDGES:
			if (is_indexed<typename mesh_traits<MESH>::Edge>(m))
			{
				init_lines<true>(m, table_indices,table_indices_emb);
				func_update_ebo(INDEX_EDGES,table_indices_emb);
			}
			else
			{
				init_lines<false>(m, table_indices,table_indices_emb);
				func_update_ebo2(INDEX_EDGES,table_indices_emb);
			}
			func_update_ebo(LINES,table_indices);
			break;
		case TRIANGLES:
		case INDEX_FACES:
			if constexpr (mesh_traits<MESH>::dimension >= 2)
			{
				if (is_indexed<typename mesh_traits<MESH>::Face>(m))
				{
					if (position == nullptr)
						init_triangles<true>(m, table_indices, table_indices_emb);
					else
						init_ear_triangles<true>(m, table_indices, table_indices_emb, position);
					func_update_ebo(INDEX_FACES,table_indices_emb);
				}
				else
				{
					if (position == nullptr)
						init_triangles<false>(m, table_indices, table_indices_emb);
					else
						init_ear_triangles<false>(m, table_indices, table_indices_emb, position);
					func_update_ebo2(INDEX_FACES,table_indices_emb);
				}
				func_update_ebo(TRIANGLES,table_indices);
			}
			break;

		case VOLUMES_VERTICES:
		case VOLUMES_EDGES:
		case VOLUMES_FACES:
		case INDEX_VOLUMES:
			if constexpr (mesh_traits<MESH>::dimension >= 2)
			{
				if (is_indexed<typename mesh_traits<MESH>::Volume>(m))
				{
					init_volumes<true>(m, table_indices, table_indices_e,table_indices_v, table_indices_emb, position);
					func_update_ebo(VOLUMES_FACES,table_indices);
					func_update_ebo(VOLUMES_EDGES,table_indices_e);
					func_update_ebo(VOLUMES_VERTICES,table_indices_v);
					func_update_ebo2(INDEX_VOLUMES, table_indices_emb);
				}
				else
				{
					init_volumes<false>(m, table_indices, table_indices_e,table_indices_v, table_indices_emb, position);
					func_update_ebo3(VOLUMES_FACES,table_indices,table_indices_emb,4);
					func_update_ebo3(VOLUMES_EDGES,table_indices_e,table_indices_emb,3);
					func_update_ebo3(VOLUMES_VERTICES,table_indices_v,table_indices_emb,2);
					func_update_ebo2(INDEX_VOLUMES, table_indices_emb);
				}
			}
			break;

		default:
			break;
		}

		auto end_timer = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = end_timer-start_timer;
		std::cout << "init primitive "<<prim<< " in "<< elapsed_seconds.count() << std::endl;

	}

	void draw(DrawingType prim, GLint binding_point=10);

	inline void bind_ebo_tb(DrawingType prim, GLint binding_point)
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
