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

#ifndef CGOGN_RENDERING_VOLUME_DRAWER_H_
#define CGOGN_RENDERING_VOLUME_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shaders/shader_explode_volumes.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_color.h>
#include <cgogn/rendering/shaders/shader_explode_volumes_line.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

using namespace geometry;

namespace rendering
{

/**
 * @brief Rendering of volumes
 *
 * Typical usage:
 *
 *  std::unique_ptr<cgogn::rendering::VolumeDrawer> volu_;	// can be shared between contexts
 *  std::unique_ptr<cgogn::rendering::VolumeDrawer::Renderer> volu_rend_; // one by context,
 *
 * init:
 *  volu_ = std::make_unique<cgogn::rendering::VolumeDrawer>();
 *  volu_rend_ = volu_->generate_renderer();
 *  volu_->update_face(map_, vertex_position_);
 *  volu_->update_edge(map_, vertex_position_);
 *
 * draw:
 *  volu_rend_->set_explode_volume(0.9);
 *  volu_rend_->draw_faces(proj, view);
 *  volu_rend_->draw_edges(proj, view);
 *
 */
class CGOGN_RENDERING_EXPORT VolumeDrawerGen
{
protected:
	using Vec3f = std::array<float32, 3>;

	std::unique_ptr<VBO> vbo_pos_;
	std::unique_ptr<VBO> vbo_col_;

	GLColor face_color_;

	std::unique_ptr<VBO> vbo_pos2_;
	GLColor edge_color_;

	float32 shrink_v_;
	float32 shrink_f_;

	void init_with_color();
	void init_without_color();
	void init_edge();

public:
	class CGOGN_RENDERING_EXPORT Renderer
	{
		friend class VolumeDrawerGen;

		std::unique_ptr<ShaderExplodeVolumes::Param> param_expl_vol_;
		std::unique_ptr<ShaderExplodeVolumesColor::Param> param_expl_vol_col_;
		std::unique_ptr<ShaderExplodeVolumesLine::Param> param_expl_vol_line_;
		VolumeDrawerGen* volume_drawer_data_;

		Renderer(VolumeDrawerGen* tr);

	public:
		~Renderer();
		void draw_faces(const GLMat4& projection, const GLMat4& modelview);
		void draw_edges(const GLMat4& projection, const GLMat4& modelview);
		void set_explode_volume(float32 x);
		void set_face_color(const GLColor& rgb);
		void set_edge_color(const GLColor& rgb);
		void set_clipping_plane(const GLVec4& pl);
		void set_clipping_plane2(const GLVec4& pl);
		void set_thick_clipping_plane(const GLVec4& p, float32 th);
	};

	using Self = VolumeDrawerGen;

	/**
	 * constructor, init all buffers (data and OpenGL) and shader
	 * @Warning need OpenGL context
	 */
	VolumeDrawerGen(bool with_color_per_face);

	/**
	 * release buffers and shader
	 */
	virtual ~VolumeDrawerGen();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(VolumeDrawerGen);

	/**
	 * @brief generate a renderer (one per context)
	 * @return pointer on renderer
	 */
	inline std::unique_ptr<Renderer> generate_renderer()
	{
		return std::unique_ptr<Renderer>(new Renderer(this));
	}

	template <typename MESH>
	void update_edge(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* position);
};

template <typename MESH>
void VolumeDrawerGen::update_edge(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Volume = typename mesh_traits<MESH>::Volume;

	std::vector<Vec3f> out_pos;
	out_pos.reserve(1024 * 1024);

	std::vector<uint32> ear_indices;
	ear_indices.reserve(256);

	foreach_cell(m, [&](Volume v) {
		Vec3 CV = geometry::centroid(m, v, position);
		foreach_incident_edge(m, v, [&](Edge e) -> bool {
			auto vs = incident_vertices(m, e); // WARNING PERFORMANCE ISSUE, SPECIAL PAIR VERSION
			const Vec3& P1 = value<Vec3>(m, position, vs[0]);
			const Vec3& P2 = value<Vec3>(m, position, vs[1]);
			out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
			out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
			out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
			return true;
		});
		return true;
	});

	uint32 nbvec = uint32(uint32(out_pos.size()));
	vbo_pos2_->allocate(nbvec, 3);
	vbo_pos2_->bind();
	vbo_pos2_->copy_data(0, nbvec * 12, out_pos[0].data());
	vbo_pos2_->release();
}

template <bool CPV>
class VolumeDrawerTpl : public VolumeDrawerGen
{
};

template <>
class CGOGN_RENDERING_EXPORT VolumeDrawerTpl<false> : public VolumeDrawerGen
{
public:
	VolumeDrawerTpl() : VolumeDrawerGen(false)
	{
	}

	~VolumeDrawerTpl() override;

	template <typename MESH>
	void update_face(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* position)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Face = typename mesh_traits<MESH>::Face;
		using Volume = typename mesh_traits<MESH>::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		foreach_cell(m, [&](Volume v) -> bool {
			Vec3 CV = geometry::centroid(m, v, position);
			foreach_incident_face(m, v, [&](Face f) -> bool {
				if (codegree(m, f) < 3)
				{
					auto vs = incident_vertices(m, f); // WARNING PERFORMANCE ISSUE, SPECIAL TRIPLET VERSION
					const Vec3& P1 = value<Vec3>(m, position, vs[0]);
					const Vec3& P2 = value<Vec3>(m, position, vs[1]);
					const Vec3& P3 = value<Vec3>(m, position, vs[2]);
					out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
					out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
					out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
					out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
				}
				else
				{
					ear_indices.clear();
					cgogn::geometry::append_ear_triangulation(m, f, position, ear_indices);
					for (std::size_t i = 0; i < uint32(ear_indices.size()); i += 3)
					{
						const Vec3& P1 = (*position)[ear_indices[i]];
						const Vec3& P2 = (*position)[ear_indices[i + 1]];
						const Vec3& P3 = (*position)[ear_indices[i + 2]];
						out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
						out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
						out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
						out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
					}
				}
				return true;
			});
			return true;
		});

		uint32 nbvec = uint32(uint32(out_pos.size()));

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();
	}
};

template <>
class CGOGN_RENDERING_EXPORT VolumeDrawerTpl<true> : public VolumeDrawerGen
{
public:
	VolumeDrawerTpl() : VolumeDrawerGen(true)
	{
	}

	~VolumeDrawerTpl() override;

	template <typename MESH>
	void update_face_color_vertex(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* position,
								  const typename mesh_traits<MESH>::template Attribute<Vec3>* color)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Face = typename mesh_traits<MESH>::Face;
		using Volume = typename mesh_traits<MESH>::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<Vec3f> out_color;
		out_color.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		foreach_cell(m, [&](Volume v) -> bool {
			Vec3 CV = geometry::centroid(m, v, position);
			foreach_incident_face(m, v, [&](Face f) -> bool {
				if (m.has_codegree(f, 3))
				{
					auto vs = incident_vertices(m, f); // WARNING PERFORMANCE ISSUE, SPECIAL TRIPLET VERSION
					const Vec3& P1 = value<Vec3>(m, position, vs[0]);
					const Vec3& C1 = value<Vec3>(m, color, vs[0]);
					const Vec3& P2 = value<Vec3>(m, position, vs[1]);
					const Vec3& C2 = value<Vec3>(m, color, vs[1]);
					const Vec3& P3 = value<Vec3>(m, position, vs[2]);
					const Vec3& C3 = value<Vec3>(m, color, vs[2]);
					out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
					out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
					out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
					out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
					out_color.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
					out_color.push_back({float32(C1[0]), float32(C1[1]), float32(C1[2])});
					out_color.push_back({float32(C2[0]), float32(C2[1]), float32(C2[2])});
					out_color.push_back({float32(C3[0]), float32(C3[1]), float32(C3[2])});
				}
				else
				{
					ear_indices.clear();
					cgogn::geometry::append_ear_triangulation(m, f, position, ear_indices);
					for (std::size_t i = 0; i < uint32(ear_indices.size()); i += 3)
					{
						const Vec3& P1 = (*position)[ear_indices[i]];
						const Vec3& C1 = (*color)[ear_indices[i]];
						const Vec3& P2 = (*position)[ear_indices[i + 1]];
						const Vec3& C2 = (*color)[ear_indices[i + 1]];
						const Vec3& P3 = (*position)[ear_indices[i + 2]];
						const Vec3& C3 = (*color)[ear_indices[i + 2]];
						out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
						out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
						out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
						out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
						out_color.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
						out_color.push_back({float32(C1[0]), float32(C1[1]), float32(C1[2])});
						out_color.push_back({float32(C2[0]), float32(C2[1]), float32(C2[2])});
						out_color.push_back({float32(C3[0]), float32(C3[1]), float32(C3[2])});
					}
				}
				return true;
			});
			return true;
		});

		std::size_t nbvec = uint32(out_pos.size());

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();

		vbo_col_->allocate(nbvec, 3);
		vbo_col_->bind();
		vbo_col_->copy_data(0, nbvec * 12, out_color[0].data());
		vbo_col_->release();
	}

	template <typename MESH>
	void update_face_color_volume(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* position,
								  const typename mesh_traits<MESH>::template Attribute<Vec3>* color)
	{
		//		using Vertex = typename mesh_traits<MESH>::Vertex;
		using Face = typename mesh_traits<MESH>::Face;
		using Volume = typename mesh_traits<MESH>::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<Vec3f> out_color;
		out_color.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		foreach_cell(m, [&](Volume v) -> bool {
			Vec3 CV = geometry::centroid(m, v, position);
			const Vec3& C = value<Vec3>(m, color, v);
			m.foreach_incident_face(v, [&](Face f) {
				if (m.has_codegree(f, 3))
				{
					auto vs = incident_vertices(m, f); // WARNING PERFORMANCE ISSUE, SPECIAL TRIPLET VERSION
					const Vec3& P1 = value<Vec3>(m, position, vs[0]);
					const Vec3& P2 = value<Vec3>(m, position, vs[1]);
					const Vec3& P3 = value<Vec3>(m, position, vs[2]);
					out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
					out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
					out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
					out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
					out_color.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
					out_color.push_back({float32(C[0]), float32(C[1]), float32(C[2])});
					out_color.push_back({float32(C[0]), float32(C[1]), float32(C[2])});
					out_color.push_back({float32(C[0]), float32(C[1]), float32(C[2])});
				}
				else
				{
					ear_indices.clear();
					cgogn::geometry::append_ear_triangulation(m, f, position, ear_indices);
					for (std::size_t i = 0; i < uint32(ear_indices.size()); i += 3)
					{
						const Vec3& P1 = value<Vec3>(m, position, ear_indices[i]);
						const Vec3& P2 = value<Vec3>(m, position, ear_indices[i + 1]);
						const Vec3& P3 = value<Vec3>(m, position, ear_indices[i + 2]);
						out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
						out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
						out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
						out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
						out_color.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
						out_color.push_back({float32(C[0]), float32(C[1]), float32(C[2])});
						out_color.push_back({float32(C[0]), float32(C[1]), float32(C[2])});
						out_color.push_back({float32(C[0]), float32(C[1]), float32(C[2])});
					}
				}
				return true;
			});
			return true;
		});

		std::size_t nbvec = uint32(out_pos.size());

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();

		vbo_col_->allocate(nbvec, 3);
		vbo_col_->bind();
		vbo_col_->copy_data(0, nbvec * 12, out_color[0].data());
		vbo_col_->release();
	}
};

using VolumeDrawer = VolumeDrawerTpl<false>;
using VolumeDrawerColor = VolumeDrawerTpl<true>;

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_VOLUME_DRAWER_H_
