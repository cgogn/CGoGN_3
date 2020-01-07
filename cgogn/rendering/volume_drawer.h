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
#include <cgogn/geometry/types/geometry_traits.h>

namespace cgogn
{

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
 *  volu_ = cgogn::make_unique<cgogn::rendering::VolumeDrawer>();
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

	template <typename MAP, typename MASK, typename VERTEX_ATTR>
	void update_edge(const MAP& m, const MASK& mask, const VERTEX_ATTR& position);

	template <typename MAP, typename VERTEX_ATTR>
	void update_edge(const MAP& m, const VERTEX_ATTR& position);
};

template <typename MAP, typename MASK, typename VERTEX_ATTR>
void VolumeDrawerGen::update_edge(const MAP& m, const MASK& mask, const VERTEX_ATTR& position)
{
	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position must be a vertex attribute");

	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
	using Vertex = typename MAP::Vertex;
	using Edge = typename MAP::Edge;
	using Volume = typename MAP::Volume;

	std::vector<Vec3f> out_pos;
	out_pos.reserve(1024 * 1024);

	std::vector<uint32> ear_indices;
	ear_indices.reserve(256);

	m.foreach_cell(
		[&](Volume v) {
			VEC3 CV = geometry::centroid(m, v, position);
			m.foreach_incident_edge(v, [&](Edge e) {
				const VEC3& P1 = position[Vertex(e.dart)];
				const VEC3& P2 = position[Vertex(m.phi1(e.dart))];
				out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
				out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
				out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
			});
		},
		mask);

	uint32 nbvec = uint32(out_pos.size());
	vbo_pos2_->allocate(nbvec, 3);
	vbo_pos2_->bind();
	vbo_pos2_->copy_data(0, nbvec * 12, out_pos[0].data());
	vbo_pos2_->release();
}

template <typename MAP, typename VERTEX_ATTR>
void VolumeDrawerGen::update_edge(const MAP& m, const VERTEX_ATTR& position)
{
	update_edge(m, AllCellsFilter(), position);
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

	template <typename MAP, typename MASK, typename VERTEX_ATTR>
	void update_face(const MAP& m, const MASK& mask, const VERTEX_ATTR& position)
	{
		static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position must be a vertex attribute");

		using VEC3 = InsideTypeOf<VERTEX_ATTR>;
		using Vertex = typename MAP::Vertex;
		using Face = typename MAP::Face;
		using Volume = typename MAP::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		m.foreach_cell(
			[&](Volume v) {
				VEC3 CV = geometry::centroid(m, v, position);
				m.foreach_incident_face(v, [&](Face f) {
					if (m.has_codegree(f, 3))
					{
						const VEC3& P1 = position[Vertex(f.dart)];
						const VEC3& P2 = position[Vertex(m.phi1(f.dart))];
						const VEC3& P3 = position[Vertex(m.phi1(m.phi1(f.dart)))];
						out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
						out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
						out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
						out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
					}
					else
					{
						ear_indices.clear();
						cgogn::geometry::append_ear_triangulation(m, f, position, ear_indices);
						for (std::size_t i = 0; i < ear_indices.size(); i += 3)
						{
							const VEC3& P1 = position[ear_indices[i]];
							const VEC3& P2 = position[ear_indices[i + 1]];
							const VEC3& P3 = position[ear_indices[i + 2]];
							out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
							out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
							out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
							out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
						}
					}
				});
			},
			mask);

		uint32 nbvec = uint32(out_pos.size());

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();
	}

	template <typename MAP, typename VERTEX_ATTR>
	void update_face(const MAP& m, const VERTEX_ATTR& position)
	{
		update_face(m, AllCellsFilter(), position);
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

	template <typename MAP, typename MASK, typename VERTEX_ATTR>
	void update_face(const MAP& m, const MASK& mask, const VERTEX_ATTR& position, const VERTEX_ATTR& color)
	{
		static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position must be a vertex attribute");

		using VEC3 = InsideTypeOf<VERTEX_ATTR>;
		using Vertex = typename MAP::Vertex;
		using Face = typename MAP::Face;
		using Volume = typename MAP::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<Vec3f> out_color;
		out_color.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		m.foreach_cell(
			[&](Volume v) {
				VEC3 CV = geometry::centroid(m, v, position);
				m.foreach_incident_face(v, [&](Face f) {
					if (m.has_codegree(f, 3))
					{
						Dart d = f.dart;
						const VEC3& P1 = position[Vertex(d)];
						const VEC3& C1 = color[Vertex(d)];
						d = m.phi1(d);
						const VEC3& P2 = position[Vertex(d)];
						const VEC3& C2 = color[Vertex(d)];
						d = m.phi1(d);
						const VEC3& P3 = position[Vertex(d)];
						const VEC3& C3 = color[Vertex(d)];
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
						for (std::size_t i = 0; i < ear_indices.size(); i += 3)
						{
							const VEC3& P1 = position[ear_indices[i]];
							const VEC3& C1 = color[ear_indices[i]];
							const VEC3& P2 = position[ear_indices[i + 1]];
							const VEC3& C2 = color[ear_indices[i + 1]];
							const VEC3& P3 = position[ear_indices[i + 2]];
							const VEC3& C3 = color[ear_indices[i + 2]];
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
				});
			},
			mask);

		std::size_t nbvec = out_pos.size();

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();

		vbo_col_->allocate(nbvec, 3);
		vbo_col_->bind();
		vbo_col_->copy_data(0, nbvec * 12, out_color[0].data());
		vbo_col_->release();
	}

	template <typename MAP, typename VERTEX_ATTR>
	void update_face(const MAP& m, const VERTEX_ATTR& position, const VERTEX_ATTR& color)
	{
		update_face(m, AllCellsFilter(), position, color);
	}

	template <typename MAP, typename MASK, typename VERTEX_ATTR, typename FACE_ATTR,
			  typename std::enable_if<is_orbit_of<FACE_ATTR, MAP::Face::ORBIT>::value, void>::type* = nullptr>
	void update_face(const MAP& m, const MASK& mask, const VERTEX_ATTR& position, const FACE_ATTR& color)
	{
		static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position must be a vertex attribute");

		using VEC3 = InsideTypeOf<VERTEX_ATTR>;
		using Vertex = typename MAP::Vertex;
		using Face = typename MAP::Face;
		using Volume = typename MAP::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<Vec3f> out_color;
		out_color.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		m.foreach_cell(
			[&](Volume v) {
				VEC3 CV = geometry::centroid(m, v, position);
				m.foreach_incident_face(v, [&](Face f) {
					const VEC3& C = color[f];

					if (m.has_codegree(f, 3))
					{
						Dart d = f.dart;
						const VEC3& P1 = position[Vertex(d)];
						d = m.phi1(d);
						const VEC3& P2 = position[Vertex(d)];
						d = m.phi1(d);
						const VEC3& P3 = position[Vertex(d)];
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
						for (std::size_t i = 0; i < ear_indices.size(); i += 3)
						{
							const VEC3& P1 = position[ear_indices[i]];
							const VEC3& P2 = position[ear_indices[i + 1]];
							const VEC3& P3 = position[ear_indices[i + 2]];
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
				});
			},
			mask);

		std::size_t nbvec = out_pos.size();

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();

		vbo_col_->allocate(nbvec, 3);
		vbo_col_->bind();
		vbo_col_->copy_data(0, nbvec * 12, out_color[0].data());
		vbo_col_->release();
	}

	template <typename MAP, typename VERTEX_ATTR, typename FACE_ATTR,
			  typename std::enable_if<is_orbit_of<FACE_ATTR, MAP::Face::ORBIT>::value, void>::type* = nullptr>
	void update_face(const MAP& m, const VERTEX_ATTR& position, const FACE_ATTR& color)
	{
		update_face(m, AllCellsFilter(), position, color);
	}

	template <typename MAP, typename MASK, typename VERTEX_ATTR, typename VOLUME_ATTR,
			  typename std::enable_if<is_orbit_of<VOLUME_ATTR, MAP::Volume::ORBIT>::value, void>::type* = nullptr>
	void update_face(const MAP& m, const MASK& mask, const VERTEX_ATTR& position, const VOLUME_ATTR& color)
	{
		static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value, "position must be a vertex attribute");

		using VEC3 = InsideTypeOf<VERTEX_ATTR>;
		using Vertex = typename MAP::Vertex;
		using Face = typename MAP::Face;
		using Volume = typename MAP::Volume;

		std::vector<Vec3f> out_pos;
		out_pos.reserve(1024 * 1024);

		std::vector<Vec3f> out_color;
		out_color.reserve(1024 * 1024);

		std::vector<uint32> ear_indices;
		ear_indices.reserve(256);

		m.foreach_cell(
			[&](Volume v) {
				VEC3 CV = geometry::centroid(m, v, position);
				const VEC3& C = color[v];
				m.foreach_incident_face(v, [&](Face f) {
					if (m.has_codegree(f, 3))
					{
						Dart d = f.dart;
						const VEC3& P1 = position[Vertex(d)];
						d = m.phi1(d);
						const VEC3& P2 = position[Vertex(d)];
						d = m.phi1(d);
						const VEC3& P3 = position[Vertex(d)];
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
						for (std::size_t i = 0; i < ear_indices.size(); i += 3)
						{
							const VEC3& P1 = position[ear_indices[i]];
							const VEC3& P2 = position[ear_indices[i + 1]];
							const VEC3& P3 = position[ear_indices[i + 2]];
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
				});
			},
			mask);

		std::size_t nbvec = out_pos.size();

		vbo_pos_->allocate(nbvec, 3);
		vbo_pos_->bind();
		vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
		vbo_pos_->release();

		vbo_col_->allocate(nbvec, 3);
		vbo_col_->bind();
		vbo_col_->copy_data(0, nbvec * 12, out_color[0].data());
		vbo_col_->release();
	}

	template <typename MAP, typename VERTEX_ATTR, typename VOLUME_ATTR,
			  typename std::enable_if<is_orbit_of<VOLUME_ATTR, MAP::Volume::ORBIT>::value, void>::type* = nullptr>
	void update_face(const MAP& m, const VERTEX_ATTR& position, const VOLUME_ATTR& color)
	{
		update_face(m, AllCellsFilter(), position, color);
	}
};

using VolumeDrawer = VolumeDrawerTpl<false>;
using VolumeDrawerColor = VolumeDrawerTpl<true>;

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_VOLUME_DRAWER_H_
