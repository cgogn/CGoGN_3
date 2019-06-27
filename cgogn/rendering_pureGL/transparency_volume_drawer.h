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

#ifndef CGOGN_RENDERING_TRANSP_VOLUME_DRAWER_H_
#define CGOGN_RENDERING_TRANSP_VOLUME_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shaders/vbo.h>
#include <cgogn/rendering/transparency_shaders/shader_transparent_volumes.h>

#include <cgogn/geometry/types/geometry_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/ear_triangulation.h>

#include <QOpenGLFunctions_3_3_Core>
#include <GLColor>
#include <QOpenGLFramebufferObject>

namespace cgogn
{

namespace rendering
{


class CGOGN_RENDERING_EXPORT VolumeTransparencyDrawer
{
protected:

	using Vec3f = std::array<float32, 3>;
	std::unique_ptr<VBO> vbo_pos_;
	GLColor face_color_;
	float32 shrink_v_;

public:
	class CGOGN_RENDERING_EXPORT Renderer
	{
		friend class VolumeTransparencyDrawer;
		std::unique_ptr<ShaderTransparentVolumes::Param> param_transp_vol_;
		VolumeTransparencyDrawer* volume_drawer_data_;

		Renderer(VolumeTransparencyDrawer* tr);

	public:
		~Renderer();

		void draw_faces(const QMatrix4x4& projection, const QMatrix4x4& modelview);

		void set_explode_volume(float32 x);

		void set_color(const GLColor& rgb);

		void set_clipping_plane(const GLVec4& pl);

		void set_clipping_plane2(const GLVec4& pl);

		void set_thick_clipping_plane(const GLVec4& p, float32 th);

		void set_back_face_culling(bool cull);

		void set_lighted(bool lighted);

	};

	using Self = VolumeTransparencyDrawer;

	/**
	 * constructor, init all buffers (data and OpenGL) and shader
	 * @Warning need OpenGL context
	 */
	VolumeTransparencyDrawer();


	CGOGN_NOT_COPYABLE_NOR_MOVABLE(VolumeTransparencyDrawer);

	/**
	 * @brief generate a renderer (one per context)
	 * @return pointer on renderer
	 */
	inline std::unique_ptr<Renderer> generate_renderer()
	{
		return std::unique_ptr<Renderer>(new Renderer(this));
	}

	template <typename MAP, typename VERTEX_ATTR>
	void update_face(const MAP& m, const VERTEX_ATTR& position);
};





template <typename MAP, typename VERTEX_ATTR>
void VolumeTransparencyDrawer::update_face(const MAP& m, const VERTEX_ATTR& position)
{
	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value,"position must be a vertex attribute");

	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
	using Vertex = typename MAP::Vertex;
	using Face = typename MAP::Face;
	using Volume = typename MAP::Volume;

	std::vector<Vec3f> out_pos;
	out_pos.reserve(1024 * 1024);

	std::vector<uint32> ear_indices;
	ear_indices.reserve(256);

	m.foreach_cell([&] (Volume v)
	{
		VEC3 CV = geometry::centroid(m, v, position);
		m.foreach_incident_face(v, [&] (Face f)
		{
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
				for(std::size_t i = 0; i < ear_indices.size(); i += 3)
				{
					const VEC3& P1 = position[ear_indices[i]];
					const VEC3& P2 = position[ear_indices[i+1]];
					const VEC3& P3 = position[ear_indices[i+2]];
					out_pos.push_back({float32(CV[0]), float32(CV[1]), float32(CV[2])});
					out_pos.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
					out_pos.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});
					out_pos.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
				}
			}
		});
	});

	uint32 nbvec = uint32(out_pos.size());

	vbo_pos_->allocate(nbvec, 3);
	vbo_pos_->bind();
	vbo_pos_->copy_data(0, nbvec * 12, out_pos[0].data());
	vbo_pos_->release();
}



} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_TRANSP_VOLUME_DRAWER_H_
