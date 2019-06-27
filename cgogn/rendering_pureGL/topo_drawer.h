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

#ifndef CGOGN_RENDERING_TOPO_DRAWER_H_
#define CGOGN_RENDERING_TOPO_DRAWER_H_

#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>

#include <cgogn/rendering_pureGL/shaders/shader_simple_color.h>
#include <cgogn/rendering_pureGL/shaders/shader_bold_line.h>
#include <cgogn/rendering_pureGL/shaders/shader_bold_line_color.h>
#include <cgogn/rendering_pureGL/shaders/shader_round_point_color.h>

#include <cgogn/geometry/algos/centroid.h>
//#include <cgogn/geometry/functions/distance.h>
#include <array>

namespace cgogn
{

namespace rendering_pgl
{
/**
 * @brief Rendering of the topology
 *
 * Typical usage:
 *
 *  std::unique_ptr<cgogn::rendering::TopoDrawer> topo_;	// can be shared between contexts
 *  std::unique_ptr<cgogn::rendering::TopoDrawer::Renderer> topo_rend_; // one by context,
 *
 * init:
 *  topo_ = cgogn::make_unique<cgogn::rendering::TopoDrawer>();
 *  topo_rend_ = topo_->generate_renderer();
 *  topo_->update(map_,vertex_position_);
 *
 * draw:
 *  topo_rend_->draw(proj,view,this);
 *
 */
class CGOGN_RENDERING_PUREGL_EXPORT TopoDrawer
{
	using Vec3f = Eigen::Vector3f;

protected:

	std::unique_ptr<VBO> vbo_darts_;
	std::unique_ptr<VBO> vbo_relations_;
	std::unique_ptr<VBO> vbo_color_darts_;

	GLColor dart_color_;
	GLColor phi2_color_;
	GLColor phi3_color_;

	float32 shrink_v_;
	float32 shrink_f_;
	float32 shrink_e_;

	std::vector<Vec3f> darts_pos_;
	std::vector<Dart> darts_id_;

public:
	template <typename MAP, typename VERTEX_ATTR>
	auto update(const MAP& m, const VERTEX_ATTR& position)
	-> typename std::enable_if<MAP::DIMENSION == 2, void>::type;

	template <typename MAP, typename VERTEX_ATTR>
	auto update(const MAP& m, const VERTEX_ATTR& position)
	-> typename std::enable_if<MAP::DIMENSION == 3, void>::type;

	class CGOGN_RENDERING_PUREGL_EXPORT Renderer
	{
		friend class TopoDrawer;

		std::unique_ptr<ShaderBoldLineColor::Param> param_bl_;
		std::unique_ptr<ShaderBoldLine::Param> param_bl2_;
		std::unique_ptr<ShaderRoundPointColor::Param> param_rp_;
		TopoDrawer* topo_drawer_data_;

		Renderer(TopoDrawer* tr);

	public:

		~Renderer();

		/**
		 * @brief draw
		 * @param projection projection matrix
		 * @param modelview model-view matrix
		 * @param with_blending
		 */
		void draw(const GLMat4& projection, const GLMat4& modelview, bool with_blending = true);

		void set_clipping_plane(const GLVec4& p);

		void set_clipping_plane2(const GLVec4& p);

		void set_thick_clipping_plane(const GLVec4& p, float32 th);
	};

	using Self = TopoDrawer;

	/**
	 * constructor, init all buffers (data and OpenGL) and shader
	 * @Warning need OpenGL context
	 */
	TopoDrawer();

	/**
	 * release buffers and shader
	 */
	~TopoDrawer();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(TopoDrawer);

	/**
	 * @brief generate a renderer (one per context)
	 * @return pointer on renderer
	 */
	inline std::unique_ptr<Renderer> generate_renderer()
	{
		return std::unique_ptr<Renderer>(new Renderer(this));
	}

	inline void set_explode_volume(float32 x) { shrink_v_ = x; }

	inline void set_explode_face(float32 x) { shrink_f_ = x; }

	inline void set_explode_edge(float32 x) { shrink_e_ = x; }

	/**
	 * @brief update colors of darts
	 * @warning positions must be updated before
	 * @param color color attribute (of any orbit)
	 */
	template <typename ATTR>
	void update_colors(const ATTR& color);

	/**
	 * @brief update color of one dart
	 * @warning O(n) perf.
	 * @param d the dart
	 * @param rgb the color
	 */
	template <typename RGB>
	void update_color(Dart d, const RGB& rgb);

	/**
	 * @brief update color of one dart
	 * @warning O(n) perf.
	 * @param d the dart
	 * @param rgb the color
	 */
	void update_color(Dart d, const GLColor& rgb);

	/**
	 * @brief pick the closest dart to a given ray
	 * @param A ray first point
	 * @param B ray second point
	 * @param plane picking plane (use std::array<int,4>{{0,0,0,0}} for none)
	 * @param dp1 first point of selected dart
	 * @param dp2 second point of selected dart
	 * @return the selected dart
	 */
	template <typename VEC3, typename VEC4>
	Dart pick(const VEC3& A, const VEC3& B, const VEC4& plane, VEC3* dp1=nullptr, VEC3* dp2=nullptr);

	template <typename VEC3, typename VEC4>
	Dart pick(const VEC3& A, const VEC3& B, const VEC4& plane1, const VEC4& plane2, VEC3* dp1=nullptr, VEC3* dp2=nullptr);

	template <typename VEC3, typename VEC4>
	Dart pick(const VEC3& A, const VEC3& B, const VEC4& plane, float32 thickness, VEC3* dp1=nullptr, VEC3* dp2=nullptr);

};


template <typename MAP, typename VERTEX_ATTR>
auto TopoDrawer::update(const MAP& m, const VERTEX_ATTR& position)
-> typename std::enable_if<MAP::DIMENSION == 2, void>::type
{
//	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value,"position must be a vertex attribute");

//	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
	using Vertex = typename MAP::Vertex;
	using Face = typename MAP::Face;
	using Scalar = typename geometry::vector_traits<VEC3>::Scalar;

	Scalar opp_shrink_e = Scalar(1.0f - shrink_e_);
	Scalar opp_shrink_f = Scalar(1.0f - shrink_f_);

	darts_pos_.clear();
	darts_pos_.reserve(1024 * 1024);

	darts_id_.clear();
	darts_id_.reserve(1024 * 1024);

	std::vector<Vec3f> out_pos2;
	out_pos2.reserve(1024 * 1024);

	std::vector<VEC3> local_vertices;
	local_vertices.reserve(256);

	std::vector<Dart> local_darts;
	local_darts.reserve(256);

	m.foreach_cell([&] (Face f)
	{
		local_vertices.clear();
		local_darts.clear();
		VEC3 center;
		center.setZero();
		uint32 count = 0u;
		m.foreach_incident_vertex(f, [&] (Vertex v)
		{
			local_vertices.push_back(position[v]);
			local_darts.push_back(v.dart);
			center += position[v];
			count++;
		});
		center /= Scalar(count);

		// phi2 mid-edge: N -> 2N-1
		for (uint32 i = 0; i < count; ++i)
			local_vertices.push_back((local_vertices[i]+local_vertices[(i+1)%count])/Scalar(2.0));

		// dart round point: 0 -> N-1
		for (uint32 i = 0; i < count; ++i)
			local_vertices[i] = local_vertices[i] * Scalar(shrink_f_) + center * (opp_shrink_f);

		//dart other extremety: 2N -> 3N-1
		for (uint32 i = 0; i < count; ++i)
			local_vertices.push_back(local_vertices[i]*(opp_shrink_e) + local_vertices[(i+1)%count]*Scalar(shrink_e_));

		//phi2 mid-dart: 3N -> 4N-1
		for (uint32 i = 0; i < count; ++i)
			local_vertices.push_back((local_vertices[i]+local_vertices[(2*count+i+1)%count])/Scalar(2.0));

		for (uint32 i = 0; i < count; ++i)
		{
			darts_id_.push_back(local_darts[i]);
			const VEC3& P1 = local_vertices[i];
			darts_pos_.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
			const VEC3& P2 = local_vertices[2*count+i];
			darts_pos_.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});

			const VEC3& P3 = local_vertices[count+i];
			out_pos2.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
			const VEC3& P4 = local_vertices[3*count+i];
			out_pos2.push_back({float32(P4[0]), float32(P4[1]), float32(P4[2])});
		}
	});

	std::vector<Vec3f> darts_col;
	darts_col.resize(darts_pos_.size());
	for (auto& c: darts_col)
	{
		c[0] = dart_color_.x();
		c[1] = dart_color_.y();
		c[2] = dart_color_.z();
	}

	uint32 nbvec = std::uint32_t(darts_pos_.size());

	vbo_darts_->allocate(nbvec, 3);
	vbo_darts_->bind();
	vbo_darts_->copy_data(0, nbvec * 12, darts_pos_[0].data());
	vbo_darts_->release();

	vbo_color_darts_->allocate(nbvec, 3);
	vbo_color_darts_->bind();
	vbo_color_darts_->copy_data(0, nbvec*12, darts_col[0].data());
	vbo_color_darts_->release();

	vbo_relations_->allocate(nbvec, 3);
	vbo_relations_->bind();
	vbo_relations_->copy_data(0, nbvec * 12, out_pos2[0].data());
	vbo_relations_->release();
}

template <typename MAP, typename VERTEX_ATTR>
auto TopoDrawer::update(const MAP& m, const VERTEX_ATTR& position)
-> typename std::enable_if<MAP::DIMENSION == 3, void>::type
{
	static_assert(is_orbit_of<VERTEX_ATTR, MAP::Vertex::ORBIT>::value,"position must be a vertex attribute");

	using VEC3 = InsideTypeOf<VERTEX_ATTR>;
	using Vertex = typename MAP::Vertex;
	using Face = typename MAP::Face;
	using Volume = typename MAP::Volume;
	using Scalar = typename geometry::vector_traits<VEC3>::Scalar;

	Scalar opp_shrink_e = Scalar(1.0f - shrink_e_);
	Scalar opp_shrink_f = Scalar(1.0f - shrink_f_);
	Scalar opp_shrink_v = Scalar(1.0f - shrink_v_);

	darts_pos_.clear();
	darts_pos_.reserve(1024 * 1024);

	darts_id_.clear();
	darts_id_.reserve(1024 * 1024);

	std::vector<Vec3f> out_pos2;
	out_pos2.reserve(1024 * 1024);

	std::vector<Vec3f> out_pos3;
	out_pos3.reserve(1024 * 1024);

	std::vector<VEC3> local_vertices;
	local_vertices.reserve(256);

	std::vector<Dart> local_darts;
	local_darts.reserve(256);

	m.foreach_cell([&] (Volume v)
	{
		VEC3 center_vol = geometry::centroid(m, v, position);
		m.foreach_incident_face(v, [&] (Face f)
		{
			local_vertices.clear();
			local_darts.clear();
			VEC3 center;
			center.setZero();
			uint32 count = 0u;
			m.foreach_incident_vertex(f, [&] (Vertex v)
			{
				local_vertices.push_back(position[v]);
				local_darts.push_back(v.dart);
				center += position[v];
				count++;
			});
			center /= Scalar(count);

			// phi2 mid-edge: N -> 2N-1
			for (uint32 i = 0; i < count; ++i)
				local_vertices.push_back((local_vertices[i]+local_vertices[(i+1)%count])/Scalar(2.0));

			//phi3: 2N -> 3N-1
			for (uint32 i = 0; i < count; ++i)
				local_vertices.push_back(local_vertices[count+i]* shrink_f_ + center * (opp_shrink_f));

			// dart round point: 0 -> N-1
			for (uint32 i = 0; i < count; ++i)
				local_vertices[i] = local_vertices[i] * shrink_f_ + center * (opp_shrink_f);

			//dart other extremety: 3N -> 4N-1
			for (uint32 i = 0; i < count; ++i)
				local_vertices.push_back(local_vertices[i]*(opp_shrink_e) + local_vertices[(i+1)%count]*shrink_e_);

			//phi2/3 mid-dart: 4N -> 5N-1
			for (uint32 i = 0; i < count; ++i)
				local_vertices.push_back((local_vertices[i]+local_vertices[(2*count+i+1)%count])/Scalar(2.0));

			for (uint32 i = 0; i < count; ++i)
			{
				darts_id_.push_back(local_darts[i]);
				VEC3 P1 = (local_vertices[i] * shrink_v_) + (center_vol * opp_shrink_v);
				darts_pos_.push_back({float32(P1[0]), float32(P1[1]), float32(P1[2])});
				VEC3 P2 = (local_vertices[3*count+i] * shrink_v_) + (center_vol * opp_shrink_v);
				darts_pos_.push_back({float32(P2[0]), float32(P2[1]), float32(P2[2])});

				const VEC3 P3 = (local_vertices[count+i] * shrink_v_) + (center_vol * opp_shrink_v);
				out_pos2.push_back({float32(P3[0]), float32(P3[1]), float32(P3[2])});
				const VEC3 P4 = (local_vertices[4*count+i] * shrink_v_) + (center_vol * opp_shrink_v);
				out_pos2.push_back({float32(P4[0]), float32(P4[1]), float32(P4[2])});
				const VEC3& P5 = local_vertices[2*count+i];
				out_pos3.push_back({float32(P5[0]), float32(P5[1]), float32(P5[2])});
				out_pos3.push_back({float32(P4[0]), float32(P4[1]), float32(P4[2])});
			}
		});
	});

	std::vector<Vec3f> darts_col;
	darts_col.resize(darts_pos_.size());
	for (auto& c: darts_col)
	{
		c[0] = dart_color_.x();
		c[1] = dart_color_.y();
		c[2] = dart_color_.z();
	}

	uint32 nbvec = uint32(darts_pos_.size());
	vbo_darts_->allocate(nbvec, 3);
	vbo_darts_->bind();
	vbo_darts_->copy_data(0, nbvec * 12, darts_pos_[0].data());
	vbo_darts_->release();

	vbo_color_darts_->allocate(nbvec, 3);
	vbo_color_darts_->bind();
	vbo_color_darts_->copy_data(0, nbvec*12, darts_col[0].data());
	vbo_color_darts_->release();

	vbo_relations_->allocate(2 * nbvec, 3);
	vbo_relations_->bind();
	vbo_relations_->copy_data(0, nbvec * 12, out_pos2[0].data());
	vbo_relations_->copy_data(nbvec * 12, nbvec*12, out_pos3[0].data());
	vbo_relations_->release();
}



template <typename VERTEX_ATTR>
void TopoDrawer::update_colors(const VERTEX_ATTR& color)
{
	using VEC3 = InsideTypeOf<VERTEX_ATTR>;

	std::vector<Vec3f> darts_col;
	darts_col.reserve(2*darts_id_.size());
//	darts_col.clear();

	for (Dart d: darts_id_)
	{
		const VEC3& col = color[d];
		darts_col.push_back({float32(col[0]),float32(col[1]),float32(col[2])});
		darts_col.push_back({float32(col[0]),float32(col[1]),float32(col[2])});
	}

	uint32 nbvec = darts_col.size();
	vbo_color_darts_->allocate(nbvec, 3);
	vbo_color_darts_->bind();
	vbo_color_darts_->copy_data(0, nbvec*12, darts_col[0].data());
	vbo_color_darts_->release();

}

template <typename RGB>
void TopoDrawer::update_color(Dart d, const RGB& rgb)
{
	auto it = std::find(darts_id_.begin(), darts_id_.end(), d);
	if (it != darts_id_.end())
	{
		std::size_t x = it - darts_id_.begin();

		vbo_color_darts_->bind();
		float32 rgbf[6] = {float32(rgb[0]),float32(rgb[1]),float32(rgb[2]),
						  float32(rgb[0]),float32(rgb[1]),float32(rgb[2])};
		vbo_color_darts_->copy_data(uint32(x)*24u, 24u, rgbf);
		vbo_color_darts_->release();
	}
}

template <typename VEC3, typename VEC4>
Dart TopoDrawer::pick(const VEC3& xA, const VEC3& xB, const VEC4& plane, VEC3* xdp1, VEC3* xdp2)
{
//	using LVEC = geometry::Vec_T<Vec3f>;

	VEC3 xAB = xB-xA;
	Vec3f A(xA[0],xA[1],xA[2]);
	Vec3f AB(xAB[0],xAB[1],xAB[2]);

	float32 dmax = std::numeric_limits<float32>::max();
	float32 AB2 = AB.dot(AB);

	std::size_t isel = INVALID_INDEX;

	for(std::size_t i=0, nb_d = darts_id_.size(); i<nb_d; ++i)
	{
		const Vec3f& PP = darts_pos_[2*i];
		const Vec3f& QQ = darts_pos_[2*i+1];

		float32 prod1 = PP[0]*float32(plane[0]);
		prod1 += PP[1]*float32(plane[1]);
		prod1 += PP[2]*float32(plane[2]);
		prod1 += float32(plane[3]);


		float32 prod2 = QQ[0]*float32(plane[0]);
		prod2 += QQ[1]*float32(plane[1]);
		prod2 += QQ[2]*float32(plane[2]);
		prod2 += float32(plane[3]);

		if ((prod1<=0.0f)||(prod2<=0.0f))
		{
			float32 d2 = geometry::squared_distance_line_seg(A, AB, AB2, PP, QQ);
			if (d2 < dmax)
			{
				dmax = d2;
				isel = i;
			}
		}
	}

	if (isel != INVALID_INDEX)
	{
		if (xdp1 && xdp2)
		{
			Vec3f dp1 = darts_pos_[2*isel];
			Vec3f dp2 = darts_pos_[2*isel+1];
			*xdp1 = VEC3(dp1[0],dp1[1],dp1[2]);
			*xdp2 = VEC3(dp2[0],dp2[1],dp2[2]);
		}
		return darts_id_[isel];
	}

	return Dart(INVALID_INDEX);
}



template <typename VEC3, typename VEC4>
Dart TopoDrawer::pick(const VEC3& xA, const VEC3& xB, const VEC4& plane1, const VEC4& plane2, VEC3* xdp1, VEC3* xdp2)
{
//	using LVEC = geometry::Vec_T<Vec3f>;

	VEC3 xAB = xB-xA;
	Vec3f A(xA[0],xA[1],xA[2]);
	Vec3f AB(xAB[0],xAB[1],xAB[2]);

	float32 dmax = std::numeric_limits<float32>::max();
	float32 AB2 = AB.dot(AB);

	std::size_t isel = INVALID_INDEX;

	for(std::size_t i=0, nb_d = darts_id_.size(); i<nb_d; ++i)
	{
		const Vec3f& PP = darts_pos_[2*i];
		const Vec3f& QQ = darts_pos_[2*i+1];

		float32 prod1 = PP[0]*float32(plane1[0]);
		prod1 += PP[1]*float32(plane1[1]);
		prod1 += PP[2]*float32(plane1[2]);
		prod1 += float32(plane1[3]);

		float32 prod2 = QQ[0]*float32(plane1[0]);
		prod2 += QQ[1]*float32(plane1[1]);
		prod2 += QQ[2]*float32(plane1[2]);
		prod2 += float32(plane1[3]);

		float32 prod3 = PP[0]*float32(plane2[0]);
		prod3 += PP[1]*float32(plane2[1]);
		prod3 += PP[2]*float32(plane2[2]);
		prod3 += float32(plane2[3]);

		float32 prod4 = QQ[0]*float32(plane2[0]);
		prod4 += QQ[1]*float32(plane2[1]);
		prod4 += QQ[2]*float32(plane2[2]);
		prod4 += float32(plane2[3]);

		if (((prod1<=0.0f)||(prod2<=0.0f)) && ((prod3<=0.0f)||(prod4<=0.0f)))
		{
			float32 d2 = geometry::squared_distance_line_seg(A, AB, AB2, PP, QQ);
			if (d2 < dmax)
			{
				dmax = d2;
				isel = i;
			}
		}
	}

	if (isel != INVALID_INDEX)
	{
		if (xdp1 && xdp2)
		{
			Vec3f dp1 = darts_pos_[2*isel];
			Vec3f dp2 = darts_pos_[2*isel+1];
			*xdp1 = VEC3(dp1[0],dp1[1],dp1[2]);
			*xdp2 = VEC3(dp2[0],dp2[1],dp2[2]);
		}
		return darts_id_[isel];
	}

	return Dart(INVALID_INDEX);
}

template <typename VEC3, typename VEC4>
Dart TopoDrawer::pick(const VEC3& xA, const VEC3& xB, const VEC4& plane, float32 thickness, VEC3* xdp1, VEC3* xdp2)
{
	VEC4 p1 = plane;
	p1[3] -= thickness/2.0f;
	VEC4 p2 = -plane;
	p2[3] -= thickness/2.0f;
	return pick(xA,xB,p1,p2,xdp1,xdp2);
}


} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_TOPO_DRAWER_H_
