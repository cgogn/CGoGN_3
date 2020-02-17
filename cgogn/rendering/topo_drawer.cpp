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

#include <cgogn/rendering/topo_drawer.h>

namespace cgogn
{

namespace rendering
{

TopoDrawer::TopoDrawer()
	: dart_color_(1, 1, 1, 1), phi2_color_(1, 0, 0, 1), phi3_color_(1, 1, 0, 1), shrink_v_(0.8f), shrink_f_(0.85f),
	  shrink_e_(0.95f)
{
	vbo_darts_ = std::make_unique<VBO>(3);
	vbo_relations_ = std::make_unique<VBO>(3);
	vbo_color_darts_ = std::make_unique<VBO>(3);
}

TopoDrawer::~TopoDrawer()
{
}

TopoDrawer::Renderer::Renderer(TopoDrawer* tr) : topo_drawer_data_(tr)
{
	param_bl_ = ShaderBoldLineColor::generate_param();
	param_bl_->set_vbos({tr->vbo_darts_.get(), tr->vbo_color_darts_.get()});

	param_bl2_ = ShaderBoldLine::generate_param();
	param_bl2_->set_vbos({tr->vbo_relations_.get()});
	param_bl2_->color_ = tr->phi2_color_;

	param_rp_ = ShaderRoundPointColor::generate_param();

	param_rp_->bind_vao();
	tr->vbo_darts_->associate(1, 2, 0);
	tr->vbo_color_darts_->associate(2, 2, 0);
	param_rp_->release_vao();
}

TopoDrawer::Renderer::~Renderer()
{
}

void TopoDrawer::Renderer::draw(const GLMat4& projection, const GLMat4& modelview, bool with_blending)
{
	float32 lw = 2.0f;
	if (with_blending)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		lw = 3.0f;
	}

	param_bl_->width_ = lw;
	param_bl2_->width_ = lw;
	param_rp_->size_ = 2.0f * lw;

	param_bl_->bind(projection, modelview);
	
	glDrawArrays(GL_LINES, 0, topo_drawer_data_->vbo_darts_->size());
	
	param_bl_->release();
	

	param_rp_->bind(projection, modelview);
	

	glDrawArrays(GL_POINTS, 0, topo_drawer_data_->vbo_darts_->size() / 2);
	param_rp_->release();

	param_bl2_->color_ = topo_drawer_data_->phi2_color_;
	param_bl2_->bind(projection, modelview);
	glDrawArrays(GL_LINES, 0, topo_drawer_data_->vbo_darts_->size());
	param_bl2_->release();

	if (topo_drawer_data_->vbo_relations_->size() > topo_drawer_data_->vbo_darts_->size())
	{
		param_bl2_->color_ = topo_drawer_data_->phi3_color_;
		param_bl2_->bind(projection, modelview);
		glDrawArrays(GL_LINES, topo_drawer_data_->vbo_darts_->size(), topo_drawer_data_->vbo_darts_->size());
		param_bl2_->release();
	}

	if (with_blending)
		glDisable(GL_BLEND);
}

void TopoDrawer::Renderer::set_clipping_plane(const GLVec4& p)
{
	param_bl_->plane_clip_ = p;
	param_bl2_->plane_clip_ = p;
	param_rp_->plane_clip_ = p;
}

void TopoDrawer::Renderer::set_clipping_plane2(const GLVec4& p)
{
	param_bl_->plane_clip2_ = p;
	param_bl2_->plane_clip2_ = p;
	param_rp_->plane_clip2_ = p;
}

void TopoDrawer::Renderer::set_thick_clipping_plane(const GLVec4& p, float32 th)
{
	GLVec4 p1 = p;
	p1[3] -= th / 2.0f;
	set_clipping_plane(p1);

	GLVec4 p2 = -p;
	p2[3] -= th / 2.0f;
	set_clipping_plane2(p2);
}

void TopoDrawer::update_color(Dart d, const GLColor& rgb)
{
	auto it = std::find(darts_id_.begin(), darts_id_.end(), d);
	if (it != darts_id_.end())
	{
		std::size_t x = it - darts_id_.begin();

		vbo_color_darts_->bind();
		float32 rgbf[6] = {rgb.x(), rgb.y(), rgb.z(), rgb.x(), rgb.y(), rgb.z()};
		vbo_color_darts_->copy_data(uint32(x) * 24u, 24u, rgbf);
		vbo_color_darts_->release();
	}
}


void TopoDrawer::update_color(Dart d, const Vec3& rgb)
{
	auto it = std::find(darts_id_.begin(), darts_id_.end(), d);
	if (it != darts_id_.end())
	{
		std::size_t x = it - darts_id_.begin();

		vbo_color_darts_->bind();
		float32 rgbf[6] = {float32(rgb[0]), float32(rgb[1]), float32(rgb[2]),
						   float32(rgb[0]), float32(rgb[1]), float32(rgb[2])};
		vbo_color_darts_->copy_data(uint32(x) * 24u, 24u, rgbf);
		vbo_color_darts_->release();
	}
}

Dart TopoDrawer::pick(const Vec3& A, const Vec3& B, const Vec4& plane, Vec3* dp1, Vec3* dp2)
{
	Vec3 AB = B - A;

	Scalar dmax = std::numeric_limits<Scalar>::max();
	Scalar AB2 = AB.dot(AB);

	std::size_t isel = INVALID_INDEX;

	for (std::size_t i = 0, nb_d = darts_id_.size(); i < nb_d; ++i)
	{
		const Vec3& PP = darts_pos_[2 * i].cast<Scalar>();
		const Vec3& QQ = darts_pos_[2 * i + 1].cast<Scalar>();

		Scalar prod1 = PP[0] * Scalar(plane[0]);
		prod1 += PP[1] * Scalar(plane[1]);
		prod1 += PP[2] * Scalar(plane[2]);
		prod1 += Scalar(plane[3]);

		Scalar prod2 = QQ[0] * Scalar(plane[0]);
		prod2 += QQ[1] * Scalar(plane[1]);
		prod2 += QQ[2] * Scalar(plane[2]);
		prod2 += Scalar(plane[3]);

		if ((prod1 <= 0) || (prod2 <= 0))
		{
			Scalar d2 = geometry::squared_distance_line_seg(A, AB, AB2, PP, QQ);
			if (d2 < dmax)
			{
				dmax = d2;
				isel = i;
			}
		}
	}

	if (isel != INVALID_INDEX)
	{
		if (dp1 && dp2)
		{
			Vec3f fdp1 = darts_pos_[2 * isel];
			Vec3f fdp2 = darts_pos_[2 * isel + 1];
			*dp1 = Vec3(fdp1[0], fdp1[1], fdp1[2]);
			*dp2 = Vec3(fdp2[0], fdp2[1], fdp2[2]);
		}
		return darts_id_[isel];
	}

	return Dart(INVALID_INDEX);
}

Dart TopoDrawer::pick(const Vec3& A, const Vec3& B, const Vec4& plane1, const Vec4& plane2, Vec3* dp1, Vec3* dp2)
{
	//	using LVEC = geometry::Vec_T<Vec3f>;

	Vec3 AB = B - A;

	Scalar dmax = std::numeric_limits<Scalar>::max();
	Scalar AB2 = AB.dot(AB);

	std::size_t isel = INVALID_INDEX;

	for (std::size_t i = 0, nb_d = darts_id_.size(); i < nb_d; ++i)
	{
		const Vec3& PP = darts_pos_[2 * i].cast<Scalar>();
		const Vec3& QQ = darts_pos_[2 * i + 1].cast<Scalar>();

		Scalar prod1 = PP[0] * Scalar(plane1[0]);
		prod1 += PP[1] * Scalar(plane1[1]);
		prod1 += PP[2] * Scalar(plane1[2]);
		prod1 += Scalar(plane1[3]);

		Scalar prod2 = QQ[0] * Scalar(plane1[0]);
		prod2 += QQ[1] * Scalar(plane1[1]);
		prod2 += QQ[2] * Scalar(plane1[2]);
		prod2 += Scalar(plane1[3]);

		Scalar prod3 = PP[0] * Scalar(plane2[0]);
		prod3 += PP[1] * Scalar(plane2[1]);
		prod3 += PP[2] * Scalar(plane2[2]);
		prod3 += Scalar(plane2[3]);

		Scalar prod4 = QQ[0] * Scalar(plane2[0]);
		prod4 += QQ[1] * Scalar(plane2[1]);
		prod4 += QQ[2] * Scalar(plane2[2]);
		prod4 += Scalar(plane2[3]);

		if (((prod1 <= 0) || (prod2 <= 0)) && ((prod3 <= 0) || (prod4 <= 0)))
		{
			Scalar d2 = geometry::squared_distance_line_seg(A, AB, AB2, PP, QQ);
			if (d2 < dmax)
			{
				dmax = d2;
				isel = i;
			}
		}
	}

	if (isel != INVALID_INDEX)
	{
		if (dp1 && dp2)
		{
			Vec3f fdp1 = darts_pos_[2 * isel];
			Vec3f fdp2 = darts_pos_[2 * isel + 1];
			*dp1 = Vec3(fdp1[0], fdp1[1], fdp1[2]);
			*dp2 = Vec3(fdp2[0], fdp2[1], fdp2[2]);
		}
		return darts_id_[isel];
	}

	return Dart(INVALID_INDEX);
}


Dart TopoDrawer::pick(const Vec3& A, const Vec3& B, const Vec4& plane, float32 thickness, Vec3* dp1, Vec3* dp2)
{
	Vec4 p1 = plane;
	p1[3] -= thickness / 2.0f;
	Vec4 p2 = -plane;
	p2[3] -= thickness / 2.0f;
	return pick(A, B, p1, p2, dp1, dp2);
}










} // namespace rendering

} // namespace cgogn
