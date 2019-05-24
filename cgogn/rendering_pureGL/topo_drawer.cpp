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


#include <cgogn/rendering_pureGL/topo_drawer.h>


namespace cgogn
{

namespace rendering_pgl
{

TopoDrawer::TopoDrawer() :
	dart_color_(1,1,1,1),
	phi2_color_(1,0,0,1),
	phi3_color_(1,1,0,1),
	shrink_v_(0.6f),
	shrink_f_(0.85f),
	shrink_e_(0.95f)
{
	vbo_darts_ = cgogn::make_unique<VBO>(3);
	vbo_relations_ = cgogn::make_unique<VBO>(3);
	vbo_color_darts_ =  cgogn::make_unique<VBO>(3);
}

TopoDrawer::~TopoDrawer()
{}

TopoDrawer::Renderer::Renderer(TopoDrawer* tr) :
	topo_drawer_data_(tr)
{
	param_bl_ = ShaderBoldLineColor::generate_param();
	param_bl_->set_vbos(tr->vbo_darts_.get(), tr->vbo_color_darts_.get());

	param_bl2_ = ShaderBoldLine::generate_param();
	param_bl2_->set_vbos(tr->vbo_relations_.get());
	param_bl2_->color_= tr->phi2_color_;

	param_rp_ = ShaderRoundPointColor::generate_param();

	param_rp_->bind_vao();
	tr->vbo_darts_->associate(1,2,0);
	tr->vbo_color_darts_->associate(2,2,0);
	param_rp_->release_vao();
}

TopoDrawer::Renderer::~Renderer()
{}

void TopoDrawer::Renderer::draw(const GLMat4& projection, const GLMat4& modelview, bool with_blending)
{
	float32 lw = 2.0f;
	if(with_blending)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		lw = 3.0f;
	}

	param_bl_->width_ = lw;
	param_bl2_->width_ = lw;
	param_rp_->size_ = 2.0f * lw;

	param_rp_->bind(projection, modelview);
	glDrawArrays(GL_POINTS, 0, topo_drawer_data_->vbo_darts_->size()/2);
	param_rp_->release();

	param_bl_->bind(projection, modelview);
	glDrawArrays(GL_LINES, 0, topo_drawer_data_->vbo_darts_->size());
	param_bl_->release();

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
	p1[3] -= th/2.0f;
	set_clipping_plane(p1);

	GLVec4 p2 = -p;
	p2[3] -= th/2.0f;
	set_clipping_plane2(p2);
}


void TopoDrawer::update_color(Dart d, const GLColor& rgb)
{
	auto it = std::find(darts_id_.begin(), darts_id_.end(), d);
	if (it != darts_id_.end())
	{
		std::size_t x = it - darts_id_.begin();

		vbo_color_darts_->bind();
		float32 rgbf[6] = {rgb.x(),rgb.y(),rgb.z(),rgb.x(),rgb.y(),rgb.z()};
		vbo_color_darts_->copy_data(uint32(x)*24u, 24u, rgbf);
		vbo_color_darts_->release();
	}
}

} // namespace rendering

} // namespace cgogn
