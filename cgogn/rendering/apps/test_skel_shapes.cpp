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

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/rendering/skelshape.h>
#include <cgogn/modeling/skeleton_sampling.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

using namespace ::cgogn;
using namespace ::cgogn::rendering;

class MySkelRender : public ui::ViewModule
{
public:
	MySkelRender(const ui::App& app)
		: ViewModule(app, "Test skel Shape"), app_(app)
	{
	}

	~MySkelRender()
	{
	}

	void init() override
	{
		app_.current_view()->set_scene_radius(4.0f);
		app_.current_view()->set_scene_center(GLVec3(0, 0, 0));

		skel_drawer_.set_color({1, 0, 1, 1});
		skel_drawer_.set_subdiv(40);

		GLVec4 C00 = {-2.5f, -2.0f, 0.9f, 0.3f};
		GLVec4 C0 = {-2.0f, 0.0f, 0.4f, 0.2f};

		GLVec4 C1 = {-0.5f, 0.0f, 0.0f, 0.4f};
		GLVec4 C2 = {0.4f, 1.0f, 0.2f, 0.7f};
		GLVec4 C3 = {1.0f, -1.0f, -0.5f, 0.2f};


		skel_drawer_.add_vertex(C00.topRows<3>(), C00[3]); // just to test version with radius parameter
		skel_drawer_.add_vertex(GLVec3{-2.0f, 0.0f, 0.4f}, 0.2f);
		skel_drawer_.add_vertex(C1);
		skel_drawer_.add_vertex(C2);
		skel_drawer_.add_vertex(C3);
		skel_drawer_.add_edge(C00, C0);
		skel_drawer_.add_edge(C0, C1);
		skel_drawer_.add_edge(C1, C2);
		skel_drawer_.add_edge(C2, C3);
		skel_drawer_.add_edge(C3, C1);
		skel_drawer_.add_triangle(C1, C2, C3);
		skel_drawer_.update();

		
		skel_sampler_.add_vertex(C00.topRows<3>(), C00[3]); // just to test version with radius parameter
		skel_sampler_.add_vertex(GLVec3{-2.0f, 0.0f, 0.4f}, 0.2f);
		skel_sampler_.add_vertex(C1);
		skel_sampler_.add_vertex(C2);
		skel_sampler_.add_vertex(C3);
		skel_sampler_.add_edge(C00, C0);
		skel_sampler_.add_edge(C0, C1);
		skel_sampler_.add_edge(C1, C2);
		skel_sampler_.add_edge(C2, C3);
		skel_sampler_.add_edge(C3, C1);
		skel_sampler_.add_triangle(C1, C2, C3);

		
		 // for (int i=0; i<100; i++)
		//{
		//	float alpha = 0.34f*i;
		//	float h = -2.0f+0.04f*i;
		//	float r = (i%2==1)?0.3f:0.2f;
		//	float rr = (i%2==0)?0.3f:0.2f;
		// skel_sampler_.add_edge(GLVec4(1.0f * std::cos(alpha), 1.0f * std::sin(alpha), h, r),
		//					   GLVec4(3.5f * std::cos(alpha), 3.5f * std::sin(alpha), h + 0.5, rr));
		// skel_sampler_.add_vertex(GLVec4(1.0f * std::cos(alpha), 1.0f * std::sin(alpha), h, r));
		// skel_sampler_.add_vertex(GLVec4(3.5f * std::cos(alpha), 3.5f * std::sin(alpha), h + 0.5, rr));

		//	skel_drawer_.add_edge(GLVec4(1.0f * std::cos(alpha), 1.0f * std::sin(alpha), h, r),
		//					   GLVec4(3.5f * std::cos(alpha), 3.5f * std::sin(alpha), h + 0.5, rr));
		// skel_drawer_.add_vertex(GLVec4(1.0f * std::cos(alpha), 1.0f * std::sin(alpha), h, r));
		// skel_drawer_.add_vertex(GLVec4(3.5f * std::cos(alpha), 3.5f * std::sin(alpha), h + 0.5, rr));
		//}
		//skel_drawer_.update();

		
				
		//	skel_sampler_.add_vertex(GLVec4(0, 0, 0, 1));


		GLVec3 bbw = skel_sampler_.BBwidth();
		float step = std::min(std::min(bbw.x(),bbw.y()),bbw.z())/200;
		skel_sampler_.sample(step);	

		// for(auto p:skel_sampler_.samples())
		//  	std::cout << p.norm()<< std::endl;

		param_point_sprite_ = ShaderPointSprite::generate_param();
		param_point_sprite_->color_ = rendering::GLColor(1, 1, 0, 1);
		param_point_sprite_->point_size_ = 0.0025;
		param_point_sprite_->set_vbos({&vbo_samples_});

		// std::vector<GLVec3> pp;
		// for (float z=-4.0f;z<=4.0f; z += 0.2)
		// 	for (float y=-4.0f;y<=4.0f; y += 0.2)
		// 		for (float x=-4.0f;x<=4.0f; x += 0.2)
		// 		pp.emplace_back(x,y,z);

		update_vbo(skel_sampler_.samples(), &vbo_samples_);
//		update_vbo(pp, &vbo_samples_);

		//
		//   for (const auto& p : skel_sampler_.samples())
		//   	std::cout << p.transpose()<< " => "<<p.norm()<<std::endl;		 
	}

	void draw(ui::View* view) override
	{
		const GLMat4& proj_matrix = view->projection_matrix();
		const GLMat4& view_matrix = view->modelview_matrix();
		skel_drawer_.draw(proj_matrix, view_matrix);

		if (param_point_sprite_->attributes_initialized())
		{
			std::cout << "DR "<< skel_sampler_.samples().size()<<std::endl;
			param_point_sprite_->bind(proj_matrix, view_matrix);
			glDrawArrays(GL_POINTS,0,vbo_samples_.size());//skel_sampler_.samples().size());
			param_point_sprite_->release();
		}						
						
	}

private:
	const ui::App& app_;
	SkelShapeDrawer skel_drawer_;

	cgogn::modeling::SkeletonSampler<GLVec4, GLVec3, float> skel_sampler_;
	VBO vbo_samples_;
	std::unique_ptr<ShaderPointSprite::Param> param_point_sprite_;
public:
};


int main(int argc, char** argv)
{
	cgogn::ui::App app;
	app.set_window_title("Simple Shape Test");
	app.set_window_size(1000, 800);

	// declare a local module and link it
	MySkelRender sr(app);
	app.init_modules();
	app.current_view()->link_module(&sr);

	//no need gui here
	//app.show_gui(false);
	return app.launch();
}
