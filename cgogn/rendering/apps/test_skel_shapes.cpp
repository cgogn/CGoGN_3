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


using namespace ::cgogn;
using namespace ::cgogn::rendering;

class SkelShapeRender : public ui::ViewModule
{
public:
	SkelShapeRender(const ui::App& app) : ViewModule(app, "Test skel Shape"), app_(app),
		sphere_drawer_(nullptr),
		cone_drawer_(nullptr)
	{}

	~SkelShapeRender()
	{
	}

public:
	void init() override
	{
		app_.current_view()->set_scene_radius(1.5f);
		app_.current_view()->set_scene_center(GLVec3(0,0,0));
		sphere_drawer_ = SkelSphereDrawer::instance();
		cone_drawer_ = SkelConeDrawer::instance();
		sphere_drawer_->set_subdiv(40);
		cone_drawer_->set_subdiv(40);


		//GLVec4 C3 = {3.0f, 0.0f, 0.0f, 0.4f};
		//sphere_drawer_->set_spheres({C1, C2, C3});

		// std::vector<GLVec4> sphs;
		// static const int NB = 5;
		// sphs.reserve((2 * NB + 1) * (2 * NB + 1) * (2 * NB + 1));
		// for (int k = -NB; k < NB; ++k)
		// 	for (int j = -NB; j < NB; ++j)
		// 		for (int i = -NB; i < NB; ++i)
		// 			sphs.push_back({k / float(NB), j / float(NB), i / float(NB), // 0.2f});
		// 							(float32(rand()) / RAND_MAX) / (2*NB)});
		// sphere_drawer_->set_spheres(sphs);

		GLVec4 C1 = {0.0f, 0.0f, 0.0f, 0.3f};
		GLVec4 C2 = {0.4f, 2.0f, 0.2f, 1.1f};
		GLVec4 C3 = {3.0f, 0.0f, 0.0f, 0.4f};
		sphere_drawer_->set_spheres({C1, C2});

// 		cone_drawer_->compute_skel_cone(C1,C2,P1,P2,N1,N2);
// 		GLVec4 P1,P2;
// 		GLVec3 N1,N2;
// 		GLVec4 P1b,P2b;
// 		GLVec3 N1b,N2b;
// 		cone_drawer_->compute_skel_cone(C2,C3,P1b,P2b,N1b,N2b);
// //		cone_drawer_->set_cones({P1},{P2} ,{N1},{N2});

		// GLVec3 N1 = {1.0f,0.0f,0.0f};
		// GLVec3 N2 = {0.0f,0.0f,1.0f};
		// cone_drawer_->set_cones({C1},{C2} ,{N1},{N2});
 		GLVec4 P1,P2;
 		GLVec3 N1,N2;
		cone_drawer_->compute_skel_cone(C1,C2,P1,P2,N1,N2);
		cone_drawer_->set_cones({P1},{P2} ,{N1},{N2});

	}

	void draw(ui::View * view) override
	{
		const GLMat4& proj_matrix = view->projection_matrix();
		const GLMat4& view_matrix = view->modelview_matrix();
		sphere_drawer_->draw(proj_matrix, view_matrix);
		cone_drawer_->draw(proj_matrix, view_matrix);
	}

private:
	const ui::App& app_;
	SkelSphereDrawer* sphere_drawer_;
	SkelConeDrawer* cone_drawer_;
};

int main(int argc, char** argv)
{
	cgogn::ui::App app;
	app.set_window_title("Simple Shape Test");
	app.set_window_size(1000, 800);

	// declare a local module and link it
	SkelShapeRender sr(app);
	app.init_modules();
	app.current_view()->link_module(&sr);

	//no need gui here
	//app.show_gui(false);
	return app.launch();
}
