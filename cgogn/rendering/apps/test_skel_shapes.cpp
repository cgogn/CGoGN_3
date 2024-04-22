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

		// 2 first in double others in float to check on the fly conversion
		GLVec4d C00 = {-2.5, -2.0, 0.9, 0.3};
		GLVec4d C0 = {-2.0, 0.0, 0.4, 0.2};

		GLVec4 C1 = {-0.5f, 0.0f, 0.0f, 0.4f};
		GLVec4 C2 = {0.4f, 1.0f, 0.2f, 0.7f};
		GLVec4 C3 = {1.0f, -1.0f, -0.5f, 0.2f};
		skel_drawer_.add_vertex(C00.topRows<3>(), C00[3]); // just to test version with radius parameter
		skel_drawer_.add_vertex(GLVec3d{-2.0, 0.0, 0.4}, 0.2);
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
	}

	void draw(ui::View* view) override
	{
		const GLMat4& proj_matrix = view->projection_matrix();
		const GLMat4& view_matrix = view->modelview_matrix();
		skel_drawer_.draw(proj_matrix, view_matrix);
	}

private:
	const ui::App& app_;
	SkelShapeDrawer skel_drawer_;

	SkeletonSampler<GLVec4, GLVec3, float> sksampler_;
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
