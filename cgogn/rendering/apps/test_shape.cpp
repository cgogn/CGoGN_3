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

#include <cgogn/rendering/shape_drawer.h>
using namespace ::cgogn;
using namespace ::cgogn::rendering;

class ShapeRender : public ui::ViewModule
{
public:
	ShapeRender(const ui::App & app)
		: ViewModule(app, "Test Shape"), app_(app), shape_(nullptr)
	{}

	~ShapeRender()	{}

public:
	void init() override
	{
		app_.current_view()->set_scene_radius(5);
		app_.current_view()->set_scene_center(GLVec3(0,0,0));
		shape_ = ShapeDrawer::instance();
		shape_->update_subdivision(24);
		shape_->color(ShapeDrawer::CYLINDER) = GLColor(1,0,0,1);
		shape_->color(ShapeDrawer::SPHERE) = GLColor(0, 1, 0, 1);
		shape_->color(ShapeDrawer::CONE) = GLColor(0, 0, 1, 1);
		shape_->color(ShapeDrawer::CUBE) = GLColor(1, 0, 1, 1);
	}

	void draw(ui::View * view) override
	{
		const GLMat4& proj_matrix = view->projection_matrix();
		const GLMat4& view_matrix = view->modelview_matrix();
		Eigen::Affine3f transfo = Eigen::Translation3f(GLVec3(4, 0, 0)) * Eigen::Scaling(2.0f);
		shape_->draw(ShapeDrawer::CYLINDER, proj_matrix, view_matrix * transfo.matrix());
		transfo = Eigen::Translation3f(GLVec3(-4, 0, 0)) * Eigen::Scaling(2.0f);
		shape_->draw(ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());
		transfo = Eigen::Translation3f(GLVec3(0, -4, 0)) * Eigen::Scaling(2.0f);
		shape_->draw(ShapeDrawer::CONE, proj_matrix, view_matrix * transfo.matrix());
		transfo = Eigen::Translation3f(GLVec3(0, 4, 0)) * Eigen::Scaling(GLVec3(1.5f,0.5f,2));
		shape_->draw(ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo.matrix());

		shape_->drawCylinder(proj_matrix, view_matrix, GLVec3(-3.75f, -4, -4.1f), GLVec3(4.1f, 4.3f, 3.95f), 0.5f);

		shape_->drawSphere(proj_matrix, view_matrix, GLVec3(-3.75f, -4, -4.1f), 0.4f);
		shape_->drawSphere(proj_matrix, view_matrix, GLVec3(4.1f, 4.3f, 3.95f), 0.4f);

		shape_->drawCone(proj_matrix, view_matrix, GLVec3(0,0,2), GLVec3(-2,1,3.5), 0.75f);
		std::cout << "DRAW" << std::endl;
	}

private:
	const ui::App& app_;
	ShapeDrawer* shape_;
};

int main(int argc, char** argv)
{
	cgogn::ui::App app;
	app.set_window_title("Simple Shape Test");
	app.set_window_size(1000, 800);

	// declare a local module and link it
	ShapeRender sr(app);
	app.init_modules();
	app.current_view()->link_module(&sr);
	//no need gui here
	//app.show_gui(false);
	return app.launch();
}
